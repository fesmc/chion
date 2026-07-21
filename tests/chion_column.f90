program chion_column
    ! WP15 -- single-column driver, `chion_column.x`.
    !
    !     libchion/bin/chion_column.x par/chion_column.nml
    !
    ! Exactly one command-line argument, the parameter file; output written to
    ! the current working directory. That is the runme contract
    ! (.runme/info.json: par_path_as_argument, and the run directory is CWD).
    !
    ! The forcing is SYNTHETIC and fully specified by the &ctrl group of the
    ! parameter file, so that a Chion.jl run can be driven from the same
    ! numbers and the two output files diffed variable by variable (WP16). The
    ! &ctrl group belongs to this driver, not to the chion library: chion never
    ! reads it and it is not validated against input/chion_defaults.nml.
    !
    ! ---------------------------------------------------------------------
    ! The synthetic annual cycle, in full
    ! ---------------------------------------------------------------------
    ! With d = day of year (1-based, fractional) and Y = YEAR_LENGTH = 365 d:
    !
    !   air_temperature = t2m_mean + t2m_amplitude*cos(2*pi*(d - t2m_phase)/Y)
    !   shortwave_down  = max(sw_mean + sw_amplitude*cos(2*pi*(d - t2m_phase)/Y), 0)
    !   precipitation   = pr_mean                        (constant)
    !   snowfall_rate   = precipitation  if air_temperature <  t_snow_max
    !   rainfall_rate   = precipitation  if air_temperature >= t_snow_max
    !   wind_speed      = wind_speed                     (constant)
    !   latitude_deg    = latitude
    !   surface_height  = surface_height
    !   H_ice           = H_ice                          (ITM only)
    !   solar_longitude_deg = 360*(d - 1)/Y
    !
    ! Shortwave shares the temperature phase, which is deliberate: it is a
    ! driver simplification, not a solar-geometry calculation, and keeping one
    ! phase parameter means one number to match in Chion.jl.
    !
    ! YEAR_LENGTH is 365 days, hard-coded here rather than added to &ctrl. It
    ! is fixed by par/chion_column.nml's own comment ("time_end = 3650.0 !
    ! 10 years"), and a second, independently settable year length would be a
    ! way for the two to disagree silently.
    !
    ! ---------------------------------------------------------------------
    ! PDDs, and why it is computed once
    ! ---------------------------------------------------------------------
    ! ITM's calc_albedo_surface interpolates the critical snow depth between a
    ! "desert" and a "forest" end member using the ANNUAL positive-degree-day
    ! total. smbpal recomputes it once per year and holds it fixed for every
    ! step of that year (smbpal.f90:378-385); a host that writes a per-step
    ! degree-day increment there instead gets the desert branch always and a
    ! systematically different albedo (chion_defs, chion_forcing_class).
    !
    ! The forcing here is periodic, so the annual total is the same every year
    ! and is summed once, analytically over the 365 daily values, before the
    ! time loop.
    !
    ! ---------------------------------------------------------------------
    ! Output
    ! ---------------------------------------------------------------------
    !   <file_out>             every dt_out days, the model's full variable
    !                          table including the layer state (WP14)
    !   chion_column_restart.nc  written once at the end of the run
    !
    ! Both in the CWD. The restart file name is fixed rather than being a
    ! &ctrl parameter, so that the parameter file did not have to be extended.

    use chion
    use chion_io, only : chion_write_init, chion_write_step, chion_restart_write
    use nml,      only : nml_read

    implicit none

    ! --- &ctrl ------------------------------------------------------------
    integer  :: ncol
    real(wp) :: time_init, time_end, dt, dt_out
    real(wp) :: latitude, surface_height, H_ice
    real(wp) :: t2m_mean, t2m_amplitude, t2m_phase
    real(wp) :: pr_mean, t_snow_max
    real(wp) :: sw_mean, sw_amplitude, wind_speed
    character(len=512) :: file_out

    ! --- Driver state -----------------------------------------------------
    type(chion_class) :: chn

    character(len=512) :: path_par
    character(len=*), parameter :: file_restart = "chion_column_restart.nc"

    real(wp), parameter :: YEAR_LENGTH = 365.0_wp
    real(wp), parameter :: PI          = 3.14159265358979323846_wp

    integer  :: k, nstep, n_out, n_year
    real(wp) :: time, doy, t_air, sw_dn, precip, pdds_annual
    real(wp) :: time_next_out

    real(wp), allocatable :: smb(:)
    real(wp), allocatable :: thickness(:), wet_mass(:), bulk_density(:), liquid_water(:)

    real(wp)     :: depth_peak, depth_now
    integer      :: n_lay_peak
    real(wp_acc) :: smb_year

    ! =====================================================================
    ! Setup
    ! =====================================================================

    call chion_load_command_line_args(path_par)
    call chion_check_file(path_par)

    ! The &ctrl group is read WITHOUT a defaults file: it is this program's
    ! own schema, so the group must be complete, and a missing entry is
    ! nml's ERROR_NO_PARAM rather than a silent zero.
    call nml_read(path_par,"ctrl","ncol",           ncol)
    call nml_read(path_par,"ctrl","time_init",      time_init)
    call nml_read(path_par,"ctrl","time_end",       time_end)
    call nml_read(path_par,"ctrl","dt",             dt)
    call nml_read(path_par,"ctrl","dt_out",         dt_out)
    call nml_read(path_par,"ctrl","latitude",       latitude)
    call nml_read(path_par,"ctrl","surface_height", surface_height)
    call nml_read(path_par,"ctrl","H_ice",          H_ice)
    call nml_read(path_par,"ctrl","t2m_mean",       t2m_mean)
    call nml_read(path_par,"ctrl","t2m_amplitude",  t2m_amplitude)
    call nml_read(path_par,"ctrl","t2m_phase",      t2m_phase)
    call nml_read(path_par,"ctrl","pr_mean",        pr_mean)
    call nml_read(path_par,"ctrl","t_snow_max",     t_snow_max)
    call nml_read(path_par,"ctrl","sw_mean",        sw_mean)
    call nml_read(path_par,"ctrl","sw_amplitude",   sw_amplitude)
    call nml_read(path_par,"ctrl","wind_speed",     wind_speed)
    call nml_read(path_par,"ctrl","file_out",       file_out)

    if (ncol .lt. 1) then
        write(io_unit_err,*) "chion_column:: Error: ctrl:ncol must be positive."
        write(io_unit_err,*) "ncol = ", ncol
        stop "Program stopped."
    end if

    if (dt .le. 0.0_wp .or. dt_out .le. 0.0_wp .or. time_end .le. time_init) then
        write(io_unit_err,*) "chion_column:: Error: ctrl:dt, dt_out must be positive &
                             &and time_end > time_init."
        write(io_unit_err,*) "dt, dt_out, time_init, time_end = ", dt, dt_out, time_init, time_end
        stop "Program stopped."
    end if

    call chion_init(chn,path_par,ncol)
    call chion_init_state(chn)

    allocate(smb(ncol))
    allocate(thickness(ncol))
    allocate(wet_mass(ncol))
    allocate(bulk_density(ncol))
    allocate(liquid_water(ncol))

    ! --- Time-invariant forcing -------------------------------------------

    chn%forc%latitude_deg   = latitude
    chn%forc%surface_height = surface_height
    chn%forc%wind_speed     = wind_speed
    chn%forc%H_ice          = H_ice

    ! Annual positive-degree-day total of the prescribed cycle. Same every
    ! year, so summed once (see the header).
    pdds_annual = 0.0_wp
    do k = 1, nint(YEAR_LENGTH)
        t_air = column_air_temperature(real(k,wp))
        pdds_annual = pdds_annual + max(t_air - chn%c%T0,0.0_wp)
    end do
    chn%forc%PDDs = pdds_annual

    ! =====================================================================
    ! Run
    ! =====================================================================

    nstep = nint((time_end - time_init)/dt)

    write(*,"(a)")        "== chion_column ==========================================="
    write(*,"(a,a)")      " par file   : ", trim(path_par)
    write(*,"(a,a)")      " output     : ", trim(file_out)
    write(*,"(a,a)")      " restart    : ", trim(file_restart)
    write(*,"(a,i0,a,f10.2,a,f10.2)") " steps      : ", nstep, "   from ", time_init, " to ", time_end
    write(*,"(a,f10.4,a,f10.4)")      " dt, dt_out : ", dt, "  ", dt_out
    write(*,"(a,f10.3,a)")            " PDDs       : ", pdds_annual, " [K d yr-1]"
    write(*,"(a)")        "==========================================================="

    call chion_write_init(chn,file_out,time_init,"days")
    call chion_write_step(chn,file_out,time_init)

    time          = time_init
    time_next_out = time_init + dt_out
    n_out         = 1

    depth_peak = 0.0_wp
    n_lay_peak = 0
    smb_year   = 0.0_wp_acc
    n_year     = 0

    do k = 1, nstep

        time = time_init + real(k,wp)*dt

        ! Forcing is evaluated at the MIDPOINT of the step, which is the
        ! centred choice for a step that integrates over [time-dt, time].
        doy   = column_day_of_year(time - 0.5_wp*dt)
        t_air = column_air_temperature(doy)
        sw_dn = max(sw_mean + sw_amplitude*cos(2.0_wp*PI*(doy - t2m_phase)/YEAR_LENGTH),0.0_wp)

        precip = pr_mean

        chn%forc%air_temperature = t_air
        chn%forc%shortwave_down  = sw_dn

        if (t_air .lt. t_snow_max) then
            chn%forc%snowfall_rate = precip
            chn%forc%rainfall_rate = 0.0_wp
        else
            chn%forc%snowfall_rate = 0.0_wp
            chn%forc%rainfall_rate = precip
        end if

        chn%forc%day_of_year         = doy
        chn%forc%solar_longitude_deg = 360.0_wp*(doy - 1.0_wp)/YEAR_LENGTH

        call chion_update(chn,dt)

        ! --- Running diagnostics --------------------------------------

        call chion_get_smb(chn,smb)

        depth_now = column_snow_depth(chn,thickness,wet_mass,bulk_density,liquid_water)
        if (depth_now .gt. depth_peak) depth_peak = depth_now

        if (trim(chn%par%model) .eq. "bessi") then
            if (maxval(chn%bsi%now%n_lay) .gt. n_lay_peak) n_lay_peak = maxval(chn%bsi%now%n_lay)
        end if

        ! Ice-facing mass flux over the FINAL year, accumulated in wp_acc for
        ! the same reason the model accumulators are (docs/porting_notes.md D1).
        if (time .gt. time_end - YEAR_LENGTH) then
            smb_year = smb_year + real(smb(1),wp_acc)*real(dt,wp_acc) &
                                  *real(chn%c%seconds_per_day,wp_acc)
            n_year   = n_year + 1
        end if

        ! --- Output ---------------------------------------------------

        if (time .ge. time_next_out - 0.5_wp*dt) then
            call chion_write_step(chn,file_out,time)
            n_out         = n_out + 1
            time_next_out = time_next_out + dt_out
        end if

    end do

    ! =====================================================================
    ! Summary
    ! =====================================================================

    call chion_restart_write(chn,file_restart,time)

    write(*,"(a)") ""
    write(*,"(a)") "== chion_column summary ==================================="
    write(*,"(a,a)")     " model              : ", trim(chn%par%model)
    write(*,"(a,i0)")    " output records     : ", n_out
    write(*,"(a,f14.4)") " final time    [d]  : ", time
    write(*,"(a,f14.4,a,i0,a)") " annual SMB [kg m-2 yr-1] : ", real(smb_year,wp), &
                                "   (final year, ", n_year, " steps)"
    write(*,"(a,f14.4)") " peak snow depth    : ", depth_peak

    select case(trim(chn%par%model))

        case("bessi")
            write(*,"(a)")       "   [peak snow depth is thickness, m]"
            write(*,"(a,i0)")    " peak layer count   : ", n_lay_peak
            write(*,"(a,i0)")    " final layer count  : ", chn%bsi%now%n_lay(1)
            write(*,"(a,f14.4)") " final thickness [m]: ", thickness(1)
            write(*,"(a,f14.4)") " final bulk rho     : ", bulk_density(1)
            write(*,"(a,f14.4)") " cum mass_base      : ", real(chn%bsi%now%mass_base(1),wp)
            write(*,"(a,f14.4)") " cum smb_ice        : ", real(chn%bsi%now%smb_ice(1),wp)
            write(*,"(a,f14.4)") " cum runoff         : ", real(chn%bsi%now%runoff(1),wp)
            write(*,"(a,f14.4)") " cum melt           : ", real(chn%bsi%now%melt(1),wp)
            write(*,"(a,f14.4)") " cum refreezing     : ", real(chn%bsi%now%refreezing(1),wp)
            write(*,"(a,f14.4)") " cum sublimation    : ", real(chn%bsi%now%sublimation(1),wp)
            write(*,"(a,f14.4)") " final albedo       : ", chn%bsi%now%albedo(1)
            write(*,"(a,f14.4)") " final Tsrf     [K] : ", chn%bsi%now%t_srf(1)

        case("pdd")
            write(*,"(a)")       "   [peak snow depth is snowpack_swe, kg m-2]"
            write(*,"(a,f14.4)") " final swe [kg m-2] : ", chn%pdd%now%snowpack_swe(1)
            write(*,"(a,f14.4)") " cum smb_ice        : ", real(chn%pdd%now%smb_ice(1),wp)
            write(*,"(a,f14.4)") " cum runoff         : ", real(chn%pdd%now%runoff(1),wp)
            write(*,"(a,f14.4)") " cum pdd_sum  [K d] : ", real(chn%pdd%now%pdd_sum(1),wp)

        case("itm")
            write(*,"(a)")       "   [peak snow depth is H_snow, mm w.e.]"
            write(*,"(a,f14.4)") " final H_snow       : ", chn%itm%now%H_snow(1)
            write(*,"(a,f14.4)") " final albedo       : ", chn%itm%now%alb_s(1)
            write(*,"(a,f14.4)") " final Tsrf     [K] : ", chn%itm%now%tsrf(1)
            write(*,"(a,f14.4)") " cum smb (column)   : ", real(chn%itm%now%smb_cum(1),wp)
            write(*,"(a,f14.4)") " cum smb_ice        : ", real(chn%itm%now%smbi_cum(1),wp)
            write(*,"(a,f14.4)") " cum runoff         : ", real(chn%itm%now%runoff_cum(1),wp)
            write(*,"(a,f14.4)") " cum melt           : ", real(chn%itm%now%melt_cum(1),wp)
            write(*,"(a,f14.4)") " cum refreezing     : ", real(chn%itm%now%refrz_cum(1),wp)

    end select

    write(*,"(a)") "==========================================================="
    write(*,"(a)") ""

    call chion_end(chn)

    deallocate(smb)
    deallocate(thickness)
    deallocate(wet_mass)
    deallocate(bulk_density)
    deallocate(liquid_water)

contains

    function column_day_of_year(t) result(doy)
        ! 1-based fractional day of year of the model time `t` [d].

        implicit none

        real(wp), intent(IN) :: t
        real(wp) :: doy

        doy = modulo(t,YEAR_LENGTH) + 1.0_wp

        return

    end function column_day_of_year

    function column_air_temperature(doy) result(t_air)

        implicit none

        real(wp), intent(IN) :: doy
        real(wp) :: t_air

        t_air = t2m_mean + t2m_amplitude*cos(2.0_wp*PI*(doy - t2m_phase)/YEAR_LENGTH)

        return

    end function column_air_temperature

    function column_snow_depth(chn,thickness,wet_mass,bulk_density,liquid_water) result(depth)
        ! The model's own notion of "how much snow is there", for the peak
        ! diagnostic. Deliberately NOT unified across models: BESSI resolves a
        ! layered column and reports a thickness in m, while PDD and ITM carry
        ! a bulk reservoir in kg m-2 / mm w.e. Reporting one number under three
        ! meanings is why the summary labels which one it printed.

        implicit none

        type(chion_class), intent(IN)  :: chn
        real(wp),          intent(OUT) :: thickness(:)
        real(wp),          intent(OUT) :: wet_mass(:)
        real(wp),          intent(OUT) :: bulk_density(:)
        real(wp),          intent(OUT) :: liquid_water(:)
        real(wp) :: depth

        thickness    = 0.0_wp
        wet_mass     = 0.0_wp
        bulk_density = 0.0_wp
        liquid_water = 0.0_wp

        select case(trim(chn%par%model))

            case("bessi")
                call summarize_domain_state(chn%bsi%now%mass,chn%bsi%now%mass_w, &
                                            chn%bsi%now%density,chn%bsi%now%n_lay, &
                                            thickness,wet_mass,bulk_density,liquid_water)
                depth = maxval(thickness)

            case("pdd")
                depth = maxval(chn%pdd%now%snowpack_swe)

            case("itm")
                depth = maxval(chn%itm%now%H_snow)

            case DEFAULT
                depth = 0.0_wp

        end select

        return

    end function column_snow_depth

end program chion_column
