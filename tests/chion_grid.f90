program chion_grid
    ! Gridded driver, `chion_grid.x`.
    !
    !     libchion/bin/chion_grid.x par/chion_grid.nml
    !
    ! Exactly one command-line argument, the parameter file; output written to
    ! the current working directory (the runme contract).
    !
    ! Two forcing sources, selected by ctrl:forcing_source:
    !
    !   "file"   (default, the original behaviour) -- a NetCDF forcing file with
    !            the common (time,y,x) layout, mirroring Chion.jl's
    !            load_forcing_file (src/dataloaders.jl:185). netCDF's dimension
    !            order is C-order, so a variable declared (time,y,x) in the file
    !            arrives in Fortran as (nx,ny,ntime) and one time slice is read
    !            with start=[1,1,it], count=[nx,ny,1]. No transposition anywhere.
    !            The driver walks the file's time axis once.
    !
    !   "domain" -- a standardized monthly climatology assembled straight from
    !            the raw ~/models/ice_data datasets by the driver-layer domain
    !            loader (libs/domains/chion_domain.f90). No preprocessing step.
    !            The 12 monthly means are turned into a daily series by the
    !            mean-preserving monthly->daily interface (chion_forcing_monthly)
    !            and the SAME annual cycle is repeated for ctrl:n_years to spin
    !            the snowpack up toward a steady seasonal state. This is the
    !            Greenland/Antarctica steady-state-snowpack path.
    !
    ! Both sources feed ONE run loop: they differ only in how the per-column
    ! forcing is produced each step. Column building, the model, output, restart
    ! and timing are shared.
    !
    ! ---------------------------------------------------------------------
    ! Columns are built from a mask
    ! ---------------------------------------------------------------------
    ! chion's state is a packed LIST of independent columns, not a 2-D field
    ! (docs/PLAN.md section 1). Cell (ix,iy) becomes a column iff
    ! mask(ix,iy) > mask_threshold. Masked-out cells are never allocated, never
    ! stepped and never written -- they appear in the output as the missing
    ! value. The column ORDER is x-fastest within y, Fortran's own memory order
    ! over the mask, so the list has the best achievable locality for the loop.
    !
    ! ---------------------------------------------------------------------
    ! The &ctrl group this program needs
    ! ---------------------------------------------------------------------
    ! The &ctrl group is the DRIVER's own schema (chion never reads it) and is
    ! not validated against the library defaults. Keys are required WITHIN the
    ! forcing_source branch that uses them, so a "file" par file need not carry
    ! the "domain" keys and vice versa.
    !
    !   &ctrl
    !       forcing_source = "file"        ! "file" | "domain"
    !       file_out       = "chion_grid.nc"
    !       mask_threshold = 0.0           ! cells with mask > threshold -> column
    !       wind_default   = 5.0           ! [m s-1] no wind variable is read
    !       dust_dep_default = 0.0         ! [kg m-2 s-1] uniform dust deposition
    !       rh_default     = 0.0           ! [1] uniform relative humidity, 0 = off
    !       dt_out         = 30.0          ! [d] output interval
    !
    !     ! --- forcing_source = "file" ---
    !       file_forcing   = "forcing.nc"
    !       name_x  = "x"   name_y = "y"   name_time = "time"
    !       name_t2m = "TT" name_sf = "SF" name_rf = "RF" name_swd = "SWD"
    !       name_mask = "mask"  name_lat = "LAT"  name_zs = "SH"
    !       t2m_in_celsius     = .FALSE.   ! convert TT by +273.15
    !       precip_in_mmwe_day = .FALSE.   ! convert SF/RF by /86400
    !       dt             = -1.0          ! [d] <=0 -> infer from the time axis
    !
    !     ! --- forcing_source = "domain" ---
    !       domain         = "greenland"   ! "greenland" | "antarctica"
    !       grid_name      = "GRL-16KM"
    !       path_ice_data  = "/path/to/ice_data"
    !       path_insol     = "libs/insol/input"
    !       path_racmo     = "/path/to/racmo"  ! antarctica only (RACMO climatology root)
    !       n_years        = 50            ! annual cycles to repeat
    !       swd_source     = "file"        ! "file" | "transmissivity" | "transmissivity_seasonal"
    !       trans_a        = 0.46          ! tau = trans_a + trans_b*z_srf (+ trans_c*tcc)
    !       trans_b        = 6.0e-5        ! [m-1]
    !       trans_c        = 0.0           ! [1] cloud term, used by transmissivity_seasonal
    !       H_ice_default  = 1000.0        ! [m] ice thickness set on every column (ITM)
    !   /

    use chion
    use chion_io,     only : chion_write_init, chion_write_step, chion_restart_write
    use chion_domain, only : chion_domain_class, chion_domain_load
    use ncio
    use nml, only : nml_read

    implicit none

    ! --- &ctrl : shared ---------------------------------------------------
    character(len=512) :: file_out
    character(len=56)  :: forcing_source
    real(wp) :: mask_threshold, wind_default, dt_out
    real(wp) :: dust_dep_default
    real(wp) :: rh_default

    ! --- &ctrl : file source ----------------------------------------------
    character(len=512) :: file_forcing
    character(len=56)  :: name_x, name_y, name_time
    character(len=56)  :: name_t2m, name_sf, name_rf, name_swd
    character(len=56)  :: name_mask, name_lat, name_zs
    real(wp) :: dt
    logical  :: t2m_in_celsius, precip_in_mmwe_day

    ! --- &ctrl : domain source --------------------------------------------
    character(len=56)  :: domain, grid_name, swd_source
    character(len=512) :: path_ice_data, path_insol, path_racmo
    integer  :: n_years
    real(wp) :: trans_a, trans_b, trans_c, H_ice_default

    ! --- Driver state -----------------------------------------------------
    type(chion_class) :: chn

    character(len=512) :: path_par
    character(len=*), parameter :: file_restart = "chion_grid_restart.nc"

    logical  :: is_domain
    integer  :: nx, ny, nt, ncol, nmon, nday_year
    integer  :: ix, iy, it, i, n_out, iyear, doy
    integer  :: clock0, clock1, clock_rate

    real(wp) :: year_length
    real(wp) :: time, time_first, time_next_out, dt_use

    real(wp), allocatable :: xc(:), yc(:), times(:)
    real(wp), allocatable :: mask2D(:,:), lat2D(:,:), zs2D(:,:)
    real(wp), allocatable :: t2m(:,:), sf(:,:), rf(:,:), swd(:,:)
    integer,  allocatable :: col_is(:), col_js(:)
    real(wp), allocatable :: maskT(:,:)

    ! Domain / monthly path
    type(chion_domain_class)     :: dom
    type(monthly_to_daily_class) :: md
    real(wp), allocatable :: t2m_c(:,:), sf_c(:,:), rf_c(:,:), swd_c(:,:), tcc_c(:,:)
    real(wp), allocatable :: S_toa_c(:,:)
    real(wp), allocatable :: fday(:), smb_step(:), export_year(:), mass_prev(:)
    real(wp), allocatable :: mass_now(:)

    ! =====================================================================
    ! Setup: parameters
    ! =====================================================================

    call chion_load_command_line_args(path_par)
    call chion_check_file(path_par)

    nmon      = 12
    nday_year = 360

    ! Shared keys.
    call nml_read(path_par,"ctrl","forcing_source", forcing_source)
    call nml_read(path_par,"ctrl","file_out",       file_out)
    call nml_read(path_par,"ctrl","mask_threshold", mask_threshold)
    call nml_read(path_par,"ctrl","wind_default",   wind_default)
    call nml_read(path_par,"ctrl","dust_dep_default", dust_dep_default)
    call nml_read(path_par,"ctrl","rh_default",     rh_default)
    call nml_read(path_par,"ctrl","dt_out",         dt_out)

    is_domain = (trim(forcing_source) .eq. "domain")

    if (.not. is_domain .and. trim(forcing_source) .ne. "file") then
        write(io_unit_err,*) "chion_grid:: Error: forcing_source must be 'file' or 'domain', got '"// &
                             trim(forcing_source)//"'."
        stop "Program stopped."
    end if

    if (is_domain) then
        call nml_read(path_par,"ctrl","domain",        domain)
        call nml_read(path_par,"ctrl","grid_name",     grid_name)
        call nml_read(path_par,"ctrl","path_ice_data", path_ice_data)
        call nml_read(path_par,"ctrl","path_insol",    path_insol)
        call nml_read(path_par,"ctrl","n_years",       n_years)
        ! path_racmo: the RACMO climatology root, required only for Antarctica.
        ! Read conditionally so Greenland/other par files need not carry it.
        path_racmo = ""
        if (trim(domain) .eq. "antarctica") &
            call nml_read(path_par,"ctrl","path_racmo", path_racmo)
        call nml_read(path_par,"ctrl","swd_source",    swd_source)
        call nml_read(path_par,"ctrl","trans_a",       trans_a)
        call nml_read(path_par,"ctrl","trans_b",       trans_b)
        call nml_read(path_par,"ctrl","trans_c",       trans_c)
        call nml_read(path_par,"ctrl","H_ice_default", H_ice_default)
    else
        call nml_read(path_par,"ctrl","file_forcing",      file_forcing)
        call nml_read(path_par,"ctrl","name_x",            name_x)
        call nml_read(path_par,"ctrl","name_y",            name_y)
        call nml_read(path_par,"ctrl","name_time",         name_time)
        call nml_read(path_par,"ctrl","name_t2m",          name_t2m)
        call nml_read(path_par,"ctrl","name_sf",           name_sf)
        call nml_read(path_par,"ctrl","name_rf",           name_rf)
        call nml_read(path_par,"ctrl","name_swd",          name_swd)
        call nml_read(path_par,"ctrl","name_mask",         name_mask)
        call nml_read(path_par,"ctrl","name_lat",          name_lat)
        call nml_read(path_par,"ctrl","name_zs",           name_zs)
        call nml_read(path_par,"ctrl","t2m_in_celsius",    t2m_in_celsius)
        call nml_read(path_par,"ctrl","precip_in_mmwe_day",precip_in_mmwe_day)
        call nml_read(path_par,"ctrl","dt",                dt)
    end if

    ! =====================================================================
    ! Setup: geometry, static fields, timing
    ! =====================================================================

    if (is_domain) then

        call chion_domain_load(dom, trim(domain), trim(grid_name), &
                               trim(path_ice_data), trim(path_insol), &
                               nmon=nmon, nday_year=nday_year, &
                               path_racmo=trim(path_racmo))

        nx = dom%nx
        ny = dom%ny

        allocate(xc(nx), yc(ny))
        xc = dom%xc
        yc = dom%yc

        allocate(mask2D(nx,ny), lat2D(nx,ny), zs2D(nx,ny))
        mask2D = dom%mask
        lat2D  = dom%lat2D
        zs2D   = dom%z_srf

        year_length = real(nday_year,wp)
        dt_use      = 1.0_wp                       ! daily steps
        nt          = n_years*nday_year
        time_first  = 0.0_wp

    else

        call chion_check_file(file_forcing)

        nx = nc_size(file_forcing,trim(name_x))
        ny = nc_size(file_forcing,trim(name_y))
        nt = nc_size(file_forcing,trim(name_time))

        allocate(xc(nx), yc(ny), times(nt))
        call nc_read(file_forcing,trim(name_x),   xc)
        call nc_read(file_forcing,trim(name_y),   yc)
        call nc_read(file_forcing,trim(name_time),times)

        ! Timestep. Preferred source is the time axis itself (Chion.jl's
        ! infer_dt_days); ctrl:dt overrides when positive.
        if (dt .gt. 0.0_wp) then
            dt_use = dt
        else if (nt .ge. 2) then
            dt_use = times(2) - times(1)
        else
            write(io_unit_err,*) "chion_grid:: Error: cannot infer dt from a single time record."
            write(io_unit_err,*) "Set ctrl:dt explicitly."
            stop "Program stopped."
        end if

        if (dt_use .le. 0.0_wp) then
            write(io_unit_err,*) "chion_grid:: Error: inferred dt is not positive: ", dt_use
            stop "Program stopped."
        end if

        year_length = 365.0_wp
        time_first  = times(1)

        ! (y,x) in the file arrives as (nx,ny) in Fortran.
        allocate(mask2D(nx,ny), lat2D(nx,ny), zs2D(nx,ny))

        if (trim(name_mask) .eq. "None") then
            mask2D = 1.0_wp
        else
            call nc_read(file_forcing,trim(name_mask),mask2D)
        end if

        lat2D = 0.0_wp
        if (trim(name_lat) .ne. "None") call nc_read(file_forcing,trim(name_lat),lat2D)

        zs2D = 0.0_wp
        if (trim(name_zs) .ne. "None") call nc_read(file_forcing,trim(name_zs),zs2D)

    end if

    ! =====================================================================
    ! Setup: column list (shared)
    ! =====================================================================

    ncol = 0
    do iy = 1, ny
    do ix = 1, nx
        if (mask2D(ix,iy) .gt. mask_threshold) ncol = ncol + 1
    end do
    end do

    if (ncol .eq. 0) then
        write(io_unit_err,*) "chion_grid:: Error: the mask selects no columns."
        write(io_unit_err,*) "mask_threshold = ", mask_threshold
        stop "Program stopped."
    end if

    allocate(col_is(ncol), col_js(ncol))
    i = 0
    do iy = 1, ny
    do ix = 1, nx
        if (mask2D(ix,iy) .gt. mask_threshold) then
            i = i + 1
            col_is(i) = ix
            col_js(i) = iy
        end if
    end do
    end do

    ! =====================================================================
    ! Setup: chion
    ! =====================================================================

    call chion_init(chn,path_par,ncol)
    call chion_init_state(chn)

    ! chion_grid_class carries the mask as (ny,nx), the Chion.jl orientation;
    ! everything above is (nx,ny). Transposed once, here.
    allocate(maskT(ny,nx))
    maskT = transpose(mask2D)
    call chion_set_grid(chn,xc,yc,col_js,col_is,mask=maskT)

    ! Time-invariant per-column forcing.
    do i = 1, ncol
        chn%forc%latitude_deg(i)   = lat2D(col_is(i),col_js(i))
        chn%forc%surface_height(i) = zs2D(col_is(i),col_js(i))
    end do
    chn%forc%wind_speed = wind_default

    ! Uniform dust deposition, for SEMIX-albedo dust sensitivity experiments.
    ! A spatially/temporally varying dust field would come from a domain loader
    ! (paleo) or from the coupler; this is the single-value sensitivity knob.
    ! Left flagged off at zero so the albedo skips the dust path entirely.
    chn%forc%dust_dep     = dust_dep_default
    chn%forc%has_dust_dep = (dust_dep_default .gt. 0.0_wp)

    ! Uniform relative humidity. No domain loader carries a humidity field yet,
    ! so without this knob has_relative_humidity is false everywhere and the
    ! turbulent latent flux is identically zero -- under BOTH surface schemes.
    ! That is the standing state of every GRL benchmark to date, and it is why
    ! the seb_scheme comparison in docs/semix_port_scope.md is a
    ! sensible-heat-only result.
    !
    ! Left flagged OFF at zero so existing par files are unaffected. Note the
    ! two schemes read the same number differently: BESSI takes it relative to
    ! saturation over WATER (energy_flux.jl:57-60), SEMIX over ICE
    ! (semi.f90:201). Both are their own source's reading, so a run that varies
    ! rh_default across seb_scheme is not a controlled comparison of the
    ! turbulent exchange alone.
    chn%forc%relative_humidity     = rh_default
    chn%forc%has_relative_humidity = (rh_default .gt. 0.0_wp)

    ! =====================================================================
    ! Setup: per-source forcing buffers
    ! =====================================================================

    if (is_domain) then

        ! Scatter monthly gridded fields onto the column list, then build the
        ! mean-preserving control values ONCE (every field, so PDD sums and
        ! precip totals are both honoured). S_toa is already daily.
        call monthly_to_daily_init(md, nmon, nday_year/nmon)

        allocate(t2m_c(ncol,nmon), sf_c(ncol,nmon), rf_c(ncol,nmon))
        allocate(swd_c(ncol,nmon), tcc_c(ncol,nmon))
        allocate(S_toa_c(ncol,nday_year))

        do i = 1, ncol
            t2m_c(i,:) = dom%t2m(col_is(i),col_js(i),:)
            sf_c(i,:)  = dom%sf (col_is(i),col_js(i),:)
            rf_c(i,:)  = dom%rf (col_is(i),col_js(i),:)
            swd_c(i,:) = dom%swd(col_is(i),col_js(i),:)
            tcc_c(i,:) = dom%tcc(col_is(i),col_js(i),:)
            S_toa_c(i,:) = dom%S_toa(col_is(i),col_js(i),:)
        end do

        ! Replace monthly means with mean-preserving control values. A scratch
        ! copy avoids aliasing the intent(IN)/intent(OUT) pair.
        block
            real(wp), allocatable :: mtmp(:,:)
            allocate(mtmp(ncol,nmon))
            mtmp = t2m_c ; call monthly_to_daily_controls(md, mtmp, t2m_c)
            mtmp = sf_c  ; call monthly_to_daily_controls(md, mtmp, sf_c)
            mtmp = rf_c  ; call monthly_to_daily_controls(md, mtmp, rf_c)
            mtmp = swd_c ; call monthly_to_daily_controls(md, mtmp, swd_c)
            mtmp = tcc_c ; call monthly_to_daily_controls(md, mtmp, tcc_c)
            deallocate(mtmp)
        end block

        allocate(fday(ncol), smb_step(ncol), export_year(ncol), mass_prev(ncol))
        export_year = 0.0_wp
        mass_prev   = chion_column_mass(chn)     ! cold start: 0

        ! ITM ice thickness (BESSI and PDD ignore it). Annual PDDs is set once,
        ! from the repeating climatology, below.
        chn%forc%H_ice = H_ice_default
        call domain_set_annual_pdds(chn, md, t2m_c, chn%c%T0, dt_use)

    else

        allocate(t2m(nx,ny), sf(nx,ny), rf(nx,ny), swd(nx,ny))

    end if

    ! =====================================================================
    ! Header
    ! =====================================================================

    write(*,"(a)")        "== chion_grid ============================================="
    write(*,"(a,a)")      " par file      : ", trim(path_par)
    write(*,"(a,a)")      " forcing source: ", trim(forcing_source)
    if (is_domain) then
        write(*,"(a,a,a,a)") " domain / grid : ", trim(domain), " / ", trim(grid_name)
        write(*,"(a,a)")     " swd source    : ", trim(swd_source)
        write(*,"(a,i0)")    " years         : ", n_years
    else
        write(*,"(a,a)")     " forcing file  : ", trim(file_forcing)
    end if
    write(*,"(a,a)")      " output        : ", trim(file_out)
    write(*,"(a,i0,a,i0,a,i0)") " grid          : nx=", nx, "  ny=", ny, "  steps=", nt
    write(*,"(a,i0,a,i0,a,f6.2,a)") " columns       : ncol=", ncol, " of ", nx*ny, &
                                    "  (", 100.0_wp*real(ncol,wp)/real(nx*ny,wp), " %)"
    write(*,"(a,f10.4,a,f10.4)") " dt, dt_out    : ", dt_use, "  ", dt_out
    write(*,"(a)")        "==========================================================="

    ! =====================================================================
    ! Output record 1: the INITIAL state, before any step.
    ! =====================================================================

    call chion_write_init(chn,file_out,time_first,"days")
    call chion_write_step(chn,file_out,time_first)

    time_next_out = time_first + dt_out
    n_out         = 1

    call system_clock(clock0,clock_rate)

    ! =====================================================================
    ! Run
    ! =====================================================================

    time = time_first

    do it = 1, nt

        ! --- Populate the per-column forcing for this step ---------------
        if (is_domain) then

            doy = mod(it-1, nday_year) + 1

            call interp_monthly_to_day(md, t2m_c, doy, fday)
            chn%forc%air_temperature = fday
            call interp_monthly_to_day(md, sf_c,  doy, fday)
            chn%forc%snowfall_rate   = fday
            call interp_monthly_to_day(md, rf_c,  doy, fday)
            chn%forc%rainfall_rate   = fday

            select case(trim(swd_source))
                case("file")
                    call interp_monthly_to_day(md, swd_c, doy, fday)
                    chn%forc%shortwave_down = fday
                case("transmissivity")
                    ! tau(z_srf) * S_toa(day), the ITM clear-sky form. z_srf is
                    ! the per-column surface_height set above.
                    chn%forc%shortwave_down = (trans_a + trans_b*chn%forc%surface_height) &
                                              * S_toa_c(:,doy)
                case("transmissivity_seasonal")
                    ! Cloud-aware transmissivity calibrated on the melt season:
                    !   tau = trans_a + trans_b*z_srf + trans_c*tcc,  clamped to
                    ! [0,1], times S_toa. tcc is the daily-interpolated cloud
                    ! cover (fday is reused as the tcc-of-day scratch here). The
                    ! diagnostic (diagnostics/transmissivity.jl) shows a single
                    ! annual (a,b,c) is misleading -- season confounds it -- so
                    ! these default to the JJA fit; see docs/steady_state_snowpack.md.
                    call interp_monthly_to_day(md, tcc_c, doy, fday)
                    chn%forc%shortwave_down = max(0.0_wp, min(1.0_wp, &
                        trans_a + trans_b*chn%forc%surface_height + trans_c*fday)) &
                        * S_toa_c(:,doy)
                case DEFAULT
                    write(io_unit_err,*) "chion_grid:: Error: swd_source must be 'file', "// &
                        "'transmissivity' or 'transmissivity_seasonal', got '"// &
                        trim(swd_source)//"'."
                    stop "Program stopped."
            end select

            chn%forc%day_of_year         = real(doy,wp)
            chn%forc%solar_longitude_deg = 360.0_wp*(real(doy,wp) - 1.0_wp)/year_length

        else

            time = times(it)

            call nc_read(file_forcing,trim(name_t2m),t2m,start=[1,1,it],count=[nx,ny,1])
            call nc_read(file_forcing,trim(name_sf), sf, start=[1,1,it],count=[nx,ny,1])
            call nc_read(file_forcing,trim(name_rf), rf, start=[1,1,it],count=[nx,ny,1])
            call nc_read(file_forcing,trim(name_swd),swd,start=[1,1,it],count=[nx,ny,1])

            if (t2m_in_celsius)     t2m = t2m + 273.15_wp
            if (precip_in_mmwe_day) then
                sf = sf/86400.0_wp
                rf = rf/86400.0_wp
            end if

            do i = 1, ncol
                chn%forc%air_temperature(i) = t2m(col_is(i),col_js(i))
                chn%forc%snowfall_rate(i)   = sf(col_is(i),col_js(i))
                chn%forc%rainfall_rate(i)   = rf(col_is(i),col_js(i))
                chn%forc%shortwave_down(i)  = swd(col_is(i),col_js(i))
            end do

            chn%forc%day_of_year         = modulo(time,year_length) + 1.0_wp
            chn%forc%solar_longitude_deg = 360.0_wp*(chn%forc%day_of_year - 1.0_wp)/year_length

        end if

        ! --- Step -------------------------------------------------------
        call chion_update(chn,dt_use)

        ! --- Annual surface-SMB accounting (domain mode) ----------------
        ! chion_get_smb is the ICE-FACING flux (firn exported to the ice
        ! sheet), which is ~0 while a cold-start column is still filling. The
        ! MAR-comparable SURFACE mass balance is the column-storage tendency
        ! plus that export:  SMB_surf = d(column_mass)/dt + smb_ice.
        ! Integrated over a year this is a mass in kg m-2 == mm w.e. yr-1.
        if (is_domain) then
            call chion_get_smb(chn,smb_step)
            export_year = export_year + smb_step*chn%c%seconds_per_day*dt_use
            doy = mod(it-1, nday_year) + 1
            if (doy .eq. nday_year) then
                iyear    = it/nday_year
                mass_now = chion_column_mass(chn)
                write(*,"(a,i4,a,f9.2,a,f9.2,a,f10.2)") " year ", iyear, &
                    "  surface SMB [mm/yr] = ", &
                    mean_col((mass_now - mass_prev) + export_year), &
                    "   ice export [mm/yr] = ", mean_col(export_year), &
                    "   column mass [kg/m2] = ", mean_col(mass_now)
                mass_prev   = mass_now
                export_year = 0.0_wp
            end if
        end if

        ! --- Output at dt_out cadence -----------------------------------
        ! State now in chn is valid at time+dt_use (END of the step); stamp it
        ! there, not at `time`, so the series is not shifted one step early.
        if (time + dt_use .ge. time_next_out - 0.5_wp*dt_use) then
            call chion_write_step(chn,file_out,time+dt_use)
            n_out         = n_out + 1
            time_next_out = time_next_out + dt_out
        end if

        if (is_domain) time = time + dt_use

    end do

    call system_clock(clock1)

    ! State after the final step, stamped at the end of that step.
    call chion_restart_write(chn,file_restart,time+dt_use)

    write(*,"(a)") ""
    write(*,"(a)") "== chion_grid summary ====================================="
    write(*,"(a,a)")     " model            : ", trim(chn%par%model)
    write(*,"(a,i0)")    " columns          : ", ncol
    write(*,"(a,i0)")    " steps            : ", nt
    write(*,"(a,i0)")    " output records   : ", n_out
    write(*,"(a,f12.4)") " wall time    [s] : ", real(clock1-clock0,wp)/real(clock_rate,wp)
    write(*,"(a,es12.4)")" per column-step  : ", real(clock1-clock0,wp)/real(clock_rate,wp) &
                                                 /real(ncol,wp)/real(nt,wp)
    write(*,"(a)") "==========================================================="
    write(*,"(a)") ""

    call chion_end(chn)

contains

    function mean_col(v) result(m)
        real(wp), intent(IN) :: v(:)
        real(wp) :: m
        m = sum(v)/real(size(v),wp)
    end function mean_col

    function chion_column_mass(chn) result(mass)
        ! Total column snow mass [kg m-2] per column, for a coarse steady-state
        ! drift indicator. BESSI sums its layer masses; PDD/ITM report their
        ! bulk snow reservoir.
        type(chion_class), intent(IN) :: chn
        real(wp) :: mass(chn%grd%ncol)
        integer  :: j

        select case(trim(chn%par%model))
            case("bessi")
                do j = 1, chn%grd%ncol
                    mass(j) = sum(chn%bsi%now%mass(1:chn%bsi%now%n_lay(j),j)) &
                              + sum(chn%bsi%now%mass_w(1:chn%bsi%now%n_lay(j),j))
                end do
            case("pdd")
                mass = chn%pdd%now%snowpack_swe
            case("itm")
                mass = chn%itm%now%H_snow
            case DEFAULT
                mass = 0.0_wp
        end select
    end function chion_column_mass

    subroutine domain_set_annual_pdds(chn, md, t2m_ctrl, T0, dt_use)
        ! Annual positive degree days per column, integrated from the daily
        ! temperature of the repeating climatology. ITM reads forc%PDDs to
        ! interpolate its critical snow depth; it is a whole-year quantity, set
        ! once (smbpal's calc_pdds), not a per-step increment. BESSI and PDD
        ! ignore it.
        type(chion_class),            intent(INOUT) :: chn
        type(monthly_to_daily_class), intent(IN)    :: md
        real(wp), intent(IN) :: t2m_ctrl(:,:)      ! (ncol,nmon) control values
        real(wp), intent(IN) :: T0
        real(wp), intent(IN) :: dt_use

        real(wp), allocatable :: tday(:)
        integer :: d, ncol

        ncol = size(t2m_ctrl,1)
        allocate(tday(ncol))

        chn%forc%PDDs = 0.0_wp
        do d = 1, md%nday_year
            call interp_monthly_to_day(md, t2m_ctrl, d, tday)
            where (tday .gt. T0) chn%forc%PDDs = chn%forc%PDDs + (tday - T0)*dt_use
        end do

        deallocate(tday)
    end subroutine domain_set_annual_pdds

end program chion_grid
