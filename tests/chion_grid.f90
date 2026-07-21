program chion_grid
    ! WP15 -- gridded driver, `chion_grid.x`.
    !
    !     libchion/bin/chion_grid.x par/chion_grid.nml
    !
    ! Exactly one command-line argument, the parameter file; output written to
    ! the current working directory (the runme contract).
    !
    ! Forcing comes from a NetCDF file, mirroring Chion.jl's load_forcing_file
    ! (src/dataloaders.jl:185). The common `(time, y, x)` layout is supported:
    ! netCDF's dimension order is C-order, so a variable declared `(time,y,x)`
    ! in the file arrives in Fortran as `(nx,ny,ntime)` and one time slice is
    ! read with start=[1,1,it], count=[nx,ny,1]. No transposition is needed
    ! anywhere.
    !
    ! ---------------------------------------------------------------------
    ! Columns are built from a mask
    ! ---------------------------------------------------------------------
    ! chion's state is a packed LIST of independent columns, not a 2-D field
    ! (docs/PLAN.md section 1). This driver builds that list from the mask
    ! variable: cell (ix,iy) becomes a column iff mask(ix,iy) > mask_threshold.
    ! Masked-out cells are never allocated, never stepped and never written --
    ! they appear in the output as the missing value, which is the whole point
    ! of the column list. Set mask_name = "None" to take every cell.
    !
    ! The column ORDER is x-fastest within y, which is Fortran's own memory
    ! order over the mask, so a column list built this way has the best
    ! achievable locality for the OpenMP loop.
    !
    ! ---------------------------------------------------------------------
    ! The &ctrl group this program needs
    ! ---------------------------------------------------------------------
    ! There is no par/chion_grid.nml in the repository yet -- this driver's
    ! &ctrl group is listed here so one can be written. Every entry is
    ! required: the group is read without a defaults file, so it is its own
    ! schema and a missing entry is a hard error rather than a silent zero.
    !
    !   &ctrl
    !       file_forcing   = "forcing.nc"  ! Forcing file
    !       file_out       = "chion_grid.nc"
    !
    !       name_x         = "x"           ! Coordinate variables
    !       name_y         = "y"
    !       name_time      = "time"
    !
    !       name_t2m       = "TT"          ! (time,y,x) forcing variables
    !       name_sf        = "SF"
    !       name_rf        = "RF"
    !       name_swd       = "SWD"
    !
    !       name_mask      = "mask"        ! (y,x) static, or "None"
    !       mask_threshold = 0.0
    !       name_lat       = "LAT"         ! (y,x) static, or "None"
    !       name_zs        = "SH"          ! (y,x) static, or "None"
    !
    !       t2m_in_celsius     = False     ! Convert TT by +273.15
    !       precip_in_mmwe_day = False     ! Convert SF/RF by /86400
    !       wind_default       = 5.0       ! [m s-1] no wind variable is read
    !
    !       dt             = -1.0          ! [d] <= 0 -> infer from the time axis
    !       dt_out         = 30.0          ! [d]
    !   /
    !
    ! Defaults for the names follow load_forcing_file's keyword defaults, so a
    ! file written for Chion.jl is readable here unchanged.

    use chion
    use chion_io, only : chion_write_init, chion_write_step, chion_restart_write
    use ncio
    use nml, only : nml_read

    implicit none

    ! --- &ctrl ------------------------------------------------------------
    character(len=512) :: file_forcing, file_out
    character(len=56)  :: name_x, name_y, name_time
    character(len=56)  :: name_t2m, name_sf, name_rf, name_swd
    character(len=56)  :: name_mask, name_lat, name_zs
    real(wp) :: mask_threshold, wind_default, dt, dt_out
    logical  :: t2m_in_celsius, precip_in_mmwe_day

    ! --- Driver state -----------------------------------------------------
    type(chion_class) :: chn

    character(len=512) :: path_par
    character(len=*), parameter :: file_restart = "chion_grid_restart.nc"

    real(wp), parameter :: YEAR_LENGTH = 365.0_wp

    integer :: nx, ny, nt, ncol
    integer :: ix, iy, it, i, n_out

    real(wp), allocatable :: xc(:), yc(:), times(:)
    real(wp), allocatable :: mask2D(:,:), lat2D(:,:), zs2D(:,:)
    real(wp), allocatable :: t2m(:,:), sf(:,:), rf(:,:), swd(:,:)
    integer,  allocatable :: col_is(:), col_js(:)
    real(wp), allocatable :: maskT(:,:)

    real(wp) :: time, time_next_out, dt_use
    integer  :: clock0, clock1, clock_rate

    ! =====================================================================
    ! Setup
    ! =====================================================================

    call chion_load_command_line_args(path_par)
    call chion_check_file(path_par)

    call nml_read(path_par,"ctrl","file_forcing",      file_forcing)
    call nml_read(path_par,"ctrl","file_out",          file_out)
    call nml_read(path_par,"ctrl","name_x",            name_x)
    call nml_read(path_par,"ctrl","name_y",            name_y)
    call nml_read(path_par,"ctrl","name_time",         name_time)
    call nml_read(path_par,"ctrl","name_t2m",          name_t2m)
    call nml_read(path_par,"ctrl","name_sf",           name_sf)
    call nml_read(path_par,"ctrl","name_rf",           name_rf)
    call nml_read(path_par,"ctrl","name_swd",          name_swd)
    call nml_read(path_par,"ctrl","name_mask",         name_mask)
    call nml_read(path_par,"ctrl","mask_threshold",    mask_threshold)
    call nml_read(path_par,"ctrl","name_lat",          name_lat)
    call nml_read(path_par,"ctrl","name_zs",           name_zs)
    call nml_read(path_par,"ctrl","t2m_in_celsius",    t2m_in_celsius)
    call nml_read(path_par,"ctrl","precip_in_mmwe_day",precip_in_mmwe_day)
    call nml_read(path_par,"ctrl","wind_default",      wind_default)
    call nml_read(path_par,"ctrl","dt",                dt)
    call nml_read(path_par,"ctrl","dt_out",            dt_out)

    call chion_check_file(file_forcing)

    ! --- Geometry ---------------------------------------------------------

    nx = nc_size(file_forcing,trim(name_x))
    ny = nc_size(file_forcing,trim(name_y))
    nt = nc_size(file_forcing,trim(name_time))

    allocate(xc(nx))
    allocate(yc(ny))
    allocate(times(nt))

    call nc_read(file_forcing,trim(name_x),   xc)
    call nc_read(file_forcing,trim(name_y),   yc)
    call nc_read(file_forcing,trim(name_time),times)

    ! Timestep. Preferred source is the time axis itself, exactly as
    ! Chion.jl's infer_dt_days; ctrl:dt overrides it when positive, which is
    ! what a file with a calendar time axis (months of unequal length) needs.
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
        write(io_unit_err,*) "chion_grid:: Error: inferred dt is not positive."
        write(io_unit_err,*) "dt = ", dt_use
        stop "Program stopped."
    end if

    ! --- Static fields ----------------------------------------------------
    ! (y,x) in the file arrives as (nx,ny) in Fortran.

    allocate(mask2D(nx,ny))
    allocate(lat2D(nx,ny))
    allocate(zs2D(nx,ny))

    if (trim(name_mask) .eq. "None") then
        mask2D = 1.0_wp
    else
        call nc_read(file_forcing,trim(name_mask),mask2D)
    end if

    lat2D = 0.0_wp
    if (trim(name_lat) .ne. "None") call nc_read(file_forcing,trim(name_lat),lat2D)

    zs2D = 0.0_wp
    if (trim(name_zs) .ne. "None") call nc_read(file_forcing,trim(name_zs),zs2D)

    ! --- The column list --------------------------------------------------

    ncol = 0
    do iy = 1, ny
    do ix = 1, nx
        if (mask2D(ix,iy) .gt. mask_threshold) ncol = ncol + 1
    end do
    end do

    if (ncol .eq. 0) then
        write(io_unit_err,*) "chion_grid:: Error: the mask selects no columns."
        write(io_unit_err,*) "name_mask, mask_threshold = ", trim(name_mask), mask_threshold
        stop "Program stopped."
    end if

    allocate(col_is(ncol))
    allocate(col_js(ncol))

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

    ! --- chion ------------------------------------------------------------

    call chion_init(chn,path_par,ncol)
    call chion_init_state(chn)

    ! chion_grid_class carries the mask as (ny,nx), the Chion.jl orientation;
    ! everything above is (nx,ny), the Fortran-native one. Transposed once,
    ! here, rather than at every use.
    allocate(maskT(ny,nx))
    maskT = transpose(mask2D)

    call chion_set_grid(chn,xc,yc,col_js,col_is,mask=maskT)

    ! Time-invariant per-column forcing.
    do i = 1, ncol
        chn%forc%latitude_deg(i)   = lat2D(col_is(i),col_js(i))
        chn%forc%surface_height(i) = zs2D(col_is(i),col_js(i))
    end do
    chn%forc%wind_speed = wind_default

    allocate(t2m(nx,ny))
    allocate(sf(nx,ny))
    allocate(rf(nx,ny))
    allocate(swd(nx,ny))

    write(*,"(a)")        "== chion_grid ============================================="
    write(*,"(a,a)")      " par file   : ", trim(path_par)
    write(*,"(a,a)")      " forcing    : ", trim(file_forcing)
    write(*,"(a,a)")      " output     : ", trim(file_out)
    write(*,"(a,i0,a,i0,a,i0)") " grid       : nx=", nx, "  ny=", ny, "  nt=", nt
    write(*,"(a,i0,a,i0,a,f6.2,a)") " columns    : ncol=", ncol, " of ", nx*ny, &
                                    "  (", 100.0_wp*real(ncol,wp)/real(nx*ny,wp), " %)"
    write(*,"(a,f10.4,a,f10.4)") " dt, dt_out : ", dt_use, "  ", dt_out
    write(*,"(a)")        "==========================================================="

    call chion_write_init(chn,file_out,times(1),"days")
    call chion_write_step(chn,file_out,times(1))

    time_next_out = times(1) + dt_out
    n_out         = 1

    call system_clock(clock0,clock_rate)

    ! =====================================================================
    ! Run
    ! =====================================================================

    do it = 1, nt

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

        chn%forc%day_of_year         = modulo(time,YEAR_LENGTH) + 1.0_wp
        chn%forc%solar_longitude_deg = 360.0_wp*(chn%forc%day_of_year - 1.0_wp)/YEAR_LENGTH

        call chion_update(chn,dt_use)

        if (time .ge. time_next_out - 0.5_wp*dt_use) then
            call chion_write_step(chn,file_out,time)
            n_out         = n_out + 1
            time_next_out = time_next_out + dt_out
        end if

    end do

    call system_clock(clock1)

    call chion_restart_write(chn,file_restart,time)

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

end program chion_grid
