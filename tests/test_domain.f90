program test_domain
    ! Smoke test for the driver-layer domain loader. Not a unit test with
    ! tolerances -- it loads the real Greenland GRL-16KM datasets, exercises the
    ! conservative ERA5 regridding and the TOA insolation, and reports ranges so
    ! the assembled forcing can be eyeballed for sanity. Writes a small NetCDF
    ! (domain_greenland_check.nc) with the static and annual-mean fields.
    !
    !   test_domain.x <path_ice_data> <path_insol>
    !
    ! e.g. test_domain.x /Users/<you>/models/ice_data libs/insol/input

    use chion_domain
    use chion_defs, only : wp
    use ncio

    implicit none

    type(chion_domain_class) :: dom
    character(len=512) :: path_ice_data, path_insol
    character(len=*), parameter :: file_out = "domain_greenland_check.nc"
    integer :: narg
    integer :: nx, ny
    real(wp) :: msk_thr
    logical, allocatable :: ice(:,:)

    narg = command_argument_count()
    if (narg .lt. 2) then
        write(*,*) "usage: test_domain.x <path_ice_data> <path_insol>"
        stop 1
    end if
    call get_command_argument(1, path_ice_data)
    call get_command_argument(2, path_insol)

    call chion_domain_load(dom, "greenland", "GRL-16KM", trim(path_ice_data), &
                           trim(path_insol), nmon=12, nday_year=360)

    nx = dom%nx
    ny = dom%ny
    msk_thr = 50.0_wp

    allocate(ice(nx,ny))
    ice = dom%mask .gt. msk_thr

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion domain loader smoke test: greenland GRL-16KM"
    write(*,"(a)") "=========================================================="
    write(*,"(a,i0,a,i0)")   " grid            : nx=", nx, "  ny=", ny
    write(*,"(a,i0,a,i0)")   " ice cells (>", int(msk_thr), ") : ", count(ice)
    write(*,*)
    call range_line("z_srf  [m]     ", dom%z_srf, ice)
    call range_line("lat2D  [deg]   ", dom%lat2D, ice)
    write(*,*)
    call range3("t2m    [K]     ", dom%t2m, ice)
    call range3("sf     [kg/m2/s]", dom%sf, ice)
    call range3("rf     [kg/m2/s]", dom%rf, ice)
    call range3("swd    [W/m2]   ", dom%swd, ice)
    call range3("tcc    [1]      ", dom%tcc, ice)
    call range3("S_toa  [W/m2]   ", dom%S_toa, ice)
    write(*,*)
    call range3("smb_ref[kg/m2/s]", dom%smb_ref, ice)
    write(*,"(a,f10.2)") " smb_ref ice mean [mm/yr] : ", &
        360.0*86400.0*mean3(dom%smb_ref, ice)
    write(*,"(a,f10.2)") " swd     ice mean [W/m2]  : ", mean3(dom%swd, ice)
    write(*,"(a,f10.2)") " S_toa   ann max  [W/m2]  : ", maxval(dom%S_toa)
    write(*,*)

    ! Write a check file: static, annual means, AND monthly swd/tcc/S_toa so the
    ! transmissivity diagnostic (diagnostics/transmissivity.jl) can fit a
    ! seasonal tau = swd/S_toa. S_toa is daily; here it is averaged to months.
    call nc_create(file_out)
    call nc_write_dim(file_out, "xc",    x=dom%xc, units="km")
    call nc_write_dim(file_out, "yc",    x=dom%yc, units="km")
    call nc_write_dim(file_out, "month", x=1, dx=1, nx=dom%nmon, units="month")
    call nc_write(file_out, "lon2D", dom%lon2D, dim1="xc", dim2="yc")
    call nc_write(file_out, "lat2D", dom%lat2D, dim1="xc", dim2="yc")
    call nc_write(file_out, "mask",  dom%mask,  dim1="xc", dim2="yc")
    call nc_write(file_out, "z_srf", dom%z_srf, dim1="xc", dim2="yc")
    call nc_write(file_out, "t2m",   annmean(dom%t2m),   dim1="xc", dim2="yc")
    call nc_write(file_out, "smb",   annmean(dom%smb_ref), dim1="xc", dim2="yc")

    ! Monthly, for the transmissivity fit.
    call nc_write(file_out, "swd",   dom%swd, dim1="xc", dim2="yc", dim3="month")
    call nc_write(file_out, "tcc",   dom%tcc, dim1="xc", dim2="yc", dim3="month")
    call nc_write(file_out, "S_toa", monthly_S_toa(dom), dim1="xc", dim2="yc", dim3="month")

    write(*,"(a)") " wrote "//file_out
    write(*,"(a)") "=========================================================="

contains

    function annmean(v) result(m)
        real(wp), intent(IN) :: v(:,:,:)
        real(wp) :: m(size(v,1),size(v,2))
        m = sum(v, dim=3)/real(size(v,3))
    end function annmean

    function monthly_S_toa(dom) result(sm)
        ! Average the daily S_toa (nx,ny,nday_year) into months (nx,ny,nmon),
        ! nday_mon days each, to pair with the monthly swd.
        type(chion_domain_class), intent(IN) :: dom
        real(wp) :: sm(dom%nx, dom%ny, dom%nmon)
        integer  :: m, d0, ndm
        ndm = dom%nday_year/dom%nmon
        do m = 1, dom%nmon
            d0 = (m-1)*ndm
            sm(:,:,m) = sum(dom%S_toa(:,:,d0+1:d0+ndm), dim=3)/real(ndm,wp)
        end do
    end function monthly_S_toa

    subroutine range_line(name, v, ice)
        character(len=*), intent(IN) :: name
        real(wp), intent(IN) :: v(:,:)
        logical,  intent(IN) :: ice(:,:)
        write(*,"(a,a,es12.4,a,es12.4)") " ", name//" min/max (ice): ", &
            minval(v, mask=ice), "  ", maxval(v, mask=ice)
    end subroutine range_line

    subroutine range3(name, v, ice)
        character(len=*), intent(IN) :: name
        real(wp), intent(IN) :: v(:,:,:)
        logical,  intent(IN) :: ice(:,:)
        real(wp) :: vmin, vmax
        integer :: k
        vmin =  huge(1.0_wp)
        vmax = -huge(1.0_wp)
        do k = 1, size(v,3)
            vmin = min(vmin, minval(v(:,:,k), mask=ice))
            vmax = max(vmax, maxval(v(:,:,k), mask=ice))
        end do
        write(*,"(a,a,es12.4,a,es12.4)") " ", name//" min/max (ice): ", &
            vmin, "  ", vmax
    end subroutine range3

    function mean3(v, ice) result(m)
        real(wp), intent(IN) :: v(:,:,:)
        logical,  intent(IN) :: ice(:,:)
        real(wp) :: m, s
        integer :: k, n
        s = 0.0_wp
        n = 0
        do k = 1, size(v,3)
            s = s + sum(v(:,:,k), mask=ice)
            n = n + count(ice)
        end do
        m = s/real(n)
    end function mean3

end program test_domain
