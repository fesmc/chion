module chion_domain
    ! Driver-layer domain loader: assemble a STANDARDIZED monthly forcing set
    ! on an ice-sheet grid straight from the raw ~/models/ice_data datasets, so
    ! a driver can run chion on a real domain with no preprocessing step.
    !
    ! This lives ABOVE libchion.a, in the host/driver layer, exactly like
    ! libs/insol. It is where the per-domain knowledge lives -- which datasets
    ! carry temperature vs precipitation, what units they are in, what has to be
    ! regridded -- so libchion.a stays domain- and dataset-agnostic.
    !
    ! What "standardized" means. Every domain, whatever its native source
    ! (MAR for Greenland, RACMO/ERA-Interim for Antarctica, ...), is reduced to
    ! the SAME gridded monthly fields in SI units:
    !
    !     t2m   [K]              near-surface air temperature
    !     sf    [kg m-2 s-1]     snowfall rate
    !     rf    [kg m-2 s-1]     rainfall rate
    !     swd   [W m-2]          surface downward shortwave (from ERA5)
    !     tcc   [1]              total cloud cover        (from ERA5)
    !     S_toa [W m-2]          top-of-atmosphere daily insolation (computed)
    !     mask, z_srf, lon2D, lat2D                       static
    !
    ! plus the source model's own smb/melt/runoff, carried through untouched as
    ! a validation target. The driver builds its column list from `mask` and
    ! scatters these gridded fields onto columns; nothing here knows about
    ! columns or time stepping.
    !
    ! Adding a domain = one private loader subroutine (load_greenland below is
    ! the template) plus a case arm in chion_domain_load and in
    ! chion_domain_grid_def. The shared steps -- ERA5 regridding and the TOA
    ! insolation -- are done here once for every domain.
    !
    ! Regridding. Only ERA5 (global 0.25 lon-lat) is off-grid; the native
    ! ice-sheet fields already sit on the target grid. ERA5 is remapped
    ! CONSERVATIVELY (coords map_init method="con"), which preserves the
    ! domain-mean flux -- the quantity an ITM transmissivity is calibrated
    ! against -- rather than a point sample. The weight file is cached under
    ! maps/ and regenerated whenever it is absent.

    use, intrinsic :: iso_fortran_env, only : error_unit

    use coords,     only : grid_class, grid_init, map_class, map_init, map_field, dp
    use ncio,       only : nc_read, nc_size
    use insolation, only : calc_insol_day

    use chion_defs, only : wp, MV

    implicit none

    private

    type chion_domain_class
        character(len=256) :: domain
        character(len=256) :: grid_name

        integer :: nx, ny
        integer :: nmon
        integer :: nday_year

        ! Static geometry (nx,ny)
        real(wp), allocatable :: xc(:)          ! grid x axis [km]
        real(wp), allocatable :: yc(:)          ! grid y axis [km]
        real(wp), allocatable :: lon2D(:,:)     ! [deg E]
        real(wp), allocatable :: lat2D(:,:)     ! [deg N]
        real(wp), allocatable :: mask(:,:)      ! source land/ice mask (native scale)
        real(wp), allocatable :: z_srf(:,:)     ! surface elevation [m]

        ! Standardized monthly forcing (nx,ny,nmon), SI
        real(wp), allocatable :: t2m(:,:,:)     ! [K]
        real(wp), allocatable :: sf(:,:,:)      ! [kg m-2 s-1]
        real(wp), allocatable :: rf(:,:,:)      ! [kg m-2 s-1]
        real(wp), allocatable :: swd(:,:,:)     ! [W m-2]  (ERA5, regridded)
        real(wp), allocatable :: tcc(:,:,:)     ! [1]      (ERA5, regridded)

        ! Daily TOA insolation (nx,ny,nday_year) [W m-2]
        real(wp), allocatable :: S_toa(:,:,:)

        ! Source-model SMB, carried through for validation (nx,ny,nmon)
        real(wp), allocatable :: smb_ref(:,:,:)     ! [kg m-2 s-1]
        real(wp), allocatable :: melt_ref(:,:,:)    ! [kg m-2 s-1]
        real(wp), allocatable :: runoff_ref(:,:,:)  ! [kg m-2 s-1]
    end type

    ! mm water-equivalent per day -> kg m-2 s-1
    real(wp), parameter :: MMD_TO_KGM2S = 1.0_wp/86400.0_wp
    real(wp), parameter :: T_FREEZE     = 273.15_wp

    public :: chion_domain_class
    public :: chion_domain_load

contains

    subroutine chion_domain_load(dom, domain, grid_name, path_ice_data, &
                                 path_insol, nmon, nday_year, path_racmo)
        ! Assemble the standardized forcing for `domain` on `grid_name`, reading
        ! raw data from `path_ice_data` (the ~/models/ice_data root). Orbital
        ! tables for the insolation come from `path_insol` (libs/insol/input).
        !
        ! path_racmo is the root of the pre-built RACMO climatology (Antarctica
        ! only, e.g. ~/data/racmo). Required when domain=="antarctica", ignored
        ! otherwise -- the Greenland datasets all live under path_ice_data.

        implicit none

        type(chion_domain_class), intent(OUT) :: dom
        character(len=*), intent(IN) :: domain
        character(len=*), intent(IN) :: grid_name
        character(len=*), intent(IN) :: path_ice_data
        character(len=*), intent(IN) :: path_insol
        integer,          intent(IN) :: nmon
        integer,          intent(IN) :: nday_year
        character(len=*), intent(IN), optional :: path_racmo

        dom%domain    = trim(domain)
        dom%grid_name = trim(grid_name)
        dom%nmon      = nmon
        dom%nday_year = nday_year

        select case(trim(domain))

            case("greenland")
                call load_greenland(dom, grid_name, path_ice_data)

            case("antarctica")
                if (.not. present(path_racmo)) then
                    write(error_unit,*) "chion_domain_load:: error: domain 'antarctica' "// &
                                        "requires path_racmo (the RACMO climatology root)."
                    stop 1
                end if
                call load_antarctica(dom, grid_name, path_ice_data, path_racmo)

            case DEFAULT
                write(error_unit,*) "chion_domain_load:: error: unknown domain '"// &
                                    trim(domain)//"'."
                write(error_unit,*) "  supported: greenland, antarctica"
                stop 1

        end select

        ! Shared across all domains: TOA insolation from latitude.
        call domain_calc_insol(dom, path_insol)

        return

    end subroutine chion_domain_load

    ! ==================================================================
    ! Greenland: MARv3.11 monthly climatology + ERA5 radiation
    ! ==================================================================

    subroutine load_greenland(dom, grid_name, path_ice_data)

        implicit none

        type(chion_domain_class), intent(INOUT) :: dom
        character(len=*), intent(IN) :: grid_name
        character(len=*), intent(IN) :: path_ice_data

        ! Local variables
        character(len=512) :: file_mar
        integer :: nx, ny, nmon
        real(wp), allocatable :: tas(:,:,:)

        file_mar = trim(path_ice_data)//"/Greenland/"//trim(grid_name)//"/"// &
                   trim(grid_name)//"_MARv3.11-ERA_monmean_1961-1990.nc"

        call check_file(file_mar)

        nx   = nc_size(file_mar, "xc")
        ny   = nc_size(file_mar, "yc")
        nmon = nc_size(file_mar, "month")

        if (nmon .ne. dom%nmon) then
            write(error_unit,*) "load_greenland:: error: file has ", nmon, &
                                " months, expected ", dom%nmon
            stop 1
        end if

        dom%nx = nx
        dom%ny = ny

        call domain_alloc(dom)

        ! --- Static geometry ---------------------------------------------
        call nc_read(file_mar, "xc",    dom%xc)
        call nc_read(file_mar, "yc",    dom%yc)
        call nc_read(file_mar, "lon2D", dom%lon2D)
        call nc_read(file_mar, "lat2D", dom%lat2D)
        call nc_read(file_mar, "mask",  dom%mask)
        call nc_read(file_mar, "z_srf", dom%z_srf)

        ! --- Monthly forcing, converted to SI ----------------------------
        ! tas [degC] -> [K]
        allocate(tas(nx,ny,nmon))
        call nc_read(file_mar, "tas", tas)
        dom%t2m = tas + T_FREEZE
        deallocate(tas)

        ! sf, rf [mm d-1] -> [kg m-2 s-1]
        call nc_read(file_mar, "sf", dom%sf)
        call nc_read(file_mar, "rf", dom%rf)
        dom%sf = dom%sf * MMD_TO_KGM2S
        dom%rf = dom%rf * MMD_TO_KGM2S

        ! Reference SMB terms [mm d-1] -> [kg m-2 s-1]
        call nc_read(file_mar, "smb",    dom%smb_ref)
        call nc_read(file_mar, "melt",   dom%melt_ref)
        call nc_read(file_mar, "runoff", dom%runoff_ref)
        dom%smb_ref    = dom%smb_ref    * MMD_TO_KGM2S
        dom%melt_ref   = dom%melt_ref   * MMD_TO_KGM2S
        dom%runoff_ref = dom%runoff_ref * MMD_TO_KGM2S

        ! --- ERA5 radiation, conservatively regridded --------------------
        call regrid_era5_monthly(dom, grid_name, path_ice_data, &
                                 "mean_surface_downward_short_wave_radiation_flux", &
                                 "msdwswrf", dom%swd)
        call regrid_era5_monthly(dom, grid_name, path_ice_data, &
                                 "total_cloud_cover", "tcc", dom%tcc)

        return

    end subroutine load_greenland

    ! ==================================================================
    ! Antarctica: RACMO2.4/ANT-12 monthly climatology on the target ANT grid,
    ! with geometry from BedMachine.
    ! ==================================================================

    subroutine load_antarctica(dom, grid_name, path_ice_data, path_racmo)
        ! Forcing is the RACMO2.4 climatology already conservatively regridded
        ! onto the target ANT grid by scripts/build_ant12_climatology.sh (cdo
        ! remapcon from the native CORDEX ANT-12 rotated-pole grid). The regrid
        ! is done in preprocessing rather than here because coords cannot yet
        ! remap from a rotated-pole source -- its conservative map returns zero
        ! overlaps (fesmc/fesm-utils#8).
        ! Geometry (mask, surface elevation, lon/lat) comes from BedMachine, on
        ! the same grid. The CORDEX set carries no snow/rain split and no
        ! reference SMB, so precip is split by a freezing threshold and the
        ! *_ref validation fields are left at zero.

        implicit none

        type(chion_domain_class), intent(INOUT) :: dom
        character(len=*), intent(IN) :: grid_name
        character(len=*), intent(IN) :: path_ice_data
        character(len=*), intent(IN) :: path_racmo

        ! Local variables
        character(len=512) :: file_clim, file_topo
        integer :: nmon, nx, ny
        real(wp), allocatable :: mask_cat(:,:)
        type(grid_class) :: grid_tgt

        ! --- Climatology, already on the target grid ---------------------
        file_clim = trim(path_racmo)//"/clim/"//trim(grid_name)// &
                    "_RACMO24P_monclim_1981-2010.nc"
        call check_file(file_clim)

        nmon = nc_size(file_clim, "time")
        if (nmon .ne. dom%nmon) then
            write(error_unit,*) "load_antarctica:: error: climatology has ", nmon, &
                                " months, expected ", dom%nmon
            stop 1
        end if

        ! --- Geometry from BedMachine, defines the target grid -----------
        file_topo = trim(path_ice_data)//"/Antarctica/"//trim(grid_name)//"/"// &
                    trim(grid_name)//"_TOPO-BedMachine.nc"
        call check_file(file_topo)

        nx = nc_size(file_topo, "xc")
        ny = nc_size(file_topo, "yc")

        if (nc_size(file_clim,"xc") .ne. nx .or. nc_size(file_clim,"yc") .ne. ny) then
            write(error_unit,*) "load_antarctica:: error: climatology grid (", &
                nc_size(file_clim,"xc"), "x", nc_size(file_clim,"yc"), &
                ") does not match BedMachine grid (", nx, "x", ny, ") for "// &
                trim(grid_name)//". Re-run scripts/build_ant12_climatology.sh."
            stop 1
        end if

        dom%nx = nx
        dom%ny = ny

        call domain_alloc(dom)

        call nc_read(file_topo, "xc",    dom%xc)
        call nc_read(file_topo, "yc",    dom%yc)
        call nc_read(file_topo, "z_srf", dom%z_srf)

        ! lon2D/lat2D from the grid definition, not the TOPO file: the ANT-8KM
        ! BedMachine file ships without them. lat2D feeds the TOA insolation.
        call chion_domain_grid_def(grid_tgt, grid_name)
        dom%lon2D = real(grid_tgt%lon, wp)
        dom%lat2D = real(grid_tgt%lat, wp)

        ! BedMachine mask is categorical (0 ocean, 1 ice-free land, 2 grounded
        ! ice, 3 floating ice, 4 lake vostok). Reduce it to the Greenland-style
        ! percentage convention so the driver's mask_threshold works uniformly:
        ! every ice category (grounded + floating shelves + lake vostok) -> 100.
        allocate(mask_cat(nx,ny))
        call nc_read(file_topo, "mask", mask_cat)
        dom%mask = 0.0_wp
        where (mask_cat .gt. 1.5_wp) dom%mask = 100.0_wp
        deallocate(mask_cat)

        ! --- Standardized forcing (already on the target grid) -----------
        ! Total precip is read into rf first, then split into sf/rf below.
        call nc_read(file_clim, "t2m", dom%t2m)
        call nc_read(file_clim, "pr",  dom%rf)
        call nc_read(file_clim, "swd", dom%swd)
        call nc_read(file_clim, "tcc", dom%tcc)

        ! --- Snow/rain split by a freezing threshold ---------------------
        ! CORDEX gives only total precip. Antarctica is almost entirely below
        ! freezing, so precip defaults to snowfall; the few months/cells above
        ! 0 C become rainfall. dom%rf currently holds the total precip.
        dom%sf = 0.0_wp
        where (dom%t2m .le. T_FREEZE)
            dom%sf = dom%rf
            dom%rf = 0.0_wp
        end where

        ! smb_ref/melt_ref/runoff_ref: no reference SMB in this CORDEX set,
        ! left at the zero from domain_alloc.

        return

    end subroutine load_antarctica

    ! ==================================================================
    ! Shared: ERA5 conservative regridding
    ! ==================================================================

    subroutine regrid_era5_monthly(dom, grid_name, path_ice_data, era5_name, &
                                   var_name, var_out)
        ! Read a global ERA5 monthly-climatology field and remap it onto the
        ! target grid conservatively. The (source, target) grids and the weight
        ! map are rebuilt on each call; that is once per variable per run and
        ! the map itself is cached to maps/ by coords.

        implicit none

        type(chion_domain_class), intent(IN)  :: dom
        character(len=*), intent(IN)  :: grid_name
        character(len=*), intent(IN)  :: path_ice_data
        character(len=*), intent(IN)  :: era5_name
        character(len=*), intent(IN)  :: var_name
        real(wp),         intent(OUT) :: var_out(:,:,:)   ! (nx,ny,nmon)

        ! Local variables
        character(len=512) :: file_era5
        type(grid_class)   :: grid_src, grid_tgt
        integer :: nlon, nlat, nmon
        real(dp), allocatable :: lon(:), lat(:)
        real(dp), allocatable :: src(:,:,:)

        file_era5 = trim(path_ice_data)//"/ERA5/clim/era5_monthly-single-levels_"// &
                    trim(era5_name)//"_1961-1990.nc"
        call check_file(file_era5)

        nlon = nc_size(file_era5, "longitude")
        nlat = nc_size(file_era5, "latitude")
        nmon = nc_size(file_era5, "time")

        allocate(lon(nlon), lat(nlat))
        call nc_read(file_era5, "longitude", lon)
        call nc_read(file_era5, "latitude",  lat)

        ! ERA5 field arrives (time,lat,lon) in the file => (lon,lat,time) here.
        allocate(src(nlon,nlat,nmon))
        call nc_read(file_era5, var_name, src)

        ! Source grid: global regular lon-lat. lon180 lets coords reconcile the
        ! 0..360 ERA5 convention with the -180..180 target.
        call grid_init(grid_src, name="ERA5-025deg", mtype="latlon", &
                       units="degrees", lon180=.TRUE., x=lon, y=lat)

        call chion_domain_grid_def(grid_tgt, grid_name)

        call regrid_monthly_con(grid_src, grid_tgt, trim(var_name), src, var_out)

        deallocate(lon, lat, src)

        return

    end subroutine regrid_era5_monthly

    subroutine regrid_monthly_con(grid_src, grid_tgt, var_name, src, var_out)
        ! Conservatively remap every month of a (src_nx,src_ny,nmon) field from
        ! grid_src onto grid_tgt, writing (tgt_nx,tgt_ny,nmon). The weight map is
        ! cached under maps/ (regenerated if absent) keyed by the grid names, so
        ! each (source,target) pair is built once and reused across variables.

        implicit none

        type(grid_class), intent(INOUT) :: grid_src, grid_tgt
        character(len=*), intent(IN)    :: var_name
        real(dp),         intent(IN)    :: src(:,:,:)      ! (src_nx,src_ny,nmon)
        real(wp),         intent(OUT)   :: var_out(:,:,:)  ! (tgt_nx,tgt_ny,nmon)

        type(map_class) :: map
        integer :: nmon, m
        real(dp), allocatable :: src_m(:,:), tgt_m(:,:)

        nmon = size(src,3)

        call map_init(map, grid_src, grid_tgt, method="con", gen="coords", &
                      fldr="maps", load=.TRUE.)

        allocate(src_m(size(src,1),size(src,2)))
        allocate(tgt_m(size(var_out,1),size(var_out,2)))
        do m = 1, nmon
            src_m = src(:,:,m)
            call map_field(map, trim(var_name), src_m, tgt_m, stat="mean", &
                           missing_value=real(MV,dp))
            var_out(:,:,m) = real(tgt_m, wp)
        end do

        deallocate(src_m, tgt_m)

        return

    end subroutine regrid_monthly_con

    ! ==================================================================
    ! Shared: top-of-atmosphere daily insolation
    ! ==================================================================

    subroutine domain_calc_insol(dom, path_insol)
        ! Daily-mean TOA insolation for every day of the (nday_year) year, on
        ! the domain grid, from the Laskar orbital solution at present day.

        implicit none

        type(chion_domain_class), intent(INOUT) :: dom
        character(len=*), intent(IN) :: path_insol

        ! Local variables
        real(dp), allocatable :: lats(:,:), insol(:,:)
        integer :: day

        allocate(lats(dom%nx,dom%ny), insol(dom%nx,dom%ny))
        lats = real(dom%lat2D, dp)

        do day = 1, dom%nday_year
            insol = calc_insol_day(day, lats, 0.0_dp, day_year=dom%nday_year, &
                                   fldr=trim(path_insol))
            dom%S_toa(:,:,day) = real(insol, wp)
        end do

        deallocate(lats, insol)

        return

    end subroutine domain_calc_insol

    ! ==================================================================
    ! Per-domain grid definitions (mirrors gridding/src/control.f90)
    ! ==================================================================

    subroutine chion_domain_grid_def(grid, grid_name)
        ! Projection definition of a supported ice-sheet grid. These reproduce
        ! the grids the ice_data files were written on
        ! (gridding/src/control.f90), so a regrid target lands exactly on the
        ! native cells.

        implicit none

        type(grid_class), intent(OUT) :: grid
        character(len=*), intent(IN)  :: grid_name

        select case(trim(grid_name))

            case("GRL-16KM")
                call grid_init(grid, name="GRL-16KM", mtype="polar_stereographic", &
                               units="kilometers", lon180=.TRUE., &
                               x0=-720.0_dp, dx=16.0_dp, nx=106, &
                               y0=-3450.0_dp, dy=16.0_dp, ny=181, &
                               lambda=-45.0_dp, phi=70.0_dp)

            case("GRL-32KM")
                call grid_init(grid, name="GRL-32KM", mtype="polar_stereographic", &
                               units="kilometers", lon180=.TRUE., &
                               x0=-720.0_dp, dx=32.0_dp, nx=54, &
                               y0=-3450.0_dp, dy=32.0_dp, ny=91, &
                               lambda=-45.0_dp, phi=70.0_dp)

            case("GRL-8KM")
                call grid_init(grid, name="GRL-8KM", mtype="polar_stereographic", &
                               units="kilometers", lon180=.TRUE., &
                               x0=-720.0_dp, dx=8.0_dp, nx=211, &
                               y0=-3450.0_dp, dy=8.0_dp, ny=361, &
                               lambda=-45.0_dp, phi=70.0_dp)

            ! Antarctic grids: south polar stereographic, central meridian 0,
            ! standard parallel 71 S (the standard ANT-XXKM / BedMachine grids).
            ! Used by load_antarctica only to derive lon2D/lat2D (the RACMO
            ! forcing is pre-regridded), so the ANT-8KM BedMachine file -- which
            ! ships without lon2D/lat2D -- still gets consistent coordinates.
            case("ANT-32KM")
                call grid_init(grid, name="ANT-32KM", mtype="polar_stereographic", &
                               units="kilometers", lon180=.TRUE., &
                               x0=-3040.0_dp, dx=32.0_dp, nx=191, &
                               y0=-3040.0_dp, dy=32.0_dp, ny=191, &
                               lambda=0.0_dp, phi=-71.0_dp)

            case("ANT-16KM")
                call grid_init(grid, name="ANT-16KM", mtype="polar_stereographic", &
                               units="kilometers", lon180=.TRUE., &
                               x0=-3040.0_dp, dx=16.0_dp, nx=381, &
                               y0=-3040.0_dp, dy=16.0_dp, ny=381, &
                               lambda=0.0_dp, phi=-71.0_dp)

            case("ANT-8KM")
                call grid_init(grid, name="ANT-8KM", mtype="polar_stereographic", &
                               units="kilometers", lon180=.TRUE., &
                               x0=-3040.0_dp, dx=8.0_dp, nx=761, &
                               y0=-3040.0_dp, dy=8.0_dp, ny=761, &
                               lambda=0.0_dp, phi=-71.0_dp)

            case DEFAULT
                write(error_unit,*) "chion_domain_grid_def:: error: unknown grid '"// &
                                    trim(grid_name)//"'."
                stop 1

        end select

        return

    end subroutine chion_domain_grid_def

    ! ==================================================================
    ! Internal
    ! ==================================================================

    subroutine domain_alloc(dom)

        implicit none

        type(chion_domain_class), intent(INOUT) :: dom
        integer :: nx, ny, nmon, nday

        nx   = dom%nx
        ny   = dom%ny
        nmon = dom%nmon
        nday = dom%nday_year

        call realloc3(dom%t2m,        nx, ny, nmon)
        call realloc3(dom%sf,         nx, ny, nmon)
        call realloc3(dom%rf,         nx, ny, nmon)
        call realloc3(dom%swd,        nx, ny, nmon)
        call realloc3(dom%tcc,        nx, ny, nmon)
        call realloc3(dom%smb_ref,    nx, ny, nmon)
        call realloc3(dom%melt_ref,   nx, ny, nmon)
        call realloc3(dom%runoff_ref, nx, ny, nmon)
        call realloc3(dom%S_toa,      nx, ny, nday)

        if (allocated(dom%xc))    deallocate(dom%xc)
        if (allocated(dom%yc))    deallocate(dom%yc)
        if (allocated(dom%lon2D)) deallocate(dom%lon2D)
        if (allocated(dom%lat2D)) deallocate(dom%lat2D)
        if (allocated(dom%mask))  deallocate(dom%mask)
        if (allocated(dom%z_srf)) deallocate(dom%z_srf)
        allocate(dom%xc(nx), dom%yc(ny))
        allocate(dom%lon2D(nx,ny), dom%lat2D(nx,ny))
        allocate(dom%mask(nx,ny),  dom%z_srf(nx,ny))

        return

    end subroutine domain_alloc

    subroutine realloc3(a, n1, n2, n3)

        implicit none

        real(wp), allocatable, intent(INOUT) :: a(:,:,:)
        integer, intent(IN) :: n1, n2, n3

        if (allocated(a)) deallocate(a)
        allocate(a(n1,n2,n3))
        a = 0.0_wp

        return

    end subroutine realloc3

    subroutine check_file(filename)

        implicit none

        character(len=*), intent(IN) :: filename
        logical :: exists

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            write(error_unit,*) "chion_domain:: error: file not found:"
            write(error_unit,*) "  "//trim(filename)
            stop 1
        end if

        return

    end subroutine check_file

end module chion_domain
