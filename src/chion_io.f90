module chion_io
    ! WP14 -- NetCDF output and restart.
    !
    ! Four public routines and nothing else:
    !
    !     call chion_write_init(chn,filename,time_init,units)   ! once
    !     call chion_write_step(chn,filename,time)              ! every output step
    !
    !     call chion_restart_write(chn,filename,time)
    !     call chion_restart_read (chn,filename,time)           ! time is OUT
    !
    ! -----------------------------------------------------------------------
    ! Metadata lives in markdown tables, not in this file
    ! -----------------------------------------------------------------------
    ! One table per model:
    !     input/chion-variables-bessi.md
    !     input/chion-variables-pdd.md
    !     input/chion-variables-itm.md
    ! loaded once per model by load_var_io_table (fesm-utils variable_io) and
    ! cached in the module, exactly as yelmo does (yelmo/src/yelmo_io.f90).
    ! Columns are `id | variable | dimensions | units | long_name`, and the
    ! writer takes units and long_name from the table -- never from code -- so
    ! that the table is the single place a name or a unit can be changed.
    !
    ! Variable names and units are Chion.jl's (src/io.jl NETCDF_METADATA), so
    ! chion and Chion.jl output files are directly comparable. See the tables
    ! for the handful of chion-only additions.
    !
    ! -----------------------------------------------------------------------
    ! The `dimensions` column is LOGICAL
    ! -----------------------------------------------------------------------
    ! chion's state is a packed LIST of ncol independent columns, with an
    ! OPTIONAL 2-D mapping (chion_set_grid). The tables therefore say `column`,
    ! and this module expands it:
    !
    !     grd%has_spatial = .TRUE.    column -> xc, yc      (scatter, see below)
    !     grd%has_spatial = .FALSE.   column -> column      (bare list)
    !
    ! `layer` is used unchanged, with extent Ntot (BESSI only). `time` is
    ! appended by the writer and is ALWAYS the unlimited dimension.
    !
    ! The scatter follows Chion.jl scatter_to_grid (src/io.jl:45): every column
    ! i is placed at grid cell (js(i),is(i)) and every cell that no column maps
    ! to keeps the missing value. Chion.jl fills with NaN; chion fills with
    ! chion_defs' MV = -9999 and declares it as the NetCDF _FillValue, because
    ! a NaN in an sp field is indistinguishable from an FPE-trapped one under
    ! the project's debug flags (docs/porting_notes.md D12).
    !
    ! In-memory the scattered array is (nx,ny) and is written with
    ! dim1="xc",dim2="yc", so the file dimension order is (time,yc,xc) --
    ! CF-standard, and the same convention as yelmo. Chion.jl declares
    ! ("t","x","y") which lands in the file as (y,x,t). The dimension ORDER
    ! therefore differs from Chion.jl and cannot be made to match, because
    ! chion's `time` is unlimited and netCDF requires the unlimited dimension
    ! to be the slowest-varying one. Names, units and long names do match,
    ! which is what WP16 compares.
    !
    ! -----------------------------------------------------------------------
    ! Restart
    ! -----------------------------------------------------------------------
    ! The restart file is a DIFFERENT file from the output file and has a
    ! different job: reproduce the prognostic state exactly, not present it.
    ! Three consequences:
    !
    !   1. It NEVER scatters. Even when a spatial mapping is attached, the
    !      restart is written on the bare `column` dimension, because a scatter
    !      is lossy the moment two columns share a grid cell and because the
    !      column ORDER is what the state arrays are indexed by.
    !   2. It writes every wp_acc accumulator as NF90_DOUBLE. Demoting them to
    !      sp would silently undo the entire point of wp_acc
    !      (docs/porting_notes.md D1) at every restart.
    !   3. It enumerates the state field by field rather than looping a table.
    !      A field missed here is a silent wrong answer after a restart, so the
    !      list is written out where it can be read against the state type
    !      declarations, and tests/test_io.f90 checks it field by field.
    !
    ! It also carries chn%smb_cum_prev and chn%dt_last (docs/porting_notes.md
    ! D13). Without those two, the first chion_get_smb after a restart reports
    ! the whole cumulative accumulator as if it were one step's flux.
    !
    ! The model name, Ntot and ncol are stored as global attributes and are
    ! CHECKED on read. Loading a BESSI restart into a PDD run, or a Ntot=15
    ! restart into a Ntot=10 run, is refused with an explicit message rather
    ! than silently producing garbage or a shape error from netCDF.

    use ncio
    use variable_io, only : var_io_type, load_var_io_table

    use chion_defs, only : wp, wp_acc, MV, io_unit_err, &
                           chion_grid_class, chion_check_file, &
                           chion_grid_set_active

    use chion_api,  only : chion_class, chion_get_smb

    use snow_diagnostics, only : summarize_domain_state

    implicit none

    private

    ! Where the metadata tables live. Hard-coded, as in yelmo and as
    ! chion_api's def_file is: they are part of the package, not part of a run
    ! configuration.
    character(len=*), parameter :: table_bessi = "input/chion-variables-bessi.md"
    character(len=*), parameter :: table_pdd   = "input/chion-variables-pdd.md"
    character(len=*), parameter :: table_itm   = "input/chion-variables-itm.md"

    ! Cached tables. Loaded on first use for the model actually running, so a
    ! PDD run never touches the BESSI table. `save` is safe here because all IO
    ! is serial -- no chion_io routine is ever called from inside the OpenMP
    ! column loop.
    type(var_io_type), allocatable, save :: var_table(:)
    character(len=56),              save :: var_table_model = "none"

    public :: chion_write_init
    public :: chion_write_step
    public :: chion_restart_write
    public :: chion_restart_read

contains

    ! =====================================================================
    ! Output
    ! =====================================================================

    subroutine chion_write_init(chn,filename,time_init,units)
        ! Create the output file: dimensions, coordinate variables and the
        ! static fields. Computes nothing and writes no state -- the first
        ! record is written by the first chion_write_step.
        !
        ! units : the time units string, e.g. "days" or "years". Optional,
        !         default "days", which is the unit chion_update takes.

        implicit none

        type(chion_class),          intent(IN) :: chn
        character(len=*),           intent(IN) :: filename
        real(wp),                   intent(IN) :: time_init
        character(len=*), optional, intent(IN) :: units

        ! Local variables
        integer :: i, nx, ny
        character(len=56) :: time_units
        real(wp), allocatable :: mask2D(:,:)
        real(wp), allocatable :: col_i(:)

        time_units = "days"
        if (present(units)) time_units = trim(units)

        call nc_create(filename)

        ! --- Column / spatial dimensions ------------------------------

        if (chn%grd%has_spatial) then

            nx = size(chn%grd%x)
            ny = size(chn%grd%y)

            call nc_write_dim(filename,"xc",x=chn%grd%x,units="km")
            call nc_write_dim(filename,"yc",x=chn%grd%y,units="km")

        else

            call nc_write_dim(filename,"column",x=1,dx=1,nx=chn%grd%ncol,units="1")

        end if

        ! --- Layer dimension (BESSI only) -----------------------------
        ! Extent is the CONFIGURED Ntot. Note the snow-depth cap uses a
        ! hard-coded reference count of 15 regardless (docs/PLAN.md section 5
        ! item 11); this dimension is the array extent, not that reference.

        if (trim(chn%par%model) .eq. "bessi") then
            call nc_write_dim(filename,"layer",x=1,dx=1,nx=chn%bsi%now%Ntot,units="1")
        end if

        ! --- Time, unlimited ------------------------------------------

        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1, &
                          units=trim(time_units),unlimited=.TRUE.)

        ! --- Static fields --------------------------------------------

        if (chn%grd%has_spatial) then

            ! grd%mask is (ny,nx); the file convention is (nx,ny) written as
            ! dim1=xc,dim2=yc. Transposed here rather than at every call site.
            allocate(mask2D(nx,ny))
            mask2D = transpose(chn%grd%mask)

            call nc_write(filename,"mask",mask2D,dim1="xc",dim2="yc", &
                          units="1",long_name="Domain mask",grid_mapping="")

            ! The column -> grid mapping itself, scattered onto the grid. This
            ! is what makes an output file self-describing: a reader can
            ! recover which column produced which cell without the parameter
            ! file.
            allocate(col_i(chn%grd%ncol))
            do i = 1, chn%grd%ncol
                col_i(i) = real(i,wp)
            end do

            call chion_write_static_col(filename,"column_index",col_i,chn%grd, &
                                        "1","Index of the column mapped to this cell")

            deallocate(col_i)
            deallocate(mask2D)

        end if

        return

    end subroutine chion_write_init

    subroutine chion_write_step(chn,filename,time)
        ! Append one record. Opens the file ONCE, resolves the record index
        ! with nc_time_index, writes every variable in the model's table, and
        ! closes once -- rather than the ncio default of open/close per
        ! variable, which for a 20-variable BESSI table is 40 file operations
        ! per output step.

        implicit none

        type(chion_class), intent(IN) :: chn
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        ! Local variables
        integer :: q, n, ncid, ncol
        real(wp), allocatable :: smb(:)
        real(wp), allocatable :: thickness(:), wet_mass(:), bulk_density(:), liquid_water(:)

        call chion_load_table(chn%par%model)

        ncol = chn%grd%ncol

        ! The model-agnostic ice-facing flux for the step just completed.
        ! Computed here rather than stored, because chn is intent(IN).
        allocate(smb(ncol))
        call chion_get_smb(chn,smb)

        ! BESSI's four diagnostics are NOT updated by bessi_column_step -- in
        ! Chion.jl they are produced by summarize_domain_state on the way to
        ! output, and chion keeps that split. So they are recomputed here from
        ! the layer state. The values in chn%bsi%now%thickness etc. are
        ! whatever the last caller of summarize_domain_state left there and
        ! are deliberately not trusted.
        if (trim(chn%par%model) .eq. "bessi") then
            allocate(thickness(ncol))
            allocate(wet_mass(ncol))
            allocate(bulk_density(ncol))
            allocate(liquid_water(ncol))
            call summarize_domain_state(chn%bsi%now%mass,chn%bsi%now%mass_w, &
                                        chn%bsi%now%density,chn%bsi%now%n_lay, &
                                        thickness,wet_mass,bulk_density,liquid_water)
        end if

        call nc_open(filename,ncid,writable=.TRUE.)

        n = nc_time_index(filename,"time",time,ncid)
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],grid_mapping="",ncid=ncid)

        select case(trim(chn%par%model))

            case("bessi")
                do q = 1, size(var_table)
                    call chion_write_var_bessi(filename,var_table(q),chn,smb, &
                                               thickness,wet_mass,bulk_density,liquid_water, &
                                               n,ncid)
                end do

            case("pdd")
                do q = 1, size(var_table)
                    call chion_write_var_pdd(filename,var_table(q),chn,smb,n,ncid)
                end do

            case("itm")
                do q = 1, size(var_table)
                    call chion_write_var_itm(filename,var_table(q),chn,smb,n,ncid)
                end do

            case DEFAULT
                call chion_io_model_error("chion_write_step",chn%par%model)

        end select

        call nc_close(ncid)

        return

    end subroutine chion_write_step

    ! ---------------------------------------------------------------------
    ! Per-model variable dispatch. One case per table entry, yelmo style
    ! (yelmo_write_var_io_ytopo). A table entry with no case here is a hard
    ! error rather than a silent omission.
    ! ---------------------------------------------------------------------

    subroutine chion_write_var_bessi(filename,v,chn,smb,thickness,wet_mass, &
                                     bulk_density,liquid_water,n,ncid)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        type(chion_class), intent(IN) :: chn
        real(wp),          intent(IN) :: smb(:)
        real(wp),          intent(IN) :: thickness(:)
        real(wp),          intent(IN) :: wet_mass(:)
        real(wp),          intent(IN) :: bulk_density(:)
        real(wp),          intent(IN) :: liquid_water(:)
        integer,           intent(IN) :: n
        integer,           intent(IN) :: ncid

        select case(trim(v%varname))

            case("thickness")
                call chion_write_col(filename,v,thickness,chn%grd,n,ncid)
            case("wet_mass")
                call chion_write_col(filename,v,wet_mass,chn%grd,n,ncid)
            case("bulk_density")
                call chion_write_col(filename,v,bulk_density,chn%grd,n,ncid)
            case("liquid_water")
                call chion_write_col(filename,v,liquid_water,chn%grd,n,ncid)

            case("mass_base")
                call chion_write_col(filename,v,real(chn%bsi%now%mass_base,wp),chn%grd,n,ncid)
            case("smb_ice")
                call chion_write_col(filename,v,real(chn%bsi%now%smb_ice,wp),chn%grd,n,ncid)
            case("runoff")
                call chion_write_col(filename,v,real(chn%bsi%now%runoff,wp),chn%grd,n,ncid)
            case("melt")
                call chion_write_col(filename,v,real(chn%bsi%now%melt,wp),chn%grd,n,ncid)
            case("refreezing")
                call chion_write_col(filename,v,real(chn%bsi%now%refreezing,wp),chn%grd,n,ncid)
            case("sublimation")
                call chion_write_col(filename,v,real(chn%bsi%now%sublimation,wp),chn%grd,n,ncid)
            case("vapor_mass")
                call chion_write_col(filename,v,real(chn%bsi%now%vapor_mass,wp),chn%grd,n,ncid)
            case("latent_heat_flux_sum")
                call chion_write_col(filename,v,real(chn%bsi%now%latent_heat_flux_sum,wp),chn%grd,n,ncid)

            case("Tsrf")
                call chion_write_col(filename,v,chn%bsi%now%t_srf,chn%grd,n,ncid)
            case("albedo")
                call chion_write_col(filename,v,chn%bsi%now%albedo,chn%grd,n,ncid)
            case("N")
                ! Written as a float so that unmapped grid cells can carry MV.
                ! The restart writes n_lay as an integer.
                call chion_write_col(filename,v,real(chn%bsi%now%n_lay,wp),chn%grd,n,ncid)
            case("smb")
                call chion_write_col(filename,v,smb,chn%grd,n,ncid)

            case("mass")
                call chion_write_lay(filename,v,chn%bsi%now%mass,chn%grd,n,ncid)
            case("mass_w")
                call chion_write_lay(filename,v,chn%bsi%now%mass_w,chn%grd,n,ncid)
            case("density")
                call chion_write_lay(filename,v,chn%bsi%now%density,chn%grd,n,ncid)
            case("temperature")
                call chion_write_lay(filename,v,chn%bsi%now%temperature,chn%grd,n,ncid)

            case DEFAULT
                call chion_io_var_error("chion_write_var_bessi",v%varname,table_bessi)

        end select

        return

    end subroutine chion_write_var_bessi

    subroutine chion_write_var_pdd(filename,v,chn,smb,n,ncid)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        type(chion_class), intent(IN) :: chn
        real(wp),          intent(IN) :: smb(:)
        integer,           intent(IN) :: n
        integer,           intent(IN) :: ncid

        select case(trim(v%varname))

            case("snowpack_swe")
                call chion_write_col(filename,v,chn%pdd%now%snowpack_swe,chn%grd,n,ncid)
            case("smb_ice")
                call chion_write_col(filename,v,real(chn%pdd%now%smb_ice,wp),chn%grd,n,ncid)
            case("runoff")
                call chion_write_col(filename,v,real(chn%pdd%now%runoff,wp),chn%grd,n,ncid)
            case("pdd_sum")
                call chion_write_col(filename,v,real(chn%pdd%now%pdd_sum,wp),chn%grd,n,ncid)
            case("smb")
                call chion_write_col(filename,v,smb,chn%grd,n,ncid)

            case DEFAULT
                call chion_io_var_error("chion_write_var_pdd",v%varname,table_pdd)

        end select

        return

    end subroutine chion_write_var_pdd

    subroutine chion_write_var_itm(filename,v,chn,smb,n,ncid)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        type(chion_class), intent(IN) :: chn
        real(wp),          intent(IN) :: smb(:)
        integer,           intent(IN) :: n
        integer,           intent(IN) :: ncid

        select case(trim(v%varname))

            case("H_snow")
                call chion_write_col(filename,v,chn%itm%now%H_snow,chn%grd,n,ncid)
            case("albedo")
                call chion_write_col(filename,v,chn%itm%now%alb_s,chn%grd,n,ncid)
            case("Tsrf")
                call chion_write_col(filename,v,chn%itm%now%tsrf,chn%grd,n,ncid)

            case("smb_ice")
                call chion_write_col(filename,v,real(chn%itm%now%smbi_cum,wp),chn%grd,n,ncid)
            case("runoff")
                call chion_write_col(filename,v,real(chn%itm%now%runoff_cum,wp),chn%grd,n,ncid)
            case("melt")
                call chion_write_col(filename,v,real(chn%itm%now%melt_cum,wp),chn%grd,n,ncid)
            case("refreezing")
                call chion_write_col(filename,v,real(chn%itm%now%refrz_cum,wp),chn%grd,n,ncid)
            case("smb_total")
                call chion_write_col(filename,v,real(chn%itm%now%smb_cum,wp),chn%grd,n,ncid)
            case("smb")
                call chion_write_col(filename,v,smb,chn%grd,n,ncid)

            case DEFAULT
                call chion_io_var_error("chion_write_var_itm",v%varname,table_itm)

        end select

        return

    end subroutine chion_write_var_itm

    ! ---------------------------------------------------------------------
    ! The two writers that know about the column/grid duality
    ! ---------------------------------------------------------------------

    subroutine chion_write_col(filename,v,dat,grd,n,ncid)
        ! One (ncol) field, one time record.

        implicit none

        character(len=*),       intent(IN) :: filename
        type(var_io_type),      intent(IN) :: v
        real(wp),               intent(IN) :: dat(:)
        type(chion_grid_class), intent(IN) :: grd
        integer,                intent(IN) :: n
        integer,                intent(IN) :: ncid

        ! Local variables
        real(wp), allocatable :: dat2D(:,:)

        if (grd%has_spatial) then

            call chion_scatter_to_grid(dat2D,dat,grd)

            call nc_write(filename,trim(v%varname),dat2D, &
                          dim1="xc",dim2="yc",dim3="time",start=[1,1,n], &
                          units=trim(v%units),long_name=trim(v%long_name), &
                          missing_value=MV,grid_mapping="",ncid=ncid)

            deallocate(dat2D)

        else

            call nc_write(filename,trim(v%varname),dat, &
                          dim1="column",dim2="time",start=[1,n], &
                          units=trim(v%units),long_name=trim(v%long_name), &
                          missing_value=MV,grid_mapping="",ncid=ncid)

        end if

        return

    end subroutine chion_write_col

    subroutine chion_write_lay(filename,v,dat,grd,n,ncid)
        ! One (Ntot,ncol) layered field, one time record.
        !
        ! Layer 1 is the SURFACE and the index increases downward, matching
        ! Chion.jl. The `layer` coordinate is 1..Ntot and carries no depth
        ! information -- layer thickness is mass/density and varies per column
        ! and per step, so there is no common vertical axis to write.

        implicit none

        character(len=*),       intent(IN) :: filename
        type(var_io_type),      intent(IN) :: v
        real(wp),               intent(IN) :: dat(:,:)
        type(chion_grid_class), intent(IN) :: grd
        integer,                intent(IN) :: n
        integer,                intent(IN) :: ncid

        ! Local variables
        integer :: k, nx, ny, Ntot
        real(wp), allocatable :: dat3D(:,:,:)
        real(wp), allocatable :: dat2D(:,:)

        Ntot = size(dat,1)

        if (grd%has_spatial) then

            nx = size(grd%x)
            ny = size(grd%y)

            allocate(dat3D(nx,ny,Ntot))

            do k = 1, Ntot
                call chion_scatter_to_grid(dat2D,dat(k,:),grd)
                dat3D(:,:,k) = dat2D
                deallocate(dat2D)
            end do

            call nc_write(filename,trim(v%varname),dat3D, &
                          dim1="xc",dim2="yc",dim3="layer",dim4="time",start=[1,1,1,n], &
                          units=trim(v%units),long_name=trim(v%long_name), &
                          missing_value=MV,grid_mapping="",ncid=ncid)

            deallocate(dat3D)

        else

            call nc_write(filename,trim(v%varname),dat, &
                          dim1="layer",dim2="column",dim3="time",start=[1,1,n], &
                          units=trim(v%units),long_name=trim(v%long_name), &
                          missing_value=MV,grid_mapping="",ncid=ncid)

        end if

        return

    end subroutine chion_write_lay

    subroutine chion_write_static_col(filename,varname,dat,grd,units,long_name)
        ! A time-independent scattered field, for chion_write_init.

        implicit none

        character(len=*),       intent(IN) :: filename
        character(len=*),       intent(IN) :: varname
        real(wp),               intent(IN) :: dat(:)
        type(chion_grid_class), intent(IN) :: grd
        character(len=*),       intent(IN) :: units
        character(len=*),       intent(IN) :: long_name

        ! Local variables
        real(wp), allocatable :: dat2D(:,:)

        call chion_scatter_to_grid(dat2D,dat,grd)

        call nc_write(filename,trim(varname),dat2D,dim1="xc",dim2="yc", &
                      units=trim(units),long_name=trim(long_name),missing_value=MV,grid_mapping="")

        deallocate(dat2D)

        return

    end subroutine chion_write_static_col

    subroutine chion_scatter_to_grid(dat2D,dat,grd)
        ! Chion.jl scatter_to_grid (src/io.jl:45).
        !
        !     out = fill(missing, grid_shape)
        !     out[js[i], is[i]] = values[i]
        !
        ! chion's array is (nx,ny) rather than Julia's (ny,nx) so that it can
        ! be written dim1="xc",dim2="yc"; the mapping (column i -> cell
        ! (js(i),is(i))) is identical. Cells no column maps to keep MV.
        !
        ! ALL columns are scattered, active or not, matching Chion.jl. An
        ! inactive column holds its last state and reporting it is correct;
        ! blanking it would be indistinguishable from an unmapped cell.

        implicit none

        real(wp), allocatable,  intent(OUT) :: dat2D(:,:)
        real(wp),               intent(IN)  :: dat(:)
        type(chion_grid_class), intent(IN)  :: grd

        ! Local variables
        integer :: i, nx, ny

        nx = size(grd%x)
        ny = size(grd%y)

        allocate(dat2D(nx,ny))
        dat2D = MV

        do i = 1, grd%ncol

            if (grd%is(i) .lt. 1 .or. grd%is(i) .gt. nx .or. &
                grd%js(i) .lt. 1 .or. grd%js(i) .gt. ny) then
                write(io_unit_err,*) "chion_scatter_to_grid:: Error: column maps outside the grid."
                write(io_unit_err,*) "icol, is, js, nx, ny = ", i, grd%is(i), grd%js(i), nx, ny
                stop "Program stopped."
            end if

            dat2D(grd%is(i),grd%js(i)) = dat(i)

        end do

        return

    end subroutine chion_scatter_to_grid

    ! =====================================================================
    ! Restart
    ! =====================================================================

    subroutine chion_restart_write(chn,filename,time)
        ! Write the FULL prognostic state, so that chion_restart_read followed
        ! by chion_update reproduces the trajectory exactly.
        !
        ! Always on the bare `column` dimension -- see the module header.

        implicit none

        type(chion_class), intent(IN) :: chn
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time

        ! Local variables
        integer :: ncid, i, ncol
        integer, allocatable :: active_int(:)

        ncol = chn%grd%ncol

        ! --- File, dimensions and identity ----------------------------

        call nc_create(filename)

        call nc_write_dim(filename,"column",x=1,dx=1,nx=ncol,units="1")

        if (trim(chn%par%model) .eq. "bessi") then
            call nc_write_dim(filename,"layer",x=1,dx=1,nx=chn%bsi%now%Ntot,units="1")
        end if

        call nc_write_dim(filename,"time",x=time,dx=1.0_wp,nx=1,units="days",unlimited=.TRUE.)

        ! The three facts chion_restart_read refuses to guess.
        call nc_write_attr(filename,"chion_model",trim(chn%par%model))
        call nc_write_attr(filename,"chion_ncol", ncol)
        if (trim(chn%par%model) .eq. "bessi") then
            call nc_write_attr(filename,"chion_Ntot",chn%bsi%now%Ntot)
        else
            call nc_write_attr(filename,"chion_Ntot",0)
        end if

        call nc_open(filename,ncid,writable=.TRUE.)

        call nc_write(filename,"time",time,dim1="time",start=[1],count=[1],grid_mapping="",ncid=ncid)

        ! --- Model-agnostic bookkeeping -------------------------------
        !
        ! smb_cum_prev and dt_last are NOT optional. chion_get_smb differences
        ! the cumulative accumulator against smb_cum_prev; restoring the
        ! accumulators but not the baseline makes the first chion_get_smb after
        ! a restart report the entire run-to-date as one step's flux
        ! (docs/porting_notes.md D13).

        allocate(active_int(ncol))
        active_int = 0
        do i = 1, ncol
            if (chn%grd%active(i)) active_int(i) = 1
        end do

        call nc_write(filename,"active",active_int,dim1="column",dim2="time", &
                      start=[1,1],units="1",long_name="Column active mask",grid_mapping="",ncid=ncid)

        call nc_write(filename,"smb_cum_prev",chn%smb_cum_prev,dim1="column",dim2="time", &
                      start=[1,1],units="kg m-2", &
                      long_name="Cumulative ice-facing mass flux at the start of the last step",grid_mapping="",ncid=ncid)

        call nc_write(filename,"dt_last",chn%dt_last,dim1="time",start=[1],count=[1], &
                      units="days",long_name="Length of the last completed step",grid_mapping="",ncid=ncid)

        deallocate(active_int)

        ! --- Model state ----------------------------------------------

        select case(trim(chn%par%model))

            case("bessi")

                ! Layer count and the four (Ntot,ncol) layer arrays.
                call nc_write(filename,"n_lay",chn%bsi%now%n_lay,dim1="column",dim2="time", &
                              start=[1,1],units="1",long_name="Number of active snow layers",grid_mapping="",ncid=ncid)

                call nc_write(filename,"mass",chn%bsi%now%mass,dim1="layer",dim2="column",dim3="time", &
                              start=[1,1,1],units="kg m-2",long_name="Layer snow mass",grid_mapping="",ncid=ncid)
                call nc_write(filename,"mass_w",chn%bsi%now%mass_w,dim1="layer",dim2="column",dim3="time", &
                              start=[1,1,1],units="kg m-2",long_name="Layer liquid-water mass",grid_mapping="",ncid=ncid)
                call nc_write(filename,"density",chn%bsi%now%density,dim1="layer",dim2="column",dim3="time", &
                              start=[1,1,1],units="kg m-3",long_name="Layer density",grid_mapping="",ncid=ncid)
                call nc_write(filename,"temperature",chn%bsi%now%temperature,dim1="layer",dim2="column",dim3="time", &
                              start=[1,1,1],units="K",long_name="Layer temperature",grid_mapping="",ncid=ncid)

                ! The eight wp_acc accumulators, as NF90_DOUBLE.
                call chion_restart_write_acc(filename,"mass_base",           chn%bsi%now%mass_base,           "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"smb_ice",             chn%bsi%now%smb_ice,             "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"runoff",              chn%bsi%now%runoff,              "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"melt",                chn%bsi%now%melt,                "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"refreezing",          chn%bsi%now%refreezing,          "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"vapor_mass",          chn%bsi%now%vapor_mass,          "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"sublimation",         chn%bsi%now%sublimation,         "kg m-2",  ncid)
                call chion_restart_write_acc(filename,"latent_heat_flux_sum",chn%bsi%now%latent_heat_flux_sum,"W m-2 d", ncid)

                ! Instantaneous per-column scalars.
                call chion_restart_write_col(filename,"t_srf", chn%bsi%now%t_srf, "K",ncid)
                call chion_restart_write_col(filename,"albedo",chn%bsi%now%albedo,"1",ncid)

                ! Diagnostics. Not prognostic -- summarize_domain_state
                ! recomputes them -- but bessi_reset_columns deliberately does
                ! NOT reset them (upstream defect 23), so a deactivated column
                ! carries stale values that ARE part of what a restart has to
                ! reproduce to be exact.
                call chion_restart_write_col(filename,"thickness",   chn%bsi%now%thickness,   "m",     ncid)
                call chion_restart_write_col(filename,"wet_mass",    chn%bsi%now%wet_mass,    "kg m-2",ncid)
                call chion_restart_write_col(filename,"bulk_density",chn%bsi%now%bulk_density,"kg m-3",ncid)
                call chion_restart_write_col(filename,"liquid_water",chn%bsi%now%liquid_water,"kg m-2",ncid)

            case("pdd")

                call chion_restart_write_col(filename,"snowpack_swe",chn%pdd%now%snowpack_swe,"kg m-2",ncid)

                call chion_restart_write_acc(filename,"smb_ice",chn%pdd%now%smb_ice,"kg m-2",ncid)
                call chion_restart_write_acc(filename,"runoff", chn%pdd%now%runoff, "kg m-2",ncid)
                call chion_restart_write_acc(filename,"pdd_sum",chn%pdd%now%pdd_sum,"K d",   ncid)

            case("itm")

                call chion_restart_write_col(filename,"H_snow",  chn%itm%now%H_snow,  "mm w.e.",   ncid)
                call chion_restart_write_col(filename,"alb_s",   chn%itm%now%alb_s,   "1",         ncid)
                call chion_restart_write_col(filename,"smb",     chn%itm%now%smb,     "mm w.e. d-1",ncid)
                call chion_restart_write_col(filename,"smbi",    chn%itm%now%smbi,    "mm w.e. d-1",ncid)
                call chion_restart_write_col(filename,"melt",    chn%itm%now%melt,    "mm w.e. d-1",ncid)
                call chion_restart_write_col(filename,"runoff",  chn%itm%now%runoff,  "mm w.e. d-1",ncid)
                call chion_restart_write_col(filename,"refrz",   chn%itm%now%refrz,   "mm w.e. d-1",ncid)
                call chion_restart_write_col(filename,"tsrf",    chn%itm%now%tsrf,    "K",         ncid)
                call chion_restart_write_col(filename,"melt_net",chn%itm%now%melt_net,"mm w.e. d-1",ncid)

                call chion_restart_write_acc(filename,"smb_cum",   chn%itm%now%smb_cum,   "mm w.e.",ncid)
                call chion_restart_write_acc(filename,"smbi_cum",  chn%itm%now%smbi_cum,  "mm w.e.",ncid)
                call chion_restart_write_acc(filename,"melt_cum",  chn%itm%now%melt_cum,  "mm w.e.",ncid)
                call chion_restart_write_acc(filename,"runoff_cum",chn%itm%now%runoff_cum,"mm w.e.",ncid)
                call chion_restart_write_acc(filename,"refrz_cum", chn%itm%now%refrz_cum, "mm w.e.",ncid)

            case DEFAULT

                call nc_close(ncid)
                call chion_io_model_error("chion_restart_write",chn%par%model)

        end select

        call nc_close(ncid)

        write(*,*)
        write(*,*) "time = ", time, " : saved chion restart file: ", trim(filename)
        write(*,*)

        return

    end subroutine chion_restart_write

    subroutine chion_restart_read(chn,filename,time)
        ! Read a restart file written by chion_restart_write into an ALREADY
        ! INITIALIZED chion object (chion_init must have run, so that the
        ! model, Ntot and ncol are known and the arrays are allocated).
        !
        ! `time` is intent(OUT): the restart file carries the time it was
        ! written at, and a host that resumes from a restart should resume from
        ! that time rather than assert one. This deviates from yelmo's
        ! yelmo_restart_read, where time is IN; yelmo can afford that because
        ! its time is managed by a separate timestepping object.
        !
        ! Refuses to load a file whose model, Ntot or ncol differ from the
        ! configured ones. A mismatch is always a configuration error, and the
        ! failure modes if it were allowed through range from a netCDF shape
        ! error (visible) to silently reading a PDD reservoir into a BESSI
        ! surface layer (not visible).

        implicit none

        type(chion_class), intent(INOUT) :: chn
        character(len=*),  intent(IN)    :: filename
        real(wp),          intent(OUT)   :: time

        ! Local variables
        integer :: ncid, i, ncol, nt
        character(len=56) :: file_model
        integer :: file_ncol, file_Ntot, Ntot
        integer, allocatable :: active_int(:)
        logical, allocatable :: active(:)

        call chion_check_file(filename)

        ncol = chn%grd%ncol

        Ntot = 0
        if (trim(chn%par%model) .eq. "bessi") Ntot = chn%bsi%now%Ntot

        ! --- Identity check, before anything is read ------------------

        if (.not. nc_exists_attr(filename,"chion_model")) then
            write(io_unit_err,*) "chion_restart_read:: Error: not a chion restart file."
            write(io_unit_err,*) "The global attribute 'chion_model' is missing."
            write(io_unit_err,*) "filename = ", trim(filename)
            stop "Program stopped."
        end if

        call nc_read_attr(filename,"chion_model",file_model)
        call nc_read_attr(filename,"chion_ncol", file_ncol)
        call nc_read_attr(filename,"chion_Ntot", file_Ntot)

        if (trim(file_model) .ne. trim(chn%par%model)) then
            write(io_unit_err,*) "chion_restart_read:: Error: restart file model does not match."
            write(io_unit_err,*) "filename        = ", trim(filename)
            write(io_unit_err,*) "restart model   = ", trim(file_model)
            write(io_unit_err,*) "configured model= ", trim(chn%par%model)
            stop "Program stopped."
        end if

        if (file_ncol .ne. ncol) then
            write(io_unit_err,*) "chion_restart_read:: Error: restart file ncol does not match."
            write(io_unit_err,*) "filename        = ", trim(filename)
            write(io_unit_err,*) "restart ncol    = ", file_ncol
            write(io_unit_err,*) "configured ncol = ", ncol
            stop "Program stopped."
        end if

        if (file_Ntot .ne. Ntot) then
            write(io_unit_err,*) "chion_restart_read:: Error: restart file Ntot does not match."
            write(io_unit_err,*) "filename        = ", trim(filename)
            write(io_unit_err,*) "restart Ntot    = ", file_Ntot
            write(io_unit_err,*) "configured Ntot = ", Ntot
            stop "Program stopped."
        end if

        ! --- Read the last record -------------------------------------

        nt = nc_size(filename,"time")

        call nc_open(filename,ncid,writable=.FALSE.)

        call nc_read(filename,"time",time,start=[nt],count=[1],ncid=ncid)

        allocate(active_int(ncol))
        allocate(active(ncol))

        call nc_read(filename,"active",active_int,start=[1,nt],count=[ncol,1],ncid=ncid)
        active = (active_int .eq. 1)

        call nc_read(filename,"smb_cum_prev",chn%smb_cum_prev,start=[1,nt],count=[ncol,1],ncid=ncid)
        call nc_read(filename,"dt_last",chn%dt_last,start=[nt],count=[1],ncid=ncid)

        select case(trim(chn%par%model))

            case("bessi")

                call nc_read(filename,"n_lay",chn%bsi%now%n_lay,start=[1,nt],count=[ncol,1],ncid=ncid)

                call nc_read(filename,"mass",       chn%bsi%now%mass,       start=[1,1,nt],count=[Ntot,ncol,1],ncid=ncid)
                call nc_read(filename,"mass_w",     chn%bsi%now%mass_w,     start=[1,1,nt],count=[Ntot,ncol,1],ncid=ncid)
                call nc_read(filename,"density",    chn%bsi%now%density,    start=[1,1,nt],count=[Ntot,ncol,1],ncid=ncid)
                call nc_read(filename,"temperature",chn%bsi%now%temperature,start=[1,1,nt],count=[Ntot,ncol,1],ncid=ncid)

                call chion_restart_read_acc(filename,"mass_base",           chn%bsi%now%mass_base,           nt,ncid)
                call chion_restart_read_acc(filename,"smb_ice",             chn%bsi%now%smb_ice,             nt,ncid)
                call chion_restart_read_acc(filename,"runoff",              chn%bsi%now%runoff,              nt,ncid)
                call chion_restart_read_acc(filename,"melt",                chn%bsi%now%melt,                nt,ncid)
                call chion_restart_read_acc(filename,"refreezing",          chn%bsi%now%refreezing,          nt,ncid)
                call chion_restart_read_acc(filename,"vapor_mass",          chn%bsi%now%vapor_mass,          nt,ncid)
                call chion_restart_read_acc(filename,"sublimation",         chn%bsi%now%sublimation,         nt,ncid)
                call chion_restart_read_acc(filename,"latent_heat_flux_sum",chn%bsi%now%latent_heat_flux_sum,nt,ncid)

                call chion_restart_read_col(filename,"t_srf", chn%bsi%now%t_srf, nt,ncid)
                call chion_restart_read_col(filename,"albedo",chn%bsi%now%albedo,nt,ncid)

                call chion_restart_read_col(filename,"thickness",   chn%bsi%now%thickness,   nt,ncid)
                call chion_restart_read_col(filename,"wet_mass",    chn%bsi%now%wet_mass,    nt,ncid)
                call chion_restart_read_col(filename,"bulk_density",chn%bsi%now%bulk_density,nt,ncid)
                call chion_restart_read_col(filename,"liquid_water",chn%bsi%now%liquid_water,nt,ncid)

            case("pdd")

                call chion_restart_read_col(filename,"snowpack_swe",chn%pdd%now%snowpack_swe,nt,ncid)

                call chion_restart_read_acc(filename,"smb_ice",chn%pdd%now%smb_ice,nt,ncid)
                call chion_restart_read_acc(filename,"runoff", chn%pdd%now%runoff, nt,ncid)
                call chion_restart_read_acc(filename,"pdd_sum",chn%pdd%now%pdd_sum,nt,ncid)

            case("itm")

                call chion_restart_read_col(filename,"H_snow",  chn%itm%now%H_snow,  nt,ncid)
                call chion_restart_read_col(filename,"alb_s",   chn%itm%now%alb_s,   nt,ncid)
                call chion_restart_read_col(filename,"smb",     chn%itm%now%smb,     nt,ncid)
                call chion_restart_read_col(filename,"smbi",    chn%itm%now%smbi,    nt,ncid)
                call chion_restart_read_col(filename,"melt",    chn%itm%now%melt,    nt,ncid)
                call chion_restart_read_col(filename,"runoff",  chn%itm%now%runoff,  nt,ncid)
                call chion_restart_read_col(filename,"refrz",   chn%itm%now%refrz,   nt,ncid)
                call chion_restart_read_col(filename,"tsrf",    chn%itm%now%tsrf,    nt,ncid)
                call chion_restart_read_col(filename,"melt_net",chn%itm%now%melt_net,nt,ncid)

                call chion_restart_read_acc(filename,"smb_cum",   chn%itm%now%smb_cum,   nt,ncid)
                call chion_restart_read_acc(filename,"smbi_cum",  chn%itm%now%smbi_cum,  nt,ncid)
                call chion_restart_read_acc(filename,"melt_cum",  chn%itm%now%melt_cum,  nt,ncid)
                call chion_restart_read_acc(filename,"runoff_cum",chn%itm%now%runoff_cum,nt,ncid)
                call chion_restart_read_acc(filename,"refrz_cum", chn%itm%now%refrz_cum, nt,ncid)

            case DEFAULT

                call nc_close(ncid)
                call chion_io_model_error("chion_restart_read",chn%par%model)

        end select

        call nc_close(ncid)

        ! Restore the active mask LAST, and through chion_grid_set_active
        ! rather than chion_set_active_mask: the latter RESETS columns that it
        ! sees switching off, which would wipe the state just read.
        call chion_grid_set_active(chn%grd,active)

        deallocate(active_int)
        deallocate(active)

        write(*,*)
        write(*,*) "time = ", time, " : loaded chion restart file: ", trim(filename)
        write(*,*)

        return

    end subroutine chion_restart_read

    ! ---------------------------------------------------------------------
    ! Restart element writers/readers.
    !
    ! The wp_acc pair exists so that "write the accumulators as double" is one
    ! decision in one place rather than eight opportunities to demote one by
    ! accident. Passing a wp_acc array to nc_write selects nc_write_double_*,
    ! i.e. NF90_DOUBLE in the file, with no rounding.
    ! ---------------------------------------------------------------------

    subroutine chion_restart_write_acc(filename,varname,dat,units,ncid)

        implicit none

        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: varname
        real(wp_acc),     intent(IN) :: dat(:)
        character(len=*), intent(IN) :: units
        integer,          intent(IN) :: ncid

        call nc_write(filename,trim(varname),dat,dim1="column",dim2="time", &
                      start=[1,1],units=trim(units),grid_mapping="",ncid=ncid)

        return

    end subroutine chion_restart_write_acc

    subroutine chion_restart_write_col(filename,varname,dat,units,ncid)

        implicit none

        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: varname
        real(wp),         intent(IN) :: dat(:)
        character(len=*), intent(IN) :: units
        integer,          intent(IN) :: ncid

        call nc_write(filename,trim(varname),dat,dim1="column",dim2="time", &
                      start=[1,1],units=trim(units),grid_mapping="",ncid=ncid)

        return

    end subroutine chion_restart_write_col

    subroutine chion_restart_read_acc(filename,varname,dat,nt,ncid)

        implicit none

        character(len=*), intent(IN)    :: filename
        character(len=*), intent(IN)    :: varname
        real(wp_acc),     intent(INOUT) :: dat(:)
        integer,          intent(IN)    :: nt
        integer,          intent(IN)    :: ncid

        call nc_read(filename,trim(varname),dat,start=[1,nt],count=[size(dat),1],ncid=ncid)

        return

    end subroutine chion_restart_read_acc

    subroutine chion_restart_read_col(filename,varname,dat,nt,ncid)

        implicit none

        character(len=*), intent(IN)    :: filename
        character(len=*), intent(IN)    :: varname
        real(wp),         intent(INOUT) :: dat(:)
        integer,          intent(IN)    :: nt
        integer,          intent(IN)    :: ncid

        call nc_read(filename,trim(varname),dat,start=[1,nt],count=[size(dat),1],ncid=ncid)

        return

    end subroutine chion_restart_read_col

    ! =====================================================================
    ! Metadata tables
    ! =====================================================================

    subroutine chion_load_table(model)
        ! Load the metadata table for `model`, once. Reloads if the model
        ! changes, so that two chion instances with different models in one
        ! executable both work (at the cost of a reload per alternation, which
        ! is a file read per output step -- acceptable, and the situation is
        ! hypothetical).

        implicit none

        character(len=*), intent(IN) :: model

        if (trim(var_table_model) .eq. trim(model) .and. allocated(var_table)) return

        select case(trim(model))
            case("bessi")
                call chion_check_file(table_bessi)
                call load_var_io_table(var_table,table_bessi)
            case("pdd")
                call chion_check_file(table_pdd)
                call load_var_io_table(var_table,table_pdd)
            case("itm")
                call chion_check_file(table_itm)
                call load_var_io_table(var_table,table_itm)
            case DEFAULT
                call chion_io_model_error("chion_load_table",model)
        end select

        var_table_model = trim(model)

        return

    end subroutine chion_load_table

    ! =====================================================================
    ! Errors
    ! =====================================================================

    subroutine chion_io_model_error(routine,model)

        implicit none

        character(len=*), intent(IN) :: routine
        character(len=*), intent(IN) :: model

        write(io_unit_err,*) trim(routine)//":: Error: model not recognized."
        write(io_unit_err,*) "model should be one of: ['bessi','pdd','itm']"
        write(io_unit_err,*) "model = ", trim(model)
        stop "Program stopped."

        return

    end subroutine chion_io_model_error

    subroutine chion_io_var_error(routine,varname,table)

        implicit none

        character(len=*), intent(IN) :: routine
        character(len=*), intent(IN) :: varname
        character(len=*), intent(IN) :: table

        write(io_unit_err,*) trim(routine)//":: Error: variable is in the metadata table &
                             &but has no case in the writer."
        write(io_unit_err,*) "varname = ", trim(varname)
        write(io_unit_err,*) "table   = ", trim(table)
        write(io_unit_err,*) "Add a case to "//trim(routine)//" in src/chion_io.f90, &
                             &or remove the row from the table."
        stop "Program stopped."

        return

    end subroutine chion_io_var_error

end module chion_io
