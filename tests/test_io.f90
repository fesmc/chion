program test_io
    ! WP14 acceptance test: NetCDF output and restart.
    !
    ! Five things are checked, in the order of how much they can cost if wrong:
    !
    !   (a) RESTART ROUND-TRIP IS EXACT, for all three models. Run N steps,
    !       write a restart, read it into a FRESH object, then step both on for
    !       another M steps with identical forcing and assert the two states
    !       are BIT-IDENTICAL, field by field. The fields are enumerated
    !       explicitly rather than sampled: a restart that misses one field is
    !       a silently wrong answer, and the only way to catch a missing field
    !       is to name every field.
    !
    !   (b) A restart is REFUSED when the model or Ntot do not match. Both are
    !       fatal by design, so they are exercised by re-running this
    !       executable as a child process and checking what it printed --
    !       there is no way to catch a `stop` in-process.
    !
    !   (c) chion_get_smb IMMEDIATELY AFTER A RESTART equals what it would have
    !       been without restarting. This is exactly what smb_cum_prev and
    !       dt_last are for (docs/porting_notes.md D13), and it is the one
    !       restart failure that produces a plausible-looking wrong number
    !       rather than an obviously wrong one.
    !
    !   (d) The written NetCDF has the expected dimensions, variables, units
    !       and long_names for each model. The expected strings are HARD-CODED
    !       here, not read from input/chion-variables-*.md, so that the test
    !       fails if a table entry is changed by accident -- a test that reads
    !       its own expectations from the file under test checks nothing.
    !
    !   (e) The spatial scatter puts each column at its own (js,is) and leaves
    !       every unmapped cell at the missing value.
    !
    ! Usage:
    !   test_io.x                  run the suite
    !   test_io.x mismatch-model   child mode for (b); expected to abort
    !   test_io.x mismatch-ntot    child mode for (b); expected to abort

    use chion
    use chion_io
    use ncio

    implicit none

    ! Scratch files, all in the CWD and all removed at the end of a successful
    ! run. Named with a common prefix so a failed run leaves them identifiable.
    character(len=*), parameter :: par_a      = "test_io_a.nml"
    character(len=*), parameter :: par_b      = "test_io_b.nml"
    character(len=*), parameter :: file_rst   = "test_io_restart.nc"
    character(len=*), parameter :: file_out   = "test_io_out.nc"
    character(len=*), parameter :: file_grid  = "test_io_grid.nc"
    character(len=*), parameter :: file_child = "test_io_child.log"

    integer, parameter :: NCOL_TEST = 3
    integer, parameter :: NSTEP_1   = 24     ! steps before the restart
    integer, parameter :: NSTEP_2   = 12     ! steps after it
    real(wp), parameter :: DT_TEST  = 1.0_wp

    character(len=512) :: mode
    integer :: narg, nfail

    narg = command_argument_count()
    mode = "all"
    if (narg .ge. 1) call get_command_argument(1,mode)

    select case(trim(mode))

        case("mismatch-model")
            call run_mismatch("model")

        case("mismatch-ntot")
            call run_mismatch("ntot")

        case DEFAULT

            nfail = 0

            write(*,"(a)") "=========================================================="
            write(*,"(a)") " chion WP14 acceptance test: chion_io"
            write(*,"(a)") "=========================================================="

            call test_restart_roundtrip("bessi",nfail)
            call test_restart_roundtrip("pdd",  nfail)
            call test_restart_roundtrip("itm",  nfail)

            call test_restart_refused(nfail)

            call test_output_structure("bessi",nfail)
            call test_output_structure("pdd",  nfail)
            call test_output_structure("itm",  nfail)

            call test_spatial_scatter(nfail)

            write(*,*)
            write(*,"(a)") "=========================================================="
            if (nfail .eq. 0) then
                write(*,"(a)") " chion WP14: ALL CHECKS PASSED"
                write(*,"(a)") "=========================================================="
                call cleanup()
            else
                write(*,"(a,i0,a)") " chion WP14: ", nfail, " CHECK(S) FAILED"
                write(*,"(a)") "=========================================================="
                stop 1
            end if

    end select

contains

    ! =====================================================================
    ! (a) Restart round-trip
    ! =====================================================================

    subroutine test_restart_roundtrip(model,nfail)

        implicit none

        character(len=*), intent(IN)    :: model
        integer,          intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: chn1, chn2
        real(wp) :: time_rst
        real(wp) :: smb1(NCOL_TEST), smb2(NCOL_TEST)
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- (a) restart round-trip, model = "//trim(model)//" ---"

        call write_par(par_a,model)

        ! --- Reference run, part 1 ------------------------------------
        call chion_init(chn1,par_a,NCOL_TEST)
        call chion_init_state(chn1)

        do k = 1, NSTEP_1
            call set_forcing(chn1,k)
            call chion_update(chn1,DT_TEST)
        end do

        ! Deactivate one column before writing, so that the active mask and a
        ! reset column are part of what the restart has to reproduce.
        call set_active_one_off(chn1)

        call chion_get_smb(chn1,smb1)

        call chion_restart_write(chn1,file_rst,real(NSTEP_1,wp)*DT_TEST)

        ! --- Fresh object, restart read -------------------------------
        call chion_init(chn2,par_a,NCOL_TEST)
        call chion_init_state(chn2)
        call chion_restart_read(chn2,file_rst,time_rst)

        call check_val("restart time recovered",time_rst,real(NSTEP_1,wp)*DT_TEST,nfail)
        call check("dt_last recovered exactly",chn2%dt_last .eq. chn1%dt_last,nfail)
        call check_eq_acc("smb_cum_prev",chn2%smb_cum_prev,chn1%smb_cum_prev,nfail)
        call check("active mask recovered", &
                   all(chn2%grd%active .eqv. chn1%grd%active),nfail)
        call check("n_active recovered",chn2%grd%n_active .eq. chn1%grd%n_active,nfail)

        call compare_state("after restart read",model,chn1,chn2,nfail)

        ! --- (c) chion_get_smb immediately after the restart ----------
        call chion_get_smb(chn2,smb2)
        call check("chion_get_smb after restart == before restart", &
                   all(smb2 .eq. smb1),nfail)

        ! --- Continue both, identical forcing -------------------------
        do k = NSTEP_1+1, NSTEP_1+NSTEP_2
            call set_forcing(chn1,k)
            call set_forcing(chn2,k)
            call chion_update(chn1,DT_TEST)
            call chion_update(chn2,DT_TEST)
        end do

        call compare_state("after "//trim(itoa(NSTEP_2))//" further steps",model,chn1,chn2,nfail)

        call chion_get_smb(chn1,smb1)
        call chion_get_smb(chn2,smb2)
        call check("chion_get_smb after further steps",all(smb2 .eq. smb1),nfail)

        call chion_end(chn1)
        call chion_end(chn2)

        return

    end subroutine test_restart_roundtrip

    subroutine compare_state(label,model,a,b,nfail)
        ! Every prognostic field of the selected model, enumerated. Exact
        ! equality: the restart writes wp fields as NF90_FLOAT and wp_acc
        ! fields as NF90_DOUBLE, so there is no rounding anywhere in the round
        ! trip and anything short of bit-identical is a bug.

        implicit none

        character(len=*), intent(IN)    :: label
        character(len=*), intent(IN)    :: model
        type(chion_class), intent(IN)   :: a, b
        integer,          intent(INOUT) :: nfail

        write(*,"(a)") "  ["//trim(label)//"]"

        select case(trim(model))

            case("bessi")

                call check_eq_i1 ("n_lay",               b%bsi%now%n_lay,               a%bsi%now%n_lay,               nfail)
                call check_eq_r2 ("mass",                b%bsi%now%mass,                a%bsi%now%mass,                nfail)
                call check_eq_r2 ("mass_w",              b%bsi%now%mass_w,              a%bsi%now%mass_w,              nfail)
                call check_eq_r2 ("density",             b%bsi%now%density,             a%bsi%now%density,             nfail)
                call check_eq_r2 ("temperature",         b%bsi%now%temperature,         a%bsi%now%temperature,         nfail)
                call check_eq_acc("mass_base",           b%bsi%now%mass_base,           a%bsi%now%mass_base,           nfail)
                call check_eq_acc("smb_ice",             b%bsi%now%smb_ice,             a%bsi%now%smb_ice,             nfail)
                call check_eq_acc("runoff",              b%bsi%now%runoff,              a%bsi%now%runoff,              nfail)
                call check_eq_acc("melt",                b%bsi%now%melt,                a%bsi%now%melt,                nfail)
                call check_eq_acc("refreezing",          b%bsi%now%refreezing,          a%bsi%now%refreezing,          nfail)
                call check_eq_acc("vapor_mass",          b%bsi%now%vapor_mass,          a%bsi%now%vapor_mass,          nfail)
                call check_eq_acc("sublimation",         b%bsi%now%sublimation,         a%bsi%now%sublimation,         nfail)
                call check_eq_acc("latent_heat_flux_sum",b%bsi%now%latent_heat_flux_sum,a%bsi%now%latent_heat_flux_sum,nfail)
                call check_eq_r1 ("t_srf",               b%bsi%now%t_srf,               a%bsi%now%t_srf,               nfail)
                call check_eq_r1 ("albedo",              b%bsi%now%albedo,              a%bsi%now%albedo,              nfail)
                call check_eq_r1 ("thickness",           b%bsi%now%thickness,           a%bsi%now%thickness,           nfail)
                call check_eq_r1 ("wet_mass",            b%bsi%now%wet_mass,            a%bsi%now%wet_mass,            nfail)
                call check_eq_r1 ("bulk_density",        b%bsi%now%bulk_density,        a%bsi%now%bulk_density,        nfail)
                call check_eq_r1 ("liquid_water",        b%bsi%now%liquid_water,        a%bsi%now%liquid_water,        nfail)

            case("pdd")

                call check_eq_r1 ("snowpack_swe",b%pdd%now%snowpack_swe,a%pdd%now%snowpack_swe,nfail)
                call check_eq_acc("smb_ice",     b%pdd%now%smb_ice,     a%pdd%now%smb_ice,     nfail)
                call check_eq_acc("runoff",      b%pdd%now%runoff,      a%pdd%now%runoff,      nfail)
                call check_eq_acc("pdd_sum",     b%pdd%now%pdd_sum,     a%pdd%now%pdd_sum,     nfail)

            case("itm")

                call check_eq_r1 ("H_snow",    b%itm%now%H_snow,    a%itm%now%H_snow,    nfail)
                call check_eq_r1 ("alb_s",     b%itm%now%alb_s,     a%itm%now%alb_s,     nfail)
                call check_eq_r1 ("smb",       b%itm%now%smb,       a%itm%now%smb,       nfail)
                call check_eq_r1 ("smbi",      b%itm%now%smbi,      a%itm%now%smbi,      nfail)
                call check_eq_r1 ("melt",      b%itm%now%melt,      a%itm%now%melt,      nfail)
                call check_eq_r1 ("runoff",    b%itm%now%runoff,    a%itm%now%runoff,    nfail)
                call check_eq_r1 ("refrz",     b%itm%now%refrz,     a%itm%now%refrz,     nfail)
                call check_eq_r1 ("tsrf",      b%itm%now%tsrf,      a%itm%now%tsrf,      nfail)
                call check_eq_r1 ("melt_net",  b%itm%now%melt_net,  a%itm%now%melt_net,  nfail)
                call check_eq_acc("smb_cum",   b%itm%now%smb_cum,   a%itm%now%smb_cum,   nfail)
                call check_eq_acc("smbi_cum",  b%itm%now%smbi_cum,  a%itm%now%smbi_cum,  nfail)
                call check_eq_acc("melt_cum",  b%itm%now%melt_cum,  a%itm%now%melt_cum,  nfail)
                call check_eq_acc("runoff_cum",b%itm%now%runoff_cum,a%itm%now%runoff_cum,nfail)
                call check_eq_acc("refrz_cum", b%itm%now%refrz_cum, a%itm%now%refrz_cum, nfail)

        end select

        return

    end subroutine compare_state

    ! =====================================================================
    ! (b) Mismatched restarts are refused
    ! =====================================================================

    subroutine test_restart_refused(nfail)
        ! Both refusals end in `stop`, so they are run as child processes and
        ! judged on what they printed. gfortran's `STOP <string>` exits with
        ! status 0, so the exit code carries no information and the message is
        ! the only signal available.

        implicit none

        integer, intent(INOUT) :: nfail

        write(*,*)
        write(*,"(a)") "--- (b) mismatched restarts are refused ---"

        call check_child("mismatch-model","restart file model does not match",nfail)
        call check_child("mismatch-ntot", "restart file Ntot does not match", nfail)

        return

    end subroutine test_restart_refused

    subroutine run_mismatch(kind)
        ! Child mode. Writes a BESSI/Ntot=15 restart, then tries to load it
        ! into a deliberately mismatched configuration. chion_restart_read is
        ! expected to abort; reaching the end of this routine is the failure.

        implicit none

        character(len=*), intent(IN) :: kind

        ! Local variables
        type(chion_class) :: chn
        real(wp) :: time_rst
        integer  :: k

        call write_par(par_a,"bessi")

        call chion_init(chn,par_a,NCOL_TEST)
        call chion_init_state(chn)
        do k = 1, 3
            call set_forcing(chn,k)
            call chion_update(chn,DT_TEST)
        end do
        call chion_restart_write(chn,file_rst,3.0_wp)
        call chion_end(chn)

        select case(trim(kind))
            case("model")
                call write_par(par_b,"pdd")
            case("ntot")
                call write_par(par_b,"bessi",Ntot=10)
        end select

        call chion_init(chn,par_b,NCOL_TEST)
        call chion_init_state(chn)

        ! Expected to stop inside this call.
        call chion_restart_read(chn,file_rst,time_rst)

        write(*,"(a)") "NO-REFUSAL: chion_restart_read accepted a mismatched restart file."

        return

    end subroutine run_mismatch

    subroutine check_child(child_mode,expect,nfail)

        implicit none

        character(len=*), intent(IN)    :: child_mode
        character(len=*), intent(IN)    :: expect
        integer,          intent(INOUT) :: nfail

        ! Local variables
        character(len=512)  :: exe
        character(len=1024) :: cmd, line
        integer :: io, iostat, stat
        logical :: found

        call get_command_argument(0,exe)

        cmd = trim(exe)//" "//trim(child_mode)//" > "//file_child//" 2>&1"
        call execute_command_line(trim(cmd),wait=.TRUE.,exitstat=stat)

        found = .FALSE.

        open(newunit=io,file=file_child,status="old",action="read",iostat=iostat)
        if (iostat .eq. 0) then
            do
                read(io,"(a1024)",iostat=iostat) line
                if (iostat .ne. 0) exit
                if (index(line,trim(expect)) .gt. 0) then
                    found = .TRUE.
                    exit
                end if
            end do
            close(io)
        end if

        call check(trim(child_mode)//" -> '"//trim(expect)//"'",found,nfail)

        return

    end subroutine check_child

    ! =====================================================================
    ! (d) Output file structure
    ! =====================================================================

    subroutine test_output_structure(model,nfail)

        implicit none

        character(len=*), intent(IN)    :: model
        integer,          intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: chn
        integer :: k

        write(*,*)
        write(*,"(a)") "--- (d) output file structure, model = "//trim(model)//" ---"

        call write_par(par_a,model)

        call chion_init(chn,par_a,NCOL_TEST)
        call chion_init_state(chn)

        call chion_write_init(chn,file_out,0.0_wp,"days")
        call chion_write_step(chn,file_out,0.0_wp)

        do k = 1, 3
            call set_forcing(chn,k)
            call chion_update(chn,DT_TEST)
            call chion_write_step(chn,file_out,real(k,wp)*DT_TEST)
        end do

        ! --- Dimensions ----------------------------------------------
        call check_int("dim column",nc_size(file_out,"column"),NCOL_TEST,nfail)
        call check_int("dim time",  nc_size(file_out,"time"),  4,         nfail)

        if (trim(model) .eq. "bessi") then
            call check_int("dim layer",nc_size(file_out,"layer"),chn%bsi%now%Ntot,nfail)
        end if

        ! --- Variables, units, long_names -----------------------------
        ! Expected values written out literally, matching Chion.jl
        ! NETCDF_METADATA (src/io.jl:3-32) where a counterpart exists.

        select case(trim(model))

            case("bessi")

                call check_var("thickness",           "m",         "Snow thickness",                       nfail)
                call check_var("wet_mass",            "mmWE",      "Snow wet mass",                        nfail)
                call check_var("bulk_density",        "kg m-3",    "Bulk snow density",                    nfail)
                call check_var("liquid_water",        "kg m-2",    "Liquid water mass",                    nfail)
                call check_var("mass_base",           "mmWE",      "Firn mass exported to the ice model",  nfail)
                call check_var("smb_ice",             "mmWE",      "Net mass forcing to the ice sheet",    nfail)
                call check_var("runoff",              "mmWE",      "Cumulative runoff",                    nfail)
                call check_var("melt",                "mmWE",      "Cumulative melt",                      nfail)
                call check_var("refreezing",          "mmWE",      "Cumulative refreezing",                nfail)
                call check_var("sublimation",         "mmWE",      "Cumulative sublimation",               nfail)
                call check_var("vapor_mass",          "mmWE",      "Cumulative surface vapour mass flux",  nfail)
                call check_var("latent_heat_flux_sum","W m-2",     "Integrated turbulent latent heat flux",nfail)
                call check_var("Tsrf",                "K",         "Surface temperature",                  nfail)
                call check_var("albedo",              "1",         "Surface albedo",                       nfail)
                call check_var("N",                   "1",         "Number of active snow layers",         nfail)
                call check_var("smb",                 "kg m-2 s-1","Net mass flux to the ice sheet",       nfail)
                call check_var("mass",                "kg m-2",    "Layer snow mass",                      nfail)
                call check_var("mass_w",              "kg m-2",    "Layer liquid-water mass",              nfail)
                call check_var("density",             "kg m-3",    "Layer density",                        nfail)
                call check_var("temperature",         "K",         "Layer temperature",                    nfail)

                ! Shape of a layered variable: (layer,column,time).
                call check_shape3("mass",chn%bsi%now%Ntot,NCOL_TEST,4,nfail)

            case("pdd")

                call check_var("snowpack_swe","mmWE",      "Snowpack water equivalent",        nfail)
                call check_var("smb_ice",     "mmWE",      "Net mass forcing to the ice sheet",nfail)
                call check_var("runoff",      "mmWE",      "Cumulative runoff",                nfail)
                call check_var("pdd_sum",     "degC day",  "Cumulative positive degree days",  nfail)
                call check_var("smb",         "kg m-2 s-1","Net mass flux to the ice sheet",   nfail)

                call check("no layer dimension for pdd", &
                           .not. nc_exists_var(file_out,"mass"),nfail)

            case("itm")

                call check_var("H_snow",   "mmWE",      "Snowpack thickness",                          nfail)
                call check_var("albedo",   "1",         "Surface albedo",                              nfail)
                call check_var("Tsrf",     "K",         "Surface temperature",                         nfail)
                call check_var("smb_ice",  "mmWE",      "Net mass forcing to the ice sheet",           nfail)
                call check_var("runoff",   "mmWE",      "Cumulative runoff",                           nfail)
                call check_var("melt",     "mmWE",      "Cumulative melt",                             nfail)
                call check_var("refreezing","mmWE",     "Cumulative refreezing",                       nfail)
                call check_var("smb_total","mmWE",      "Cumulative whole-column surface mass balance",nfail)
                call check_var("smb",      "kg m-2 s-1","Net mass flux to the ice sheet",              nfail)

        end select

        call chion_end(chn)

        return

    end subroutine test_output_structure

    subroutine check_var(varname,units,long_name,nfail)

        implicit none

        character(len=*), intent(IN)    :: varname, units, long_name
        integer,          intent(INOUT) :: nfail

        ! Local variables
        character(len=256) :: got_units, got_long

        if (.not. nc_exists_var(file_out,trim(varname))) then
            write(*,"(a,a)") "  FAIL : variable missing from output file: ", trim(varname)
            nfail = nfail + 1
            return
        end if

        got_units = ""
        got_long  = ""
        call nc_read_attr(file_out,trim(varname),"units",    got_units)
        call nc_read_attr(file_out,trim(varname),"long_name",got_long)

        if (trim(got_units) .ne. trim(units)) then
            write(*,"(a,a,a,a,a,a)") "  FAIL : ", trim(varname), " units = '", &
                                     trim(got_units), "' expected '", trim(units)//"'"
            nfail = nfail + 1
        else if (trim(got_long) .ne. trim(long_name)) then
            write(*,"(a,a,a,a,a,a)") "  FAIL : ", trim(varname), " long_name = '", &
                                     trim(got_long), "' expected '", trim(long_name)//"'"
            nfail = nfail + 1
        else
            write(*,"(a,a,a,a,a)") "  ok   : ", trim(varname), "  [", trim(got_units), "]"
        end if

        return

    end subroutine check_var

    subroutine check_shape3(varname,n1,n2,n3,nfail)

        implicit none

        character(len=*), intent(IN)    :: varname
        integer,          intent(IN)    :: n1, n2, n3
        integer,          intent(INOUT) :: nfail

        ! Local variables
        character(len=32), allocatable :: names(:)
        integer,           allocatable :: dims(:)

        call nc_dims(file_out,trim(varname),names,dims)

        if (size(dims) .ne. 3) then
            write(*,"(a,a,a,i0)") "  FAIL : ", trim(varname), " rank = ", size(dims)
            nfail = nfail + 1
            return
        end if

        call check_int(trim(varname)//" dim1 (layer)", dims(1),n1,nfail)
        call check_int(trim(varname)//" dim2 (column)",dims(2),n2,nfail)
        call check_int(trim(varname)//" dim3 (time)",  dims(3),n3,nfail)

        return

    end subroutine check_shape3

    ! =====================================================================
    ! (e) Spatial scatter
    ! =====================================================================

    subroutine test_spatial_scatter(nfail)
        ! Three columns on a 4x3 grid, deliberately non-contiguous and out of
        ! order, so that a writer which quietly assumed column i -> cell i, or
        ! which swapped is and js, cannot pass.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: nx = 4, ny = 3

        type(chion_class) :: chn
        real(wp) :: xc(nx), yc(ny), mask(ny,nx)
        integer  :: js(NCOL_TEST), is(NCOL_TEST)
        real(wp) :: got(nx,ny)
        integer  :: i, ix, iy, k
        logical  :: mapped(nx,ny)
        logical  :: ok_mapped, ok_missing

        write(*,*)
        write(*,"(a)") "--- (e) spatial scatter ---"

        do ix = 1, nx
            xc(ix) = real(ix,wp)*10.0_wp
        end do
        do iy = 1, ny
            yc(iy) = real(iy,wp)*20.0_wp
        end do
        mask = 1.0_wp

        ! column 1 -> (j=1,i=2), column 2 -> (j=3,i=4), column 3 -> (j=2,i=1)
        js = [1,3,2]
        is = [2,4,1]

        call write_par(par_a,"bessi")

        call chion_init(chn,par_a,NCOL_TEST)
        call chion_init_state(chn)
        call chion_set_grid(chn,xc,yc,js,is,mask=mask)

        call check("has_spatial set by chion_set_grid",chn%grd%has_spatial,nfail)

        ! Distinct forcing per column, so the three t_srf values are distinct
        ! and a mis-scatter cannot alias.
        do k = 1, 3
            call set_forcing(chn,k)
            call chion_update(chn,DT_TEST)
        end do

        call chion_write_init(chn,file_grid,0.0_wp,"days")
        call chion_write_step(chn,file_grid,0.0_wp)

        call check_int("grid dim xc",nc_size(file_grid,"xc"),nx,nfail)
        call check_int("grid dim yc",nc_size(file_grid,"yc"),ny,nfail)

        call nc_read(file_grid,"Tsrf",got,start=[1,1,1],count=[nx,ny,1])

        mapped     = .FALSE.
        ok_mapped  = .TRUE.
        do i = 1, NCOL_TEST
            mapped(is(i),js(i)) = .TRUE.
            if (got(is(i),js(i)) .ne. chn%bsi%now%t_srf(i)) then
                write(*,"(a,i0,a,i0,a,i0,a,g14.6,a,g14.6)") &
                    "  FAIL : column ", i, " at (j=", js(i), ",i=", is(i), ") = ", &
                    got(is(i),js(i)), " expected ", chn%bsi%now%t_srf(i)
                ok_mapped = .FALSE.
            end if
        end do
        call check("each column lands at its own (j,i)",ok_mapped,nfail)

        ! And the three values really are distinct, so the check above has
        ! discriminating power.
        call check("the three column values are distinct", &
                   chn%bsi%now%t_srf(1) .ne. chn%bsi%now%t_srf(2) .and. &
                   chn%bsi%now%t_srf(2) .ne. chn%bsi%now%t_srf(3) .and. &
                   chn%bsi%now%t_srf(1) .ne. chn%bsi%now%t_srf(3),nfail)

        ok_missing = .TRUE.
        do iy = 1, ny
        do ix = 1, nx
            if (.not. mapped(ix,iy)) then
                if (got(ix,iy) .ne. MV) then
                    write(*,"(a,i0,a,i0,a,g14.6)") "  FAIL : unmapped cell (", ix, ",", iy, &
                                                   ") = ", got(ix,iy)
                    ok_missing = .FALSE.
                end if
            end if
        end do
        end do
        call check("unmapped cells hold the missing value",ok_missing,nfail)

        ! A layered variable is scattered layer by layer onto the same map.
        call check_layered_scatter(chn,is,js,nx,ny,nfail)

        call chion_end(chn)

        return

    end subroutine test_spatial_scatter

    subroutine check_layered_scatter(chn,is,js,nx,ny,nfail)

        implicit none

        type(chion_class), intent(IN)    :: chn
        integer,           intent(IN)    :: is(:), js(:)
        integer,           intent(IN)    :: nx, ny
        integer,           intent(INOUT) :: nfail

        ! Local variables
        real(wp), allocatable :: got3(:,:,:)
        integer :: i, k, Ntot
        logical :: ok

        Ntot = chn%bsi%now%Ntot

        allocate(got3(nx,ny,Ntot))
        call nc_read(file_grid,"mass",got3,start=[1,1,1,1],count=[nx,ny,Ntot,1])

        ok = .TRUE.
        do i = 1, size(is)
        do k = 1, Ntot
            if (got3(is(i),js(i),k) .ne. chn%bsi%now%mass(k,i)) ok = .FALSE.
        end do
        end do

        call check("layered variable scattered correctly, all layers",ok,nfail)

        deallocate(got3)

        return

    end subroutine check_layered_scatter

    ! =====================================================================
    ! Fixtures
    ! =====================================================================

    subroutine write_par(filename,model,Ntot)
        ! A minimal, SPARSE parameter file. Everything not named here comes
        ! from input/chion_defaults.nml, which is exactly the property WP13
        ! built (chion_api.f90 header), so the test does not have to restate
        ! the whole schema and cannot drift from it.

        implicit none

        character(len=*),  intent(IN) :: filename
        character(len=*),  intent(IN) :: model
        integer, optional, intent(IN) :: Ntot

        ! Local variables
        integer :: io

        open(newunit=io,file=filename,status="replace",action="write")

        write(io,"(a)") "&chion"
        write(io,"(a)") "    model            = """//trim(model)//""""
        write(io,"(a)") "    phys_const_file  = ""input/chion_phys_const.nml"""
        write(io,"(a)") "    phys_const       = ""Earth"""
        write(io,"(a)") "    restart          = ""None"""
        write(io,"(a)") "/"

        if (present(Ntot)) then
            write(io,"(a)")       "&bessi"
            write(io,"(a,i0)")    "    Ntot         = ", Ntot
            write(io,"(a)")       "/"
        end if

        close(io)

        return

    end subroutine write_par

    subroutine set_forcing(chn,k)
        ! Deterministic forcing, a function of the step index only, so that two
        ! runs stepped through the same k see byte-identical input. Column i is
        ! offset in temperature so that the three columns never coincide.
        !
        ! The cycle is short (40 steps) and warm enough at its peak to melt,
        ! so a 24+12-step run exercises accumulation, melt, percolation and
        ! refreezing rather than only accumulation.

        implicit none

        type(chion_class), intent(INOUT) :: chn
        integer,           intent(IN)    :: k

        ! Local variables
        integer  :: i
        real(wp) :: phase, t_air, sw

        phase = 2.0_wp*3.14159265358979_wp*real(k,wp)/40.0_wp

        t_air = 268.0_wp + 8.0_wp*sin(phase)
        sw    = max(120.0_wp + 120.0_wp*sin(phase),0.0_wp)

        do i = 1, chn%grd%ncol
            chn%forc%air_temperature(i) = t_air + 1.5_wp*real(i-1,wp)
            chn%forc%shortwave_down(i)  = sw
            chn%forc%wind_speed(i)      = 4.0_wp
            chn%forc%latitude_deg(i)    = 70.0_wp
            chn%forc%surface_height(i)  = 1500.0_wp
            chn%forc%H_ice(i)           = 1000.0_wp
            chn%forc%PDDs(i)            = 200.0_wp

            if (chn%forc%air_temperature(i) .lt. 273.15_wp) then
                chn%forc%snowfall_rate(i) = 3.0e-5_wp
                chn%forc%rainfall_rate(i) = 0.0_wp
            else
                chn%forc%snowfall_rate(i) = 0.0_wp
                chn%forc%rainfall_rate(i) = 3.0e-5_wp
            end if
        end do

        chn%forc%day_of_year         = modulo(real(k,wp),365.0_wp) + 1.0_wp
        chn%forc%solar_longitude_deg = 360.0_wp*(chn%forc%day_of_year - 1.0_wp)/365.0_wp

        return

    end subroutine set_forcing

    subroutine set_active_one_off(chn)
        ! Switch the LAST column off. chion_set_active_mask resets it on the
        ! way out (mirroring Chion.jl), so the restart then has to reproduce a
        ! mask, a reset column AND the re-baselined smb_cum_prev -- three
        ! things at once that a naive restart gets wrong independently.

        implicit none

        type(chion_class), intent(INOUT) :: chn

        ! Local variables
        logical :: active(NCOL_TEST)

        active = .TRUE.
        active(NCOL_TEST) = .FALSE.

        call chion_set_active_mask(chn,active)

        return

    end subroutine set_active_one_off

    subroutine cleanup()

        implicit none

        call unlink_if_present(par_a)
        call unlink_if_present(par_b)
        call unlink_if_present(file_rst)
        call unlink_if_present(file_out)
        call unlink_if_present(file_grid)
        call unlink_if_present(file_child)

        return

    end subroutine cleanup

    subroutine unlink_if_present(filename)

        implicit none

        character(len=*), intent(IN) :: filename

        ! Local variables
        integer :: io
        logical :: exists

        inquire(file=filename,exist=exists)
        if (.not. exists) return

        open(newunit=io,file=filename,status="old")
        close(io,status="delete")

        return

    end subroutine unlink_if_present

    ! =====================================================================
    ! Check helpers (tests/test_column_utils.f90 style)
    ! =====================================================================

    subroutine check(label,condition,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        logical,          intent(IN)    :: condition
        integer,          intent(INOUT) :: nfail

        if (condition) then
            write(*,"(a,a)") "  ok   : ", trim(label)
        else
            write(*,"(a,a)") "  FAIL : ", trim(label)
            nfail = nfail + 1
        end if

        return

    end subroutine check

    subroutine check_int(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        integer,          intent(IN)    :: value, expected
        integer,          intent(INOUT) :: nfail

        if (value .eq. expected) then
            write(*,"(a,a,a,i0)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,i0,a,i0)") "  FAIL : ", trim(label), " = ", value, &
                                       " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_int

    subroutine check_val(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value, expected
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp) :: tol

        tol = max(abs(expected)*8.0_wp*epsilon(1.0_wp),tiny(1.0_wp))

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_val

    subroutine check_eq_r1(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value(:), expected(:)
        integer,          intent(INOUT) :: nfail

        call check("    "//trim(label)//" bit-identical",all(value .eq. expected),nfail)

        return

    end subroutine check_eq_r1

    subroutine check_eq_r2(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value(:,:), expected(:,:)
        integer,          intent(INOUT) :: nfail

        call check("    "//trim(label)//" bit-identical",all(value .eq. expected),nfail)

        return

    end subroutine check_eq_r2

    subroutine check_eq_acc(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value(:), expected(:)
        integer,          intent(INOUT) :: nfail

        call check("    "//trim(label)//" bit-identical (dp)",all(value .eq. expected),nfail)

        return

    end subroutine check_eq_acc

    subroutine check_eq_i1(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        integer,          intent(IN)    :: value(:), expected(:)
        integer,          intent(INOUT) :: nfail

        call check("    "//trim(label)//" identical",all(value .eq. expected),nfail)

        return

    end subroutine check_eq_i1

    function itoa(n) result(s)

        implicit none

        integer, intent(IN) :: n
        character(len=12) :: s

        write(s,"(i0)") n
        s = adjustl(s)

        return

    end function itoa

end program test_io
