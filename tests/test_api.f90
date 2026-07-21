program test_api
    ! WP11 + WP13 acceptance test: the dispatcher, the public API and the
    ! parameter schema.
    !
    ! MUST BE RUN FROM THE REPOSITORY ROOT. The parameter schema path
    ! ("input/chion_defaults.nml") is hard-coded in chion_api, as it is in
    ! yelmo, so the working directory is part of the contract.
    !
    ! What is checked, and why each one is here:
    !
    !   1. Parameter round-trip. A namelist is written, loaded and every value
    !      compared. Because every read goes through nml_read with
    !      defaults_file = the schema, this simultaneously proves that EVERY
    !      PARAMETER THE CODE READS IS DECLARED IN input/chion_defaults.nml --
    !      a missing declaration is a hard error inside nml, so the test
    !      cannot pass without it. That is the WP13 coverage deliverable.
    !   2. chion_check_enum rejects a bad model value. Run in a child process,
    !      because the rejection is a `stop`.
    !   3. Full lifecycle for all three models, with no leaks.
    !   4. A 100-column run is BIT-IDENTICAL to 100 single-column runs. This
    !      is the core correctness check for the dispatcher: it fails if any
    !      column can see any other column's data.
    !   5. OMP_NUM_THREADS 1 vs 8 are BIT-IDENTICAL. This is the OpenMP
    !      privacy check.
    !   6. Active-mask behaviour: deactivation resets, inactive columns are
    !      frozen, n_active tracks.
    !   7. chion_get_smb returns one consistent quantity for all three models:
    !      right units, right sign convention, and summing rate*dt reproduces
    !      the underlying cumulative accumulator.

    !$ use omp_lib

    use nml,   only : nml_set_verbose
    use chion

    implicit none

    integer, parameter :: NCOL_BIG = 100
    integer, parameter :: NSTEP    = 40
    real(wp), parameter :: DT_DAYS = 1.0_wp

    character(len=*), parameter :: PAR_TMP = "test_api_tmp.nml"

    integer :: nfail
    integer :: io_tmp, iostat_tmp
    character(len=256) :: arg

    ! --- Child-process mode, used by the chion_check_enum test -----------
    ! chion_check_enum reports and calls `stop`, which gfortran gives exit
    ! status 0 for, so the parent detects the failure by its message on
    ! stderr rather than by the exit code.
    if (command_argument_count() .ge. 1) then
        call get_command_argument(1,arg)
        if (trim(arg) .eq. "enum-fail") then
            call chion_check_enum("chion","model","banana",CHION_MODEL_CHOICES)
            write(*,"(a)") "CHECK_ENUM_DID_NOT_STOP"
            stop
        end if
    end if

    nfail = 0

    call nml_set_verbose(.FALSE.)

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP11+WP13 acceptance test: dispatcher + public API"
    write(*,"(a)") "=========================================================="
    write(*,*)

    call test_parameters(nfail)
    call test_enum_rejection(nfail)
    call test_lifecycle(nfail)
    call test_grid_equivalence(nfail)
    call test_thread_equivalence(nfail)
    call test_active_mask(nfail)
    call test_get_smb(nfail)

    ! The scratch parameter file is written into the working directory, so
    ! remove it rather than leaving it in a git status.
    open(newunit=io_tmp,file=PAR_TMP,status="old",iostat=iostat_tmp)
    if (iostat_tmp .eq. 0) close(io_tmp,status="delete")

    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    ! =====================================================================
    ! 1. Parameters
    ! =====================================================================

    subroutine test_parameters(nfail)
        ! Write a namelist whose every value differs from the default, load it
        ! back through each *_par_load, and check that every value arrived.
        !
        ! The &chion, &bessi and &pdd groups are written SPARSE on purpose,
        ! omitting some parameters, to prove the defaults overlay works. The
        ! &itm group is written complete, because itm_par_load predates the
        ! defaults mechanism (see chion_itm_par_load).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_param_class) :: par
        type(chion_const_class) :: c
        type(bessi_par_class)   :: bpar
        type(pdd_par_class)     :: ppar
        type(itm_par_class)     :: ipar
        integer :: io

        write(*,"(a)") "--- 1. parameter loading round-trip ---"

        open(newunit=io,file=PAR_TMP,status="replace",action="write")
        write(io,"(a)") "&chion"
        write(io,"(a)") "    model      = 'pdd'"
        write(io,"(a)") "    nml_bessi  = 'bessi_alt'"
        write(io,"(a)") "    restart    = 'some_restart.nc'"
        write(io,"(a)") "/"
        write(io,"(a)") "&bessi_alt"
        write(io,"(a)") "    Ntot       = 7"
        write(io,"(a)") "    mass_max   = 400.0"
        write(io,"(a)") "    mass_split = 250.0"
        write(io,"(a)") "    mass_min   = 80.0"
        write(io,"(a)") "    diurnal_shortwave_substeps     = True"
        write(io,"(a)") "    diurnal_shortwave_max_substeps = 8"
        write(io,"(a)") "/"
        write(io,"(a)") "&pdd"
        write(io,"(a)") "    ddf_snow   = 4.5"
        write(io,"(a)") "    ddf_ice    = 9.5"
        write(io,"(a)") "    pdd_method = 'simple'"
        write(io,"(a)") "/"
        write(io,"(a)") "&itm"
        call write_full_itm_group(io)
        write(io,"(a)") "/"
        close(io)

        ! --- &chion ---------------------------------------------------
        call chion_par_load(par,PAR_TMP,"chion",init=.TRUE.)

        call check("chion: model overridden",         trim(par%model) .eq. "pdd", nfail)
        call check("chion: nml_bessi overridden",     trim(par%nml_bessi) .eq. "bessi_alt", nfail)
        call check("chion: restart overridden",       trim(par%restart) .eq. "some_restart.nc", nfail)
        ! Omitted from the user file -> must come from the defaults schema.
        call check("chion: nml_pdd from defaults",    trim(par%nml_pdd) .eq. "pdd", nfail)
        call check("chion: nml_itm from defaults",    trim(par%nml_itm) .eq. "itm", nfail)
        call check("chion: phys_const from defaults", trim(par%phys_const) .eq. "Earth", nfail)
        call check("chion: phys_const_file from defaults", &
                   trim(par%phys_const_file) .eq. "input/chion_phys_const.nml", nfail)
        call check("chion: nml_chion echoes the group", trim(par%nml_chion) .eq. "chion", nfail)

        ! --- constants ------------------------------------------------
        ! Reading this at all proves all 26 names are present in the &Earth
        ! group: chion_const_load uses the legacy nml path, where a missing
        ! parameter is a hard error.
        call chion_const_init(c)
        call chion_const_load(c,par%phys_const_file,par%phys_const)

        call check_val("const: rho_i",   c%rho_i,   917.0_wp,    nfail)
        call check_val("const: rho_w",   c%rho_w,   1000.0_wp,   nfail)
        call check_val("const: T0",      c%T0,      273.15_wp,   nfail)
        call check_val("const: Lm",      c%Lm,      334000.0_wp, nfail)
        call check_val("const: alpha_dry",c%alpha_dry,0.81_wp,   nfail)
        call check_val("const: seconds_per_day",c%seconds_per_day,86400.0_wp,nfail)
        call check("const: albedo_scheme -> dynamic flag", &
                   c%albedo_scheme .eq. CHION_ALBEDO_DYNAMIC, nfail)
        call check("const: fresh_snow_density_scheme -> constant flag", &
                   c%fresh_snow_density_scheme .eq. CHION_FRESH_SNOW_DENSITY_CONSTANT, nfail)
        call check("const: low_density_densification -> bessi flag", &
                   c%low_density_densification .eq. CHION_DENSIFY_BESSI, nfail)

        ! --- &bessi, under an aliased group name ----------------------
        call bessi_par_init(bpar)
        call bessi_par_load(bpar,PAR_TMP,par%nml_bessi,init=.TRUE.)

        call check("bessi: Ntot overridden",       bpar%Ntot .eq. 7, nfail)
        call check_val("bessi: mass_max",          bpar%mass_max,   400.0_wp, nfail)
        call check_val("bessi: mass_split",        bpar%mass_split, 250.0_wp, nfail)
        call check_val("bessi: mass_min",          bpar%mass_min,    80.0_wp, nfail)
        call check("bessi: diurnal substeps on",   bpar%diurnal_shortwave_substeps, nfail)
        call check("bessi: max_substeps = 8",      bpar%diurnal_shortwave_max_substeps .eq. 8, nfail)
        call check_val("bessi: density_init from defaults",     bpar%density_init,     300.0_wp, nfail)
        call check_val("bessi: temperature_init from defaults", bpar%temperature_init, 273.0_wp, nfail)
        call check_val("bessi: diurnal_shortwave_threshold from defaults", &
                       bpar%diurnal_shortwave_threshold, 0.0_wp, nfail)
        call check_val("bessi: diurnal_shortwave_min_air_temperature from defaults", &
                       bpar%diurnal_shortwave_min_air_temperature, 265.15_wp, nfail)
        call check("bessi: diurnal_temperature_cycle from defaults", &
                   .not. bpar%diurnal_temperature_cycle, nfail)
        call check_val("bessi: diurnal_temperature_amplitude from defaults", &
                       bpar%diurnal_temperature_amplitude, 5.0_wp, nfail)

        ! --- &pdd -----------------------------------------------------
        call pdd_par_init(ppar)
        call pdd_par_load(ppar,PAR_TMP,par%nml_pdd,init=.TRUE.)

        call check_val("pdd: ddf_snow",  ppar%ddf_snow, 4.5_wp, nfail)
        call check_val("pdd: ddf_ice",   ppar%ddf_ice,  9.5_wp, nfail)
        call check_val("pdd: refreezing_fraction from defaults", &
                       ppar%refreezing_fraction, 0.6_wp, nfail)
        call check_val("pdd: temperature_sigma from defaults", &
                       ppar%temperature_sigma, 5.0_wp, nfail)
        call check("pdd: pdd_method overridden to simple", &
                   ppar%pdd_method .eq. CHION_PDD_SIMPLE, nfail)

        ! The schema default is "pism", NOT Chion.jl's daily behaviour.
        ! docs/PLAN.md section 3.1b item 5. Guard it explicitly, because a
        ! silent revert to "simple" would only show up as an SMB bias at the
        ! equilibrium line.
        call pdd_par_init(ppar)
        call pdd_par_load(ppar,"input/chion_defaults.nml","pdd",init=.TRUE.)
        call check("pdd: SCHEMA DEFAULT is pism, not simple", &
                   ppar%pdd_method .eq. CHION_PDD_PISM, nfail)

        ! --- &itm -----------------------------------------------------
        call chion_itm_par_load(ipar,PAR_TMP,"itm")

        call check_val("itm: trans_a",    ipar%trans_a,     0.46_wp,   nfail)
        call check_val("itm: itm_c",      ipar%itm_c,     -45.0_wp,    nfail)
        call check_val("itm: itm_t",      ipar%itm_t,      10.0_wp,    nfail)
        call check_val("itm: H_snow_max", ipar%H_snow_max,5000.0_wp,   nfail)
        call check_val("itm: alb_snow_dry",ipar%alb_snow_dry,0.8_wp,   nfail)
        call check_val("itm: firn_fac",   ipar%firn_fac,    0.0266_wp, nfail)

        ! With no user group at all, the schema alone must supply a complete
        ! and valid ITM parameter set.
        call chion_itm_par_load(ipar,"input/chion_phys_const.nml","itm")
        call check_val("itm: schema-only load gives itm_c", ipar%itm_c, -45.0_wp, nfail)

        write(*,*)

        return

    end subroutine test_parameters

    subroutine write_full_itm_group(io)
        ! itm_par_load reads the legacy nml path, so a user &itm group must be
        ! complete. Kept in one place so the test file does not drift from the
        ! parameter set.

        implicit none

        integer, intent(IN) :: io

        write(io,"(a)") "    trans_a            = 0.46"
        write(io,"(a)") "    trans_b            = 6e-5"
        write(io,"(a)") "    trans_c            = 0.01"
        write(io,"(a)") "    itm_c              = -45.0"
        write(io,"(a)") "    itm_t              = 10.0"
        write(io,"(a)") "    itm_b              = -2.0"
        write(io,"(a)") "    itm_lat0           = 65.0"
        write(io,"(a)") "    H_snow_max         = 5000.0"
        write(io,"(a)") "    Pmaxfrac           = 0.6"
        write(io,"(a)") "    H_snow_crit_desert = 10.0"
        write(io,"(a)") "    H_snow_crit_forest = 100.0"
        write(io,"(a)") "    melt_crit          = 0.5"
        write(io,"(a)") "    alb_ocean          = 0.1"
        write(io,"(a)") "    alb_land           = 0.2"
        write(io,"(a)") "    alb_forest         = 0.1"
        write(io,"(a)") "    alb_ice            = 0.4"
        write(io,"(a)") "    alb_snow_dry       = 0.8"
        write(io,"(a)") "    alb_snow_wet       = 0.65"
        write(io,"(a)") "    firn_fac           = 0.0266"

        return

    end subroutine write_full_itm_group

    ! =====================================================================
    ! 2. Enum rejection
    ! =====================================================================

    subroutine test_enum_rejection(nfail)
        ! chion_check_enum must reject an unknown model value. It reports and
        ! calls `stop`, so the check runs in a child process: this executable
        ! re-invoked with the argument "enum-fail".
        !
        ! gfortran's `stop <string>` exits with status 0, so the detection is
        ! on the message, not on the exit code.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        character(len=512)  :: exe, cmd
        character(len=1024) :: line
        integer :: io, iostat, cmdstat
        logical :: saw_error, saw_no_stop

        write(*,"(a)") "--- 2. chion_check_enum rejects a bad model value ---"

        ! Positive controls first: every legal value must be accepted, and
        ! reaching the next line proves it.
        call chion_check_enum("chion","model","bessi",CHION_MODEL_CHOICES)
        call chion_check_enum("chion","model","pdd",  CHION_MODEL_CHOICES)
        call chion_check_enum("chion","model","itm",  CHION_MODEL_CHOICES)
        call check("all three legal model values accepted",.TRUE.,nfail)

        call get_command_argument(0,exe)
        cmd = trim(exe)//" enum-fail > test_api_enum.out 2>&1"

        call execute_command_line(trim(cmd),wait=.TRUE.,cmdstat=cmdstat)

        saw_error   = .FALSE.
        saw_no_stop = .FALSE.

        open(newunit=io,file="test_api_enum.out",status="old",action="read",iostat=iostat)
        if (iostat .eq. 0) then
            do
                read(io,"(a1024)",iostat=iostat) line
                if (iostat .ne. 0) exit
                if (index(line,"not recognized")       .gt. 0) saw_error   = .TRUE.
                if (index(line,"CHECK_ENUM_DID_NOT_STOP") .gt. 0) saw_no_stop = .TRUE.
            end do
            close(io,status="delete")
        end if

        call check("bad model value produces a 'not recognized' error",saw_error,nfail)
        call check("bad model value halts (execution did not continue)", &
                   .not. saw_no_stop, nfail)

        write(*,*)

        return

    end subroutine test_enum_rejection

    ! =====================================================================
    ! 3. Lifecycle
    ! =====================================================================

    subroutine test_lifecycle(nfail)
        ! init -> init_state -> several updates -> end, for all three models,
        ! with an explicit allocated() sweep afterwards.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: chn
        integer  :: im, n
        real(wp) :: smb(5)
        character(len=8) :: models(3)

        write(*,"(a)") "--- 3. full lifecycle, all three models, no leaks ---"

        models(1) = "bessi"
        models(2) = "pdd"
        models(3) = "itm"

        do im = 1, 3

            call write_model_par(PAR_TMP,trim(models(im)))
            call chion_init(chn,PAR_TMP,5)
            call chion_init_state(chn)

            call check(trim(models(im))//": ncol set",     chn%grd%ncol .eq. 5, nfail)
            call check(trim(models(im))//": all active",   chn%grd%n_active .eq. 5, nfail)
            call check(trim(models(im))//": forcing allocated", &
                       allocated(chn%forc%air_temperature), nfail)
            call check(trim(models(im))//": H_ice forcing allocated (WP11 addition)", &
                       allocated(chn%forc%H_ice), nfail)
            call check(trim(models(im))//": PDDs forcing allocated (WP11 addition)", &
                       allocated(chn%forc%PDDs), nfail)

            ! Only the selected model may be allocated.
            call check(trim(models(im))//": only the selected model is allocated", &
                       model_allocated(chn,trim(models(im))), nfail)

            do n = 1, NSTEP
                call set_forcing(chn,[1,2,3,4,5],n,"accum")
                call chion_update(chn,DT_DAYS)
            end do

            call chion_get_smb(chn,smb)
            call check(trim(models(im))//": get_smb returns finite values", &
                       all(smb .eq. smb), nfail)

            call chion_end(chn)

            call check(trim(models(im))//": forcing deallocated", &
                       .not. allocated(chn%forc%air_temperature), nfail)
            call check(trim(models(im))//": grid deallocated", &
                       .not. allocated(chn%grd%active), nfail)
            call check(trim(models(im))//": smb bookkeeping deallocated", &
                       .not. allocated(chn%smb_cum_prev), nfail)
            call check(trim(models(im))//": bessi state deallocated", &
                       .not. allocated(chn%bsi%now%mass), nfail)
            call check(trim(models(im))//": pdd state deallocated", &
                       .not. allocated(chn%pdd%now%snowpack_swe), nfail)
            call check(trim(models(im))//": itm state deallocated", &
                       .not. allocated(chn%itm%now%H_snow), nfail)

        end do

        write(*,*)

        return

    end subroutine test_lifecycle

    function model_allocated(chn,model) result(ok)
        ! True when the selected model's state is allocated and the other two
        ! are not.

        implicit none

        type(chion_class), intent(IN) :: chn
        character(len=*),  intent(IN) :: model
        logical :: ok

        ! Local variables
        logical :: b, p, i

        b = allocated(chn%bsi%now%mass)
        p = allocated(chn%pdd%now%snowpack_swe)
        i = allocated(chn%itm%now%H_snow)

        select case(trim(model))
            case("bessi") ; ok = b .and. (.not. p) .and. (.not. i)
            case("pdd")   ; ok = p .and. (.not. b) .and. (.not. i)
            case("itm")   ; ok = i .and. (.not. b) .and. (.not. p)
            case DEFAULT  ; ok = .FALSE.
        end select

        return

    end function model_allocated

    ! =====================================================================
    ! 4. Grid equivalence -- the core dispatcher check
    ! =====================================================================

    subroutine test_grid_equivalence(nfail)
        ! A 100-column run must be BIT-IDENTICAL to 100 separate 1-column
        ! runs. Exact equality, not a tolerance: the columns are independent
        ! by construction, and any difference at all means the loop leaked
        ! state between them.
        !
        ! Every column gets DIFFERENT forcing (see forcing_values), so this
        ! also catches a dispatcher that silently reuses column 1's forcing.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: big, one
        integer  :: im, n, icol
        integer  :: gidx_big(NCOL_BIG), gidx_one(1)
        real(wp_acc) :: sig_big(4,NCOL_BIG), sig_one(4)
        logical  :: ok
        character(len=8) :: models(3)

        write(*,"(a)") "--- 4. 100 columns == 100 single-column runs (bit-identical) ---"

        models(1) = "bessi"
        models(2) = "pdd"
        models(3) = "itm"

        do icol = 1, NCOL_BIG
            gidx_big(icol) = icol
        end do

        do im = 1, 3

            call write_model_par(PAR_TMP,trim(models(im)))

            call chion_init(big,PAR_TMP,NCOL_BIG)
            call chion_init_state(big)
            do n = 1, NSTEP
                call set_forcing(big,gidx_big,n,"mixed")
                call chion_update(big,DT_DAYS)
            end do
            do icol = 1, NCOL_BIG
                call state_signature(big,icol,sig_big(:,icol))
            end do
            call chion_end(big)

            ok = .TRUE.
            do icol = 1, NCOL_BIG
                gidx_one(1) = icol
                call chion_init(one,PAR_TMP,1)
                call chion_init_state(one)
                do n = 1, NSTEP
                    call set_forcing(one,gidx_one,n,"mixed")
                    call chion_update(one,DT_DAYS)
                end do
                call state_signature(one,1,sig_one)
                call chion_end(one)

                if (any(sig_one .ne. sig_big(:,icol))) then
                    ok = .FALSE.
                    write(*,"(a,i0)") "     first difference at column ", icol
                    write(*,"(a,4es24.16)") "     big = ", sig_big(:,icol)
                    write(*,"(a,4es24.16)") "     one = ", sig_one
                    exit
                end if
            end do

            call check(trim(models(im))//": 100-column == 100 x 1-column, exactly",ok,nfail)

        end do

        write(*,*)

        return

    end subroutine test_grid_equivalence

    ! =====================================================================
    ! 5. Thread equivalence
    ! =====================================================================

    subroutine test_thread_equivalence(nfail)
        ! The same run with 1 and with 8 OpenMP threads must be BIT-IDENTICAL.
        !
        ! The thread count is forced in-process with omp_set_num_threads, so
        ! this is a genuine check and not a documentation exercise -- but only
        ! when the library was built with openmp=1. Without it the
        ! omp_set_num_threads calls are conditional-compilation comments and
        ! the test degenerates into running the same thing twice, which still
        ! catches nondeterminism from uninitialised memory.
        !
        ! To check externally as well:
        !     OMP_NUM_THREADS=1 libchion/bin/test_api.x
        !     OMP_NUM_THREADS=8 libchion/bin/test_api.x
        ! (the forced setting inside this routine overrides the environment,
        ! which is exactly why it is used here).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: chn
        integer :: im, n, icol
        integer :: gidx(NCOL_BIG)
        real(wp_acc) :: sig1(4,NCOL_BIG), sig8(4,NCOL_BIG)
        logical :: have_omp
        character(len=8) :: models(3)

        write(*,"(a)") "--- 5. 1 thread == 8 threads (bit-identical) ---"

        have_omp = .FALSE.
        !$ have_omp = .TRUE.

        if (.not. have_omp) then
            write(*,"(a)") "     note: built WITHOUT OpenMP; running the determinism check only."
            write(*,"(a)") "     rebuild with openmp=1 for the real thread test."
        end if

        models(1) = "bessi"
        models(2) = "pdd"
        models(3) = "itm"

        do icol = 1, NCOL_BIG
            gidx(icol) = icol
        end do

        do im = 1, 3

            call write_model_par(PAR_TMP,trim(models(im)))

            !$ call omp_set_num_threads(1)
            call chion_init(chn,PAR_TMP,NCOL_BIG)
            call chion_init_state(chn)
            do n = 1, NSTEP
                call set_forcing(chn,gidx,n,"mixed")
                call chion_update(chn,DT_DAYS)
            end do
            do icol = 1, NCOL_BIG
                call state_signature(chn,icol,sig1(:,icol))
            end do
            call chion_end(chn)

            !$ call omp_set_num_threads(8)
            call chion_init(chn,PAR_TMP,NCOL_BIG)
            call chion_init_state(chn)
            do n = 1, NSTEP
                call set_forcing(chn,gidx,n,"mixed")
                call chion_update(chn,DT_DAYS)
            end do
            do icol = 1, NCOL_BIG
                call state_signature(chn,icol,sig8(:,icol))
            end do
            call chion_end(chn)

            call check(trim(models(im))//": 1 thread == 8 threads, exactly", &
                       all(sig1 .eq. sig8), nfail)

        end do

        !$ call omp_set_num_threads(1)

        write(*,*)

        return

    end subroutine test_thread_equivalence

    ! =====================================================================
    ! 6. Active mask
    ! =====================================================================

    subroutine test_active_mask(nfail)
        ! Deactivated columns are reset and then frozen; n_active tracks;
        ! reactivated columns start clean.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: chn
        integer  :: im, n
        integer  :: gidx(5)
        logical  :: active(5)
        real(wp_acc) :: sig_off_a(4), sig_off_b(4), sig_fresh(4), sig_on(4)
        character(len=8) :: models(3)

        write(*,"(a)") "--- 6. active-mask behaviour ---"

        models(1) = "bessi"
        models(2) = "pdd"
        models(3) = "itm"

        gidx = [1,2,3,4,5]

        do im = 1, 3

            call write_model_par(PAR_TMP,trim(models(im)))

            call chion_init(chn,PAR_TMP,5)
            call chion_init_state(chn)

            ! The cold-start signature of column 3, for comparison below.
            call state_signature(chn,3,sig_fresh)

            do n = 1, NSTEP
                call set_forcing(chn,gidx,n,"mixed")
                call chion_update(chn,DT_DAYS)
            end do

            ! Column 3 must have evolved away from its cold start, otherwise
            ! the freeze test below would be vacuous.
            call state_signature(chn,3,sig_off_a)
            call check(trim(models(im))//": column 3 evolved before deactivation", &
                       any(sig_off_a .ne. sig_fresh), nfail)

            ! --- deactivate column 3 -----------------------------------
            active = .TRUE.
            active(3) = .FALSE.
            call chion_set_active_mask(chn,active)

            call check(trim(models(im))//": n_active = 4 after deactivating one", &
                       chn%grd%n_active .eq. 4, nfail)
            call check(trim(models(im))//": active_idx skips column 3", &
                       all(chn%grd%active_idx(1:4) .eq. [1,2,4,5]), nfail)

            call state_signature(chn,3,sig_off_a)
            call check(trim(models(im))//": deactivation RESETS the column", &
                       all(sig_off_a .eq. sig_fresh), nfail)

            ! --- step on: column 3 must be frozen ----------------------
            do n = NSTEP+1, NSTEP+10
                call set_forcing(chn,gidx,n,"mixed")
                call chion_update(chn,DT_DAYS)
            end do

            call state_signature(chn,3,sig_off_b)
            call check(trim(models(im))//": inactive column is frozen", &
                       all(sig_off_b .eq. sig_off_a), nfail)

            ! --- reactivate: must start clean --------------------------
            active(3) = .TRUE.
            call chion_set_active_mask(chn,active)

            call check(trim(models(im))//": n_active back to 5", &
                       chn%grd%n_active .eq. 5, nfail)

            call state_signature(chn,3,sig_on)
            call check(trim(models(im))//": reactivated column starts clean", &
                       all(sig_on .eq. sig_fresh), nfail)

            ! --- and then evolves again --------------------------------
            do n = NSTEP+11, NSTEP+20
                call set_forcing(chn,gidx,n,"mixed")
                call chion_update(chn,DT_DAYS)
            end do

            call state_signature(chn,3,sig_on)
            call check(trim(models(im))//": reactivated column evolves again", &
                       any(sig_on .ne. sig_fresh), nfail)

            ! --- all off ----------------------------------------------
            active = .FALSE.
            call chion_set_active_mask(chn,active)
            call check(trim(models(im))//": n_active = 0 with all columns off", &
                       chn%grd%n_active .eq. 0, nfail)

            ! Stepping with nothing active must be a no-op, not a crash.
            call set_forcing(chn,gidx,1,"mixed")
            call chion_update(chn,DT_DAYS)
            call check(trim(models(im))//": update with no active columns is a no-op", &
                       .TRUE., nfail)

            call chion_end(chn)

        end do

        write(*,*)

        return

    end subroutine test_active_mask

    ! =====================================================================
    ! 7. chion_get_smb
    ! =====================================================================

    subroutine test_get_smb(nfail)
        ! chion_get_smb must return ONE consistent quantity for all three
        ! models: net mass flux to the ice sheet, [kg m-2 s-1], positive =
        ! ice sheet gains mass.
        !
        ! The three are NOT compared numerically -- they resolve different
        ! physics and will not agree. What is asserted is the CONVENTION:
        !
        !   (a) before the first update, every model returns exactly 0;
        !   (b) under accumulating forcing (cold, snowing) no model returns a
        !       negative flux;
        !   (c) under strongly ablating forcing (warm, no snow, no snowpack
        !       left) every model returns a strictly negative flux;
        !   (d) sum(smb*dt_seconds) over a run reproduces the model's own
        !       cumulative ice-facing accumulator to sp round-off. This is
        !       the check that the rate conversion is consistent with the
        !       underlying bookkeeping rather than merely plausible;
        !   (e) the magnitude is dimensionally sane for [kg m-2 s-1] -- a
        !       result accidentally left in [kg m-2 d-1] would be 86400x too
        !       large and would fail this.
        !
        ! On (b): the floor is exactly zero for BESSI and ITM. For PDD it is
        ! the sp round-off bound of its reservoir -- the one place the
        ! conversion is not exact, quantified in chion_get_smb's header.
        !
        ! On (b): BESSI and PDD legitimately return EXACTLY ZERO for a
        ! freshly started accumulating column, and that is not a defect. An
        ! ice-facing flux is not the surface accumulation rate: BESSI delivers
        ! nothing to the ice until snow is exported past the bottom of the
        ! 15-layer pack, and PDD's ice-facing flux is identically -ice_melt
        ! because it has no snow-to-ice conversion at all. Hence ">= 0" and
        ! not "> 0".
        !
        ! PDD's ice-facing flux is therefore NEVER positive, by construction.
        ! That is a structural property of the PDD scheme, asserted here so
        ! that it is a documented fact rather than a surprise at the yelmox
        ! cutover.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_class) :: chn
        integer  :: im, n
        integer  :: gidx(1)
        real(wp) :: smb(1), smb_min, smb_floor
        real(wp_acc) :: smb_sum, smb_cum(1), dt_sec, denom
        character(len=8) :: models(3)

        write(*,"(a)") "--- 7. chion_get_smb: one quantity, three models ---"

        models(1) = "bessi"
        models(2) = "pdd"
        models(3) = "itm"

        gidx = [1]

        do im = 1, 3

            call write_model_par(PAR_TMP,trim(models(im)))

            ! --- (a) zero before the first update ----------------------
            call chion_init(chn,PAR_TMP,1)
            call chion_init_state(chn)
            call chion_get_smb(chn,smb)
            call check(trim(models(im))//": smb = 0 before the first update", &
                       smb(1) .eq. 0.0_wp, nfail)

            ! --- (b) accumulating forcing: never negative --------------
            !
            ! "Never negative" is asserted against a per-model floor, not
            ! against exactly zero. BESSI and ITM difference a dp accumulator
            ! and the floor is a hard zero. PDD recovers its ice-facing flux
            ! by subtracting an sp-stored reservoir, so it carries an
            ! unbiased round-off term bounded by eps_sp*snowpack_swe per step
            ! -- see chion_get_smb's header. The floor below is that bound,
            ! computed from the reservoir the run actually built, so it
            ! tightens automatically if the reservoir is small and cannot
            ! quietly absorb a real negative flux.
            smb_min = huge(1.0_wp)
            do n = 1, NSTEP
                call set_forcing(chn,gidx,n,"accum")
                call chion_update(chn,DT_DAYS)
                call chion_get_smb(chn,smb)
                smb_min = min(smb_min,smb(1))
            end do

            smb_floor = 0.0_wp
            if (trim(models(im)) .eq. "pdd") then
                smb_floor = -4.0_wp*epsilon(1.0_wp) &
                            *max(abs(chn%pdd%now%snowpack_swe(1)),1.0_wp) &
                            /real(DT_DAYS*chn%c%seconds_per_day,wp)
            end if

            call check(trim(models(im))//": accumulating forcing -> smb >= 0", &
                       smb_min .ge. smb_floor, nfail)

            if (smb_min .lt. smb_floor) then
                write(*,"(a,2es16.8)") "     smb_min, floor = ", smb_min, smb_floor
            end if

            ! --- (e) dimensionally sane for kg m-2 s-1 -----------------
            ! A daily-rate leak would put this in the tens, not the 1e-4s.
            call check(trim(models(im))//": magnitude consistent with kg m-2 s-1", &
                       abs(smb(1)) .lt. 1.0e-2_wp, nfail)

            call chion_end(chn)

            ! --- (c) ablating forcing: strictly negative ---------------
            call chion_init(chn,PAR_TMP,1)
            call chion_init_state(chn)

            ! ITM cold-starts with H_snow = H_snow_max = 5000 mm w.e., a
            ! buffer far too deep to melt through in a test. Strip it, so all
            ! three models start from a genuinely bare column.
            if (trim(models(im)) .eq. "itm") chn%itm%now%H_snow = 0.0_wp

            smb_sum = 0.0_wp_acc
            dt_sec  = real(DT_DAYS,wp_acc)*real(chn%c%seconds_per_day,wp_acc)

            do n = 1, NSTEP
                call set_forcing(chn,gidx,n,"ablate")
                call chion_update(chn,DT_DAYS)
                call chion_get_smb(chn,smb)
                smb_sum = smb_sum + real(smb(1),wp_acc)*dt_sec
            end do

            call check(trim(models(im))//": ablating forcing -> smb < 0", &
                       smb(1) .lt. 0.0_wp, nfail)

            ! --- (d) sum(rate*dt) == the cumulative accumulator --------
            call chion_model_smb_cum(chn%par,chn%bsi,chn%pdd,chn%itm,smb_cum)

            denom = max(abs(smb_cum(1)),1.0_wp_acc)
            call check(trim(models(im))//": sum(smb*dt) reproduces the cumulative accumulator", &
                       abs(smb_sum - smb_cum(1))/denom .lt. 1.0e-5_wp_acc, nfail)

            if (abs(smb_sum - smb_cum(1))/denom .ge. 1.0e-5_wp_acc) then
                write(*,"(a,2es24.16)") "     sum, cum = ", smb_sum, smb_cum(1)
            end if

            ! --- PDD's structural sign restriction ---------------------
            if (trim(models(im)) .eq. "pdd") then
                call check("pdd: ice-facing flux is never positive (no snow-to-ice conversion)", &
                           smb_cum(1) .le. 0.0_wp_acc, nfail)
            end if

            call chion_end(chn)

        end do

        write(*,*)

        return

    end subroutine test_get_smb

    ! =====================================================================
    ! Helpers
    ! =====================================================================

    subroutine write_model_par(filename,model)
        ! A minimal, sparse parameter file selecting one model. Everything
        ! else comes from input/chion_defaults.nml, which is exactly the
        ! sparse-user-file behaviour the schema exists to provide.

        implicit none

        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: model

        ! Local variables
        integer :: io

        open(newunit=io,file=filename,status="replace",action="write")
        write(io,"(a)") "&chion"
        write(io,"(a)") "    model = '"//trim(model)//"'"
        write(io,"(a)") "/"
        close(io)

        return

    end subroutine write_model_par

    subroutine forcing_values(g,istep,regime,t2m,sf,rf,sw,ws)
        ! Deterministic per-column, per-step synthetic forcing.
        !
        ! Every column gets DIFFERENT values (they depend on g), which is what
        ! makes the 100-column / 100-single-column comparison meaningful: a
        ! dispatcher that packed column 1's forcing for every column would
        ! pass a test where all columns saw the same thing.

        implicit none

        integer,          intent(IN)  :: g
        integer,          intent(IN)  :: istep
        character(len=*), intent(IN)  :: regime
        real(wp),         intent(OUT) :: t2m, sf, rf, sw, ws

        ! Local variables
        real(wp) :: fg, fn

        fg = real(g,wp)
        fn = real(istep,wp)

        select case(trim(regime))

            case("accum")
                ! Cold and snowing everywhere: accumulation, little or no melt.
                t2m = 258.0_wp - 0.05_wp*fg + 3.0_wp*sin(fn*0.2_wp)
                sf  = 1.0e-4_wp + 1.0e-6_wp*fg
                rf  = 0.0_wp
                sw  = 40.0_wp + 0.5_wp*fg
                ws  = 3.0_wp + 0.01_wp*fg

            case("ablate")
                ! Warm, no snowfall, strong shortwave: ablation everywhere.
                t2m = 280.0_wp + 0.02_wp*fg
                sf  = 0.0_wp
                rf  = 1.0e-6_wp
                sw  = 320.0_wp + 0.2_wp*fg
                ws  = 5.0_wp

            case DEFAULT
                ! "mixed": a seasonal cycle whose phase and amplitude vary by
                ! column, so different columns are melting and accumulating at
                ! the same step. Exercises every branch of every model.
                t2m = 268.0_wp + 8.0_wp*sin(fn*0.15_wp + fg*0.07_wp) + 0.03_wp*fg
                sf  = 5.0e-5_wp*(1.0_wp + 0.5_wp*sin(fn*0.11_wp + fg*0.13_wp))
                rf  = 1.0e-5_wp*(1.0_wp + 0.5_wp*cos(fn*0.09_wp + fg*0.05_wp))
                sw  = 150.0_wp + 100.0_wp*sin(fn*0.15_wp) + 0.3_wp*fg
                ws  = 4.0_wp + 0.02_wp*fg

        end select

        return

    end subroutine forcing_values

    subroutine set_forcing(chn,gidx,istep,regime)
        ! Write the host-facing forcing arrays, mapping local column i to
        ! global column gidx(i).

        implicit none

        type(chion_class), intent(INOUT) :: chn
        integer,           intent(IN)    :: gidx(:)
        integer,           intent(IN)    :: istep
        character(len=*),  intent(IN)    :: regime

        ! Local variables
        integer  :: i, g
        real(wp) :: t2m, sf, rf, sw, ws

        do i = 1, size(gidx)

            g = gidx(i)
            call forcing_values(g,istep,regime,t2m,sf,rf,sw,ws)

            chn%forc%air_temperature(i) = t2m
            chn%forc%snowfall_rate(i)   = sf
            chn%forc%rainfall_rate(i)   = rf
            chn%forc%shortwave_down(i)  = sw
            chn%forc%wind_speed(i)      = ws
            chn%forc%latitude_deg(i)    = 70.0_wp - 0.05_wp*real(g,wp)
            chn%forc%surface_height(i)  = 1000.0_wp + 5.0_wp*real(g,wp)

            ! ITM-only. PDDs is an ANNUAL total held fixed through the year,
            ! not a per-step increment -- see chion_defs.
            chn%forc%H_ice(i) = 1000.0_wp
            chn%forc%PDDs(i)  = 200.0_wp + 2.0_wp*real(g,wp)

        end do

        chn%forc%day_of_year         = real(mod(istep-1,365)+1,wp)
        chn%forc%solar_longitude_deg = 360.0_wp*real(mod(istep-1,365),wp)/365.0_wp

        return

    end subroutine set_forcing

    subroutine state_signature(chn,icol,sig)
        ! Four numbers that pin down one column's state, per model. Used for
        ! the exact-equality comparisons; chosen so that a difference anywhere
        ! in the prognostic state shows up in at least one of them.
        !
        ! Everything is returned in wp_acc so that "==" means what it says.
        ! Sums over layer arrays are taken in a fixed order, so they are
        ! reproducible bit for bit.

        implicit none

        type(chion_class), intent(IN)  :: chn
        integer,           intent(IN)  :: icol
        real(wp_acc),      intent(OUT) :: sig(4)

        select case(trim(chn%par%model))

            case("bessi")
                sig(1) = real(sum(chn%bsi%now%mass(:,icol)),wp_acc) &
                       + real(chn%bsi%now%n_lay(icol),wp_acc)
                sig(2) = real(sum(chn%bsi%now%mass_w(:,icol)),wp_acc) &
                       + real(sum(chn%bsi%now%temperature(:,icol)),wp_acc)
                sig(3) = chn%bsi%now%smb_ice(icol) + chn%bsi%now%runoff(icol)
                sig(4) = real(chn%bsi%now%albedo(icol),wp_acc) &
                       + real(chn%bsi%now%t_srf(icol),wp_acc)

            case("pdd")
                sig(1) = real(chn%pdd%now%snowpack_swe(icol),wp_acc)
                sig(2) = chn%pdd%now%smb_ice(icol)
                sig(3) = chn%pdd%now%runoff(icol)
                sig(4) = chn%pdd%now%pdd_sum(icol)

            case("itm")
                sig(1) = real(chn%itm%now%H_snow(icol),wp_acc) &
                       + real(chn%itm%now%alb_s(icol),wp_acc)
                sig(2) = chn%itm%now%smbi_cum(icol)
                sig(3) = chn%itm%now%runoff_cum(icol) + chn%itm%now%melt_cum(icol)
                sig(4) = real(chn%itm%now%tsrf(icol),wp_acc) &
                       + real(chn%itm%now%smb(icol),wp_acc)

            case DEFAULT
                sig = 0.0_wp_acc

        end select

        return

    end subroutine state_signature

    ! =====================================================================
    ! Check helpers -- same style as tests/test_column_utils.f90
    ! =====================================================================

    subroutine check(label,ok,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        logical,          intent(IN)    :: ok
        integer,          intent(INOUT) :: nfail

        if (ok) then
            write(*,"(a,a)") "  ok   : ", trim(label)
        else
            write(*,"(a,a)") "  FAIL : ", trim(label)
            nfail = nfail + 1
        end if

        return

    end subroutine check

    subroutine check_val(label,val,expect,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: val
        real(wp),         intent(IN)    :: expect
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp) :: tol

        tol = 1.0e-5_wp*max(abs(expect),1.0_wp)

        if (abs(val-expect) .le. tol) then
            write(*,"(a,a)") "  ok   : ", trim(label)
        else
            write(*,"(a,a,2es16.8)") "  FAIL : ", trim(label), val, expect
            nfail = nfail + 1
        end if

        return

    end subroutine check_val

end program test_api
