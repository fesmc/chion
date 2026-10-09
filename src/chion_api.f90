module chion_api
    ! WP11 + WP13 -- the public API and its parameter loading.
    !
    ! This is the only module a host program needs to understand. The usage
    ! contract, in full:
    !
    !     use chion
    !
    !     type(chion_class) :: chn
    !
    !     call chion_init(chn,filename,ncol,group="chion",cnst=cnst)  ! params + allocate
    !     call chion_init_state(chn)                        ! cold start
    !     do while (...)
    !         ! host writes chn%forc%<field>(:) directly
    !         call chion_update(chn,dt_days)
    !     end do
    !     call chion_end(chn)
    !
    ! Conventions taken from yelmo (src/yelmo_ice.f90):
    !   * object first, parameter file second, everything else optional and
    !     keyword-addressable;
    !   * _init loads parameters and allocates, and computes NO state;
    !   * _init_state computes the initial state and nothing else;
    !   * _end deallocates and nothing else;
    !   * the namelist group name is an argument, so the package can be
    !     instantiated more than once in one executable.
    !
    ! Deviation from yelmo, deliberate (docs/PLAN.md section 2.1):
    ! chion_update takes dt_days and advances exactly one step, instead of
    ! taking a target time and sub-cycling to reach it. This matches Chion.jl's
    ! step!(integrator, dt_days, force_dt=true) and leaves the time loop with
    ! the host.
    !
    ! -----------------------------------------------------------------------
    ! Parameter loading (WP13)
    ! -----------------------------------------------------------------------
    ! Every read goes through nml_read with defaults_file = the canonical
    ! schema input/chion_defaults.nml, following ytherm_par_load
    ! (yelmo/src/yelmo_thermodynamics.f90:616). This has two consequences that
    ! are the whole point of the exercise:
    !
    !   1. A parameter the code reads but the defaults file does not declare
    !      is a HARD ERROR at run time, in nml_read_internal. So "every
    !      parameter read by the code appears in the defaults file" is not
    !      something a script has to police -- it is checked on every single
    !      run. tests/test_api.f90 asserts it explicitly as well, by loading
    !      every group of every model.
    !   2. The user's parameter file may be SPARSE. Anything it does not
    !      mention comes from the defaults. nml_validate then catches the
    !      opposite error, a user parameter that does not exist in the schema
    !      (a typo), which would otherwise be silently ignored.
    !
    ! bessi_par_load and pdd_par_load live in THIS module rather than in
    ! snow_bessi.f90 / snow_pdd.f90, where they belong, because WP11/WP13 may
    ! not modify the physics modules. itm_par_load already exists in
    ! snow_itm.f90 and is used from there, not duplicated (chion_itm_par_load
    ! below passes it the schema as defaults).

    !$ use omp_lib

    use nml, only : nml_read, nml_validate
    use phys_constants, only : phys_const_class, phys_const_load, &
                               phys_const_require, phys_const_get, sec_day

    use chion_defs, only : wp, wp_acc, io_unit_err, MV, &
                           chion_const_class, chion_param_class, &
                           chion_grid_class, chion_forcing_class, &
                           chion_const_init, chion_const_print, chion_const_validate, &
                           chion_forcing_alloc, chion_forcing_dealloc, &
                           chion_grid_init, chion_grid_dealloc, &
                           chion_grid_set_active, &
                           chion_albedo_scheme_flag, &
                           chion_semix_snow_albedo_flag, &
                           chion_seb_scheme_flag, chion_longwave_scheme_flag, &
                           chion_turbulent_flux_scheme_flag, chion_climberx_qsat_flag, &
                           chion_fresh_snow_density_scheme_flag, &
                           chion_densify_scheme_flag, &
                           chion_check_enum, chion_check_file

    use chion_model, only : CHION_MODEL_CHOICES, &
                            chion_pack_step_forcing, &
                            chion_model_alloc, chion_model_dealloc, &
                            chion_model_init_state, chion_model_reset_columns, &
                            chion_model_step, chion_model_smb_cum, &
                            chion_model_flux_cum, chion_model_cum_active, &
                            chion_model_surface, chion_model_surface_active, &
                            chion_model_summary_line, &
                            CHION_NFLUX, CHION_IMELT, CHION_IRUNOFF, &
                            CHION_IREFRZ, CHION_ISUBL

    use snow_bessi, only : bessi_class, bessi_par_class, bessi_par_init, bessi_par_validate
    use snow_pdd,   only : pdd_class, pdd_par_class, pdd_par_init, pdd_par_validate, &
                           pdd_method_flag
    use snow_itm,   only : itm_class, itm_par_class, itm_par_load

    implicit none

    private

    ! The canonical schema. Hard-coded, as in yelmo, because it is part of the
    ! package rather than part of a run configuration: a run that supplies a
    ! different schema is not running chion.
    character(len=*), parameter :: def_file = "input/chion_defaults.nml"

    ! Defaults-file group names. A user file may name its groups anything
    ! (see chion_init's `group` argument and the nml_* fields of
    ! chion_param_class); these are what those groups are validated against.
    character(len=*), parameter :: def_chion = "chion"
    character(len=*), parameter :: def_bessi = "bessi"
    character(len=*), parameter :: def_pdd   = "pdd"
    character(len=*), parameter :: def_itm   = "itm"
    character(len=*), parameter :: def_const = "chion_const"

    ! Allowed values, in one place so the error messages and the validation
    ! can never disagree. Canonical names only, as in Chion.jl since 03bb445
    ! (docs/porting_notes.md D4).
    character(len=*), parameter :: CHION_ALBEDO_CHOICES  = "constant|dynamic|prescribed|semix|aging"
    character(len=*), parameter :: CHION_SEMIX_SNOW_ALB_CHOICES = "warren_wiscombe|dang"
    character(len=*), parameter :: CHION_SEB_CHOICES     = "bessi|semix"
    character(len=*), parameter :: CHION_LONGWAVE_CHOICES = "graybody|cloud_proxy"
    character(len=*), parameter :: CHION_TURB_CHOICES    = "bessi|semix|climberx"
    character(len=*), parameter :: CHION_CLIMBERX_QSAT_CHOICES = "climberx|bessi"
    character(len=*), parameter :: CHION_RHOS_CHOICES    = "constant|parameterized"
    character(len=*), parameter :: CHION_DENSIFY_CHOICES = "bessi|htessel"
    character(len=*), parameter :: CHION_PDD_CHOICES     = "simple|pism"

    ! =====================================================================
    ! The object
    ! =====================================================================

    type chion_class
        ! docs/PLAN.md section 2.2. Defined HERE and not in chion_defs
        ! (docs/porting_notes.md D7) because it contains the three model state
        ! types, which belong to their own modules. Keeping chion_defs
        ! physics-free is what makes the Level 1 kernels independently
        ! testable.

        type(chion_param_class)   :: par    ! model choice, group names, paths
        type(chion_const_class)   :: c      ! physical constants
        type(chion_grid_class)    :: grd    ! ncol, optional (x,y,js,is,mask), active mask
        type(chion_forcing_class) :: forc   ! (ncol) arrays -- the host writes these

        type(bessi_class)         :: bsi    ! allocated only when par%model == "bessi"
        type(pdd_class)           :: pdd    ! allocated only when par%model == "pdd"
        type(itm_class)           :: itm    ! allocated only when par%model == "itm"

        ! --- Bookkeeping for chion_get_smb -----------------------------
        !
        ! All three models report an ice-facing mass flux only as a running
        ! total (BESSI and PDD) or as a total plus a rate (ITM). To hand the
        ! host a rate for the step it just took, the total at the START of
        ! that step has to be remembered. smb_cum_prev is refreshed at the top
        ! of every chion_update, so it is always exactly one step behind.
        !
        ! wp_acc, because it is differenced against a wp_acc accumulator and
        ! rounding it to sp would reintroduce precisely the cancellation the
        ! accumulators exist to avoid (docs/porting_notes.md D1).

        real(wp_acc), allocatable :: smb_cum_prev(:)   ! (ncol) [kg m-2]

        ! Same one-step-behind snapshot for the surface flux diagnostics
        ! (melt, runoff, refreezing, sublimation), so chion_get_surface_fluxes
        ! can hand back a per-step rate. Shape (CHION_NFLUX,ncol); fields a
        ! model does not resolve are held at MV. See chion_model_flux_cum.
        real(wp_acc), allocatable :: flux_cum_prev(:,:)  ! (CHION_NFLUX,ncol) [kg m-2]

        real(wp)                  :: dt_last           ! [d] length of the last step
    end type chion_class

    ! =====================================================================
    ! Public interface
    ! =====================================================================

    public :: chion_class

    ! Lifecycle
    public :: chion_init
    public :: chion_init_state
    public :: chion_update
    public :: chion_end

    ! Supporting
    public :: chion_set_grid
    public :: chion_set_active_mask
    public :: chion_get_smb
    public :: chion_get_surface
    public :: chion_get_surface_fluxes
    public :: chion_get_surface_flux_totals
    public :: chion_print_summary

    ! Parameter loading, exposed so a host or a test can load a group without
    ! constructing a chion_class.
    public :: chion_par_load
    public :: chion_const_load
    public :: chion_const_from_phys
    public :: bessi_par_load
    public :: pdd_par_load
    public :: chion_itm_par_load

contains

    ! =====================================================================
    ! Lifecycle
    ! =====================================================================

    subroutine chion_init(chn,filename,ncol,group,cnst)
        ! Load parameters and allocate. Computes NO state -- call
        ! chion_init_state next.
        !
        ! filename : the run's parameter file. May be sparse; anything it does
        !            not declare comes from input/chion_defaults.nml.
        ! ncol     : number of independent snowpack columns.
        ! group    : namelist group name for the chion block. Default "chion".
        !            The model sub-group names are read FROM that block
        !            (nml_bessi / nml_pdd / nml_itm), so a second instance can
        !            use a completely disjoint set of groups.
        ! cnst     : the host's physical constants (fesm-utils). When given,
        !            chion takes its shared constants from it, so it agrees
        !            with every other component of the program. When absent,
        !            chion loads them from &chion:phys_const_file (standalone).

        implicit none

        type(chion_class),                intent(INOUT) :: chn
        character(len=*),                 intent(IN)    :: filename
        integer,                          intent(IN)    :: ncol
        character(len=*),       optional, intent(IN)    :: group
        type(phys_const_class), optional, intent(IN)    :: cnst

        ! Local variables
        character(len=56)      :: nml_group
        integer                :: n_threads
        type(phys_const_class) :: cnst_own

        nml_group = def_chion
        if (present(group)) nml_group = trim(group)

        ! Fail early and clearly if either file is missing, rather than
        ! letting nml report a cryptic open error on the first read.
        call chion_check_file(filename)
        call chion_check_file(def_file)

        ! --- OpenMP status --------------------------------------------
        ! Recorded so chion_print_summary can state it. The column loop is
        ! correct either way; this is provenance, not a switch.
        chn%par%use_omp = .FALSE.
        !$ chn%par%use_omp = .TRUE.

        n_threads = 1
        !$ n_threads = omp_get_max_threads()

        ! --- Parameters -----------------------------------------------

        call chion_par_load(chn%par,filename,nml_group)

        ! --- Physical constants ---------------------------------------
        ! chion's own from its &chion_const group (sparse, schema defaults);
        ! the shared ones from the host's record, or, standalone, from
        ! phys_const_file (docs/porting_notes.md D28).

        call chion_const_init(chn%c)
        call chion_const_load(chn%c,filename,chn%par%nml_const)

        if (present(cnst)) then
            call chion_const_from_phys(chn%c,cnst)
            chn%par%phys_const_src = "host: "//trim(cnst%source)//" ["//trim(cnst%label)//"]"
        else
            call chion_check_file(chn%par%phys_const_file)
            call phys_const_load(cnst_own,chn%par%phys_const_file,group=chn%par%phys_const)
            call chion_const_from_phys(chn%c,cnst_own)
            chn%par%phys_const_src = trim(chn%par%phys_const_file)//" ["//trim(chn%par%phys_const)//"]"
        end if

        ! --- Model parameters -----------------------------------------
        ! Only the selected model's group is read. The other two models are
        ! never allocated, so reading their parameters would only create an
        ! opportunity for an unrelated typo to abort the run.

        select case(trim(chn%par%model))
            case("bessi")
                call bessi_par_init(chn%bsi%par)
                call bessi_par_load(chn%bsi%par,filename,chn%par%nml_bessi)
            case("pdd")
                call pdd_par_init(chn%pdd%par)
                call pdd_par_load(chn%pdd%par,filename,chn%par%nml_pdd)
            case("itm")
                call chion_itm_par_load(chn%itm%par,filename,chn%par%nml_itm)
        end select

        ! --- Grid, forcing, model state -------------------------------

        call chion_grid_init(chn%grd,ncol)
        call chion_forcing_alloc(chn%forc,ncol)
        call chion_model_alloc(chn%par,chn%bsi,chn%pdd,chn%itm,ncol)

        if (allocated(chn%smb_cum_prev)) deallocate(chn%smb_cum_prev)
        allocate(chn%smb_cum_prev(ncol))
        chn%smb_cum_prev = 0.0_wp_acc

        if (allocated(chn%flux_cum_prev)) deallocate(chn%flux_cum_prev)
        allocate(chn%flux_cum_prev(CHION_NFLUX,ncol))
        chn%flux_cum_prev = 0.0_wp_acc

        chn%dt_last      = 0.0_wp

        call chion_print_summary(chn,n_threads)

        return

    end subroutine chion_init

    subroutine chion_init_state(chn)
        ! Cold start.
        !
        ! Restart reading is deliberately NOT folded in here. chion_io depends
        ! on this module for chion_class, so chion_api cannot call chion_io
        ! without a circular dependency. Restart is therefore host-driven, the
        ! same pattern yelmox uses for smbpal:
        !
        !     call chion_init(chn,filename,ncol)
        !     if (trim(chn%par%restart) .ne. "none") then
        !         call chion_restart_read(chn,chn%par%restart,time)
        !     else
        !         call chion_init_state(chn)
        !     end if
        !
        ! `use chion` exposes both. par%restart is read from the namelist and
        ! consumed by the host or driver, not by this routine.

        implicit none

        type(chion_class), intent(INOUT) :: chn

        call chion_model_init_state(chn%par,chn%c,chn%bsi,chn%pdd,chn%itm)

        ! Bring the SMB baseline in line with the state that was just built,
        ! so a chion_get_smb before the first chion_update reports zero rather
        ! than the whole initial accumulator.
        call chion_model_smb_cum(chn%par,chn%bsi,chn%pdd,chn%itm,chn%smb_cum_prev)
        call chion_model_flux_cum(chn%par,chn%bsi,chn%pdd,chn%itm,chn%flux_cum_prev)
        chn%dt_last = 0.0_wp

        return

    end subroutine chion_init_state

    subroutine chion_update(chn,dt_days)
        ! Advance every active column by exactly one step of dt_days, using
        ! whatever the host has written into chn%forc since the last call.

        implicit none

        type(chion_class), intent(INOUT) :: chn
        real(wp),          intent(IN)    :: dt_days

        ! Snapshot BEFORE stepping. This is the only place smb_cum_prev is
        ! written during a run, which is what makes chion_get_smb mean "the
        ! step that just finished" and nothing else. Active columns only, in
        ! parallel: inactive columns are not stepped, so their snapshot stays
        ! valid (and is re-baselined on reset by chion_set_active_mask).
        call chion_model_cum_active(chn%par,chn%grd,chn%bsi,chn%pdd,chn%itm, &
                                    chn%smb_cum_prev,chn%flux_cum_prev)

        call chion_model_step(chn%par,chn%c,chn%grd,chn%forc, &
                              chn%bsi,chn%pdd,chn%itm,dt_days)

        chn%dt_last = dt_days

        return

    end subroutine chion_update

    subroutine chion_end(chn)
        ! Deallocate everything and nothing else. After this call every
        ! allocatable component of chn is unallocated, which the acceptance
        ! test checks directly.

        implicit none

        type(chion_class), intent(INOUT) :: chn

        call chion_model_dealloc(chn%bsi,chn%pdd,chn%itm)
        call chion_forcing_dealloc(chn%forc)
        call chion_grid_dealloc(chn%grd)

        if (allocated(chn%smb_cum_prev))  deallocate(chn%smb_cum_prev)
        if (allocated(chn%flux_cum_prev)) deallocate(chn%flux_cum_prev)
        chn%dt_last = 0.0_wp

        return

    end subroutine chion_end

    ! =====================================================================
    ! Supporting routines
    ! =====================================================================

    subroutine chion_set_grid(chn,x,y,js,is,mask)
        ! Attach an optional 2-D spatial mapping to the column list, for IO
        ! only (Chion.jl SnowpackGrid, src/domain.jl:20-27). The physics never
        ! reads x, y, js, is or mask.
        !
        ! Supply either all of x, y, js, is, or none. The active mask is
        ! preserved across the call -- setting the output grid must not
        ! silently switch columns back on.

        implicit none

        type(chion_class), intent(INOUT) :: chn
        real(wp),          intent(IN)    :: x(:)
        real(wp),          intent(IN)    :: y(:)
        integer,           intent(IN)    :: js(:)
        integer,           intent(IN)    :: is(:)
        real(wp), optional,intent(IN)    :: mask(:,:)

        ! Local variables
        integer :: ncol
        logical, allocatable :: active_save(:)

        ncol = chn%grd%ncol

        if (ncol .le. 0) then
            write(io_unit_err,*) "chion_set_grid:: Error: call chion_init first."
            stop "Program stopped."
        end if

        allocate(active_save(ncol))
        active_save = chn%grd%active

        if (present(mask)) then
            call chion_grid_init(chn%grd,ncol,x=x,y=y,js=js,is=is,mask=mask)
        else
            call chion_grid_init(chn%grd,ncol,x=x,y=y,js=js,is=is)
        end if

        call chion_grid_set_active(chn%grd,active_save)

        deallocate(active_save)

        return

    end subroutine chion_set_grid

    subroutine chion_set_active_mask(chn,active)
        ! Switch columns on and off at run time. Mirrors Chion.jl
        ! set_active_mask! (src/integrators.jl) including the side effect that
        ! matters: columns that are being switched OFF are RESET to their
        ! cold-start state first.
        !
        ! Resetting on deactivation rather than on reactivation is Chion.jl's
        ! choice and is preserved. It has two consequences worth stating:
        !   * a deactivated column immediately stops reporting stale
        !     accumulators into the output; and
        !   * a column that is switched back on necessarily starts clean, so
        !     there is no way to "pause" a column and resume it.
        ! Columns that were already inactive are NOT reset again.

        implicit none

        type(chion_class), intent(INOUT) :: chn
        logical,           intent(IN)    :: active(:)

        ! Local variables
        integer :: i, n_off
        integer, allocatable :: idx_off(:)

        if (size(active) .ne. chn%grd%ncol) then
            write(io_unit_err,*) "chion_set_active_mask:: Error: mask length must equal ncol."
            write(io_unit_err,*) "ncol, size(active) = ", chn%grd%ncol, size(active)
            stop "Program stopped."
        end if

        ! Collect the columns making an active -> inactive transition.
        allocate(idx_off(chn%grd%ncol))
        n_off = 0
        do i = 1, chn%grd%ncol
            if (chn%grd%active(i) .and. (.not. active(i))) then
                n_off = n_off + 1
                idx_off(n_off) = i
            end if
        end do

        if (n_off .gt. 0) then
            call chion_model_reset_columns(chn%par,chn%c,chn%bsi,chn%pdd,chn%itm, &
                                           idx_off(1:n_off))

            ! Re-baseline the SMB snapshot for the reset columns. Without
            ! this, a chion_get_smb issued between chion_set_active_mask and
            ! the next chion_update would report the reset itself as a
            ! gigantic mass flux. Same for the surface-flux snapshot.
            call chion_model_smb_cum(chn%par,chn%bsi,chn%pdd,chn%itm,chn%smb_cum_prev)
            call chion_model_flux_cum(chn%par,chn%bsi,chn%pdd,chn%itm,chn%flux_cum_prev)
        end if

        deallocate(idx_off)

        call chion_grid_set_active(chn%grd,active)

        return

    end subroutine chion_set_active_mask

    subroutine chion_get_smb(chn,smb)
        ! ===================================================================
        ! THE model-agnostic surface mass balance accessor.
        !
        ! ===================================================================
        !   RETURNS: net mass flux to the ice sheet,   [kg m-2 s-1],
        !            POSITIVE = the ice sheet GAINS mass,
        !            averaged over the step chion_update just completed.
        ! ===================================================================
        !
        ! Per unit horizontal area of the column. A RATE, not a total: the
        ! consumer is an ice-sheet model, which must apply a flux over its own
        ! timestep, and that timestep is not chion's.
        !
        ! WHY [kg m-2 s-1]. It is the SI form, and it is already the unit
        ! chion_forcing_class uses for snowfall_rate and rainfall_rate, so
        ! every mass flux crossing the chion boundary -- in and out -- is in
        ! the same unit. No model's internal unit leaks into the API: BESSI
        ! and PDD accumulate [kg m-2] and ITM works in [mm w.e. d-1], and all
        ! three are converted here.
        !
        ! HOST CONVERSIONS. [kg m-2 s-1] == [mm w.e. s-1]. So:
        !     [mm w.e. d-1]  = smb * 86400
        !     [kg m-2 yr-1]  = smb * sec_year
        !     [m i.e. yr-1]  = smb * sec_year / rho_ice     <- yelmo bnd%smb
        ! with sec_year and rho_ice the HOST's, not chion's, so that a host
        ! whose year length or ice density differ from chion's constants gets
        ! its own convention exactly. That is the reason this routine does not
        ! return m i.e. yr-1 itself.
        !
        ! Inactive columns return 0. Before the first chion_update (dt_last
        ! still zero) every column returns 0.
        !
        ! ===================================================================
        ! Why this needed deciding
        ! ===================================================================
        ! The three models do not agree, in their native output, on what
        ! "SMB" means. WP9 established the disagreement and it is recorded in
        ! docs/porting_notes.md C2 and PLAN section 3.1b item 4:
        !
        !   BESSI  smb_ice   [kg m-2, cumulative]  -- ice-only forcing:
        !          bottom export, bare-ice deposition/sublimation, bare-surface
        !          melt, ice melt. Snow held in the pack is excluded.
        !   PDD    smb_ice   [kg m-2, cumulative]  -- whole-COLUMN mass change:
        !          it is additionally credited with the change in the snowpack
        !          reservoir, so a column that merely accumulates snow reports
        !          a positive "ice" mass balance.
        !   ITM    smb       [mm w.e. d-1, a RATE] -- whole-column again
        !          (sf + rf - runoff), while its smbi is the ice-facing one.
        !
        ! This accessor picks the ICE-SHEET-FACING definition, because that is
        ! the quantity a host ice-sheet model needs and the only one all three
        ! can supply. Nothing in any model's physics or state is changed; the
        ! conversion is arithmetic on outputs.
        !
        ! ===================================================================
        ! How each model maps onto it
        ! ===================================================================
        ! Let X be the cumulative ice-facing flux, extracted by
        ! chion_model_smb_cum (see that routine for the per-model expressions):
        !
        !   BESSI  X = smb_ice                     -- used unchanged; BESSI's
        !              accumulator is already exactly this quantity.
        !   PDD    X = smb_ice - snowpack_swe      -- the reservoir term is
        !              removed. Because
        !                  d(smb_ice)      = sf - snow_melt + refrozen - ice_melt
        !                  d(snowpack_swe) = sf - snow_melt + refrozen
        !              the difference is IDENTICALLY -ice_melt, which is the
        !              only mass PDD exchanges with the ice body: it has no
        !              densification and therefore no snow-to-ice conversion.
        !              This is an algebraic identity, not a calibration.
        !   ITM    X = smbi_cum                    -- smbpal's smbi,
        !              snow_to_ice + refrz - melted_ice, integrated over the
        !              run. NOT smbpal's smb, which is the whole-column form.
        !
        ! and the returned rate is
        !
        !     smb = (X_now - X_at_start_of_step) / (dt_days*sec_day)
        !
        ! X_at_start_of_step is chn%smb_cum_prev, refreshed at the top of
        ! every chion_update. The host tracks nothing.
        !
        ! For ITM this recovers itm%now%smbi/86400 exactly, because smbi_cum
        ! is accumulated as smbi*dt with the same dt.
        !
        ! All three mappings are exact identities on the models' own
        ! bookkeeping, so summing smb*dt_seconds over a run reproduces X to
        ! round-off. The acceptance test asserts exactly that, for all three.
        !
        ! ===================================================================
        ! The one thing that is approximated, and by how much
        ! ===================================================================
        ! BESSI and ITM difference a wp_acc (dp) accumulator, so their rates
        ! carry no conversion error at all.
        !
        ! PDD does not. Its reservoir snowpack_swe is stored in wp (sp), and
        ! the ice-facing flux is recovered by subtracting it. The subtraction
        ! itself is exact in dp, but the STORED reservoir has already been
        ! rounded to sp, so what comes back is
        !
        !     -ice_melt  +  (the sp rounding of d(snowpack_swe))
        !
        ! a noise term bounded by eps_sp*snowpack_swe per step, i.e. about
        !     1e-7 * 300 kg m-2 / 86400 s  ~  4e-10 kg m-2 s-1
        ! for a 300 kg m-2 pack. It errs in BOTH directions with no bias, so
        ! it does not accumulate; it cancels in the sum over a run, which is
        ! why check (d) in the acceptance test still passes. Its only visible
        ! consequence is that a PDD column which is purely accumulating can
        ! report a rate of order -1e-10 instead of exactly 0.
        !
        ! It is bounded by the reservoir size, so it grows if snowpack_swe
        ! does -- and snowpack_swe is unbounded upstream (defect 7). At
        ! 1e5 kg m-2 the noise reaches ~1e-7 kg m-2 s-1, still four orders
        ! below any real signal.
        !
        ! The clean fix is upstream and is not available from here: give
        ! pdd_state_class a cumulative ice_melt accumulator in wp_acc, which
        ! this routine would then read directly and exactly. Reported as a
        ! follow-up; it is an addition to snow_pdd.f90, which WP11/WP13 may
        ! not modify.
        !
        ! ===================================================================
        ! What the three will and will not agree on
        ! ===================================================================
        ! Given identical forcing the three will NOT return the same numbers,
        ! and are not expected to: BESSI resolves a 15-layer firn column, PDD
        ! resolves a bulk reservoir with no energy balance, ITM resolves an
        ! insolation-temperature melt model. What IS guaranteed, and what the
        ! acceptance test asserts, is that all three return the same PHYSICAL
        ! QUANTITY in the same units with the same sign convention, so a host
        ! can swap models without touching its coupling code.
        !
        ! One consequence worth stating plainly: over a cold accumulating
        ! column BESSI and PDD both return ~0, because neither delivers mass
        ! to the ice until snow reaches the bottom of the pack (BESSI) or ice
        ! melts (PDD). That is correct for an ice-facing flux and is NOT the
        ! same as the surface accumulation rate. A host that wants the
        ! whole-column balance must take it from the model state directly.

        implicit none

        type(chion_class), intent(IN)  :: chn
        real(wp),          intent(OUT) :: smb(:)

        ! Local variables
        integer :: i, icol
        real(wp_acc), allocatable :: smb_cum(:)
        real(wp_acc) :: dt

        if (size(smb) .ne. chn%grd%ncol) then
            write(io_unit_err,*) "chion_get_smb:: Error: smb must have length ncol."
            write(io_unit_err,*) "ncol, size(smb) = ", chn%grd%ncol, size(smb)
            stop "Program stopped."
        end if

        smb = 0.0_wp

        if (chn%dt_last .le. 0.0_wp) return

        allocate(smb_cum(chn%grd%ncol))
        call chion_model_smb_cum(chn%par,chn%bsi,chn%pdd,chn%itm,smb_cum)

        ! Days -> seconds using the SAME sec_day the models used, so
        ! the round trip smb*dt_seconds -> cumulative is exact.
        dt = real(chn%dt_last,wp_acc)*real(sec_day,wp_acc)

        ! Only active columns: an inactive column was not stepped, so its
        ! difference is zero anyway, but saying so explicitly means a column
        ! deactivated mid-run cannot leak a stale value.
        do i = 1, chn%grd%n_active
            icol = chn%grd%active_idx(i)
            smb(icol) = real((smb_cum(icol) - chn%smb_cum_prev(icol))/dt,wp)
        end do

        deallocate(smb_cum)

        return

    end subroutine chion_get_smb

    subroutine chion_get_surface(chn,t_srf,albedo)
        ! ===================================================================
        ! Model-agnostic surface-boundary accessor.
        !
        !   t_srf   surface skin temperature  [K]
        !   albedo  surface broadband albedo  [1]
        !
        ! The two fields an atmosphere needs to close its radiation over the
        ! surface: OLR / longwave from t_srf, planetary albedo from albedo.
        ! Both are instantaneous diagnostics of the surface state as of the
        ! last chion_update (not step-averaged, unlike chion_get_smb), so a
        ! host wanting a period mean averages calls after each step. Either
        ! output may be omitted; the work is parallel over active columns.
        !
        ! Inactive columns, and every column when the model exposes no such
        ! field (PDD), return MV. A caller must treat MV as "not supplied" and
        ! fall back to its own value -- the same contract PDD already imposes
        ! on the melt physics it does not resolve.
        ! ===================================================================

        implicit none

        type(chion_class),  intent(IN)  :: chn
        real(wp), optional, intent(OUT) :: t_srf(:)
        real(wp), optional, intent(OUT) :: albedo(:)

        integer :: icol

        call check_surface_arg("t_srf", t_srf, chn%grd%ncol)
        call check_surface_arg("albedo",albedo,chn%grd%ncol)

        ! MV on the inactive columns, in parallel: a serial whole-array fill
        ! dominated the cost of a per-step call on a large, mostly inactive grid.
        !$omp parallel do default(shared) private(icol)
        do icol = 1, chn%grd%ncol
            if (.not. chn%grd%active(icol)) then
                if (present(t_srf))  t_srf(icol)  = MV
                if (present(albedo)) albedo(icol) = MV
            end if
        end do
        !$omp end parallel do

        call chion_model_surface_active(chn%par,chn%grd,chn%bsi,chn%pdd,chn%itm, &
                                        t_srf=t_srf,albedo=albedo)

        return

    end subroutine chion_get_surface

    subroutine check_surface_arg(name,arr,ncol)
        ! Length guard for the optional chion_get_surface outputs.
        implicit none
        character(len=*),   intent(IN) :: name
        real(wp), optional, intent(IN) :: arr(:)
        integer,            intent(IN) :: ncol
        if (present(arr)) then
            if (size(arr) .ne. ncol) then
                write(io_unit_err,*) "chion_get_surface:: Error: "//trim(name)//" must have length ncol."
                write(io_unit_err,*) "ncol, size = ", ncol, size(arr)
                stop "Program stopped."
            end if
        end if
        return
    end subroutine check_surface_arg

    subroutine chion_get_surface_fluxes(chn,melt,runoff,refrz,subl)
        ! ===================================================================
        ! Surface flux diagnostics for the step chion_update just finished, as
        ! per-column RATES [kg m-2 s-1] (== [mm w.e. s-1]):
        !
        !   melt    total melt (snow + ice)
        !   runoff  runoff
        !   refrz   refreezing
        !   subl    sublimation (>= 0)
        !
        ! These are the surface (whole-column) budget terms, NOT the ice-facing
        ! flux chion_get_smb returns. A host builds the surface mass balance --
        ! the quantity MAR/RACMO report and the right target for validation --
        ! from the precipitation it forced chion with:
        !
        !     smb_surface = snowfall + rainfall - runoff - subl
        !
        ! Unlike the ice-facing flux, this needs no firn equilibration: it is an
        ! instantaneous surface budget, so a short atmosphere spin-up suffices.
        !
        ! Differenced against the one-step-behind snapshot flux_cum_prev, the
        ! same way and for the same reason as chion_get_smb. A field the model
        ! does not resolve (e.g. sublimation under ITM, everything but runoff
        ! under PDD) is MV, as is every inactive column and every column before
        ! the first chion_update. Optional arguments: request only what is
        ! needed.
        ! ===================================================================

        implicit none

        type(chion_class),  intent(IN)  :: chn
        real(wp), optional, intent(OUT) :: melt(:)
        real(wp), optional, intent(OUT) :: runoff(:)
        real(wp), optional, intent(OUT) :: refrz(:)
        real(wp), optional, intent(OUT) :: subl(:)

        integer :: i, icol, ncol
        real(wp_acc), allocatable :: flux_cum(:,:)
        real(wp_acc) :: dt

        ncol = chn%grd%ncol
        call check_flux_arg("melt",  melt,  ncol)
        call check_flux_arg("runoff",runoff,ncol)
        call check_flux_arg("refrz", refrz, ncol)
        call check_flux_arg("subl",  subl,  ncol)

        if (present(melt))   melt   = MV
        if (present(runoff)) runoff = MV
        if (present(refrz))  refrz  = MV
        if (present(subl))   subl   = MV

        if (chn%dt_last .le. 0.0_wp) return

        allocate(flux_cum(CHION_NFLUX,ncol))
        call chion_model_flux_cum(chn%par,chn%bsi,chn%pdd,chn%itm,flux_cum)

        dt = real(chn%dt_last,wp_acc)*real(sec_day,wp_acc)

        do i = 1, chn%grd%n_active
            icol = chn%grd%active_idx(i)
            if (present(melt))   call flux_rate(flux_cum,chn%flux_cum_prev,CHION_IMELT,  icol,dt,melt(icol))
            if (present(runoff)) call flux_rate(flux_cum,chn%flux_cum_prev,CHION_IRUNOFF,icol,dt,runoff(icol))
            if (present(refrz))  call flux_rate(flux_cum,chn%flux_cum_prev,CHION_IREFRZ, icol,dt,refrz(icol))
            if (present(subl))   call flux_rate(flux_cum,chn%flux_cum_prev,CHION_ISUBL,  icol,dt,subl(icol))
        end do

        deallocate(flux_cum)

        return

    end subroutine chion_get_surface_fluxes

    subroutine chion_get_surface_flux_totals(chn,melt,runoff,refrz,subl)
        ! ===================================================================
        ! Cumulative surface flux totals since the last cold start or column
        ! reset, per column, [kg m-2] (== [mm w.e.]):
        !
        !   melt    total melt (snow + ice)
        !   runoff  runoff
        !   refrz   refreezing
        !   subl    sublimation (>= 0)
        !
        ! The same quantities as chion_get_surface_fluxes, as running totals
        ! instead of the last step's rate. A host that aggregates over many
        ! steps (e.g. an annual cycle of daily steps) differences two calls
        ! instead of fetching a rate every step. Totals restart from zero when
        ! a column is reset (chion_set_active_mask), so set the mask before
        ! the first of the two calls. Unresolved fields and inactive columns
        ! are MV, as in chion_get_surface_fluxes.
        ! ===================================================================

        implicit none

        type(chion_class),  intent(IN)  :: chn
        real(wp_acc), optional, intent(OUT) :: melt(:)
        real(wp_acc), optional, intent(OUT) :: runoff(:)
        real(wp_acc), optional, intent(OUT) :: refrz(:)
        real(wp_acc), optional, intent(OUT) :: subl(:)

        integer :: i, icol, ncol
        real(wp_acc), allocatable :: flux_cum(:,:)

        ncol = chn%grd%ncol
        call check_total_arg("melt",  melt,  ncol)
        call check_total_arg("runoff",runoff,ncol)
        call check_total_arg("refrz", refrz, ncol)
        call check_total_arg("subl",  subl,  ncol)

        if (present(melt))   melt   = real(MV,wp_acc)
        if (present(runoff)) runoff = real(MV,wp_acc)
        if (present(refrz))  refrz  = real(MV,wp_acc)
        if (present(subl))   subl   = real(MV,wp_acc)

        allocate(flux_cum(CHION_NFLUX,ncol))
        call chion_model_flux_cum(chn%par,chn%bsi,chn%pdd,chn%itm,flux_cum)

        !$omp parallel do default(shared) private(i,icol)
        do i = 1, chn%grd%n_active
            icol = chn%grd%active_idx(i)
            if (present(melt))   melt(icol)   = flux_cum(CHION_IMELT,  icol)
            if (present(runoff)) runoff(icol) = flux_cum(CHION_IRUNOFF,icol)
            if (present(refrz))  refrz(icol)  = flux_cum(CHION_IREFRZ, icol)
            if (present(subl))   subl(icol)   = flux_cum(CHION_ISUBL,  icol)
        end do
        !$omp end parallel do

        deallocate(flux_cum)

        return

    end subroutine chion_get_surface_flux_totals

    subroutine check_total_arg(name,arr,ncol)
        ! Length guard for the optional chion_get_surface_flux_totals outputs.
        implicit none
        character(len=*),       intent(IN) :: name
        real(wp_acc), optional, intent(IN) :: arr(:)
        integer,                intent(IN) :: ncol
        if (present(arr)) then
            if (size(arr) .ne. ncol) then
                write(io_unit_err,*) "chion_get_surface_flux_totals:: Error: "//trim(name)//" must have length ncol."
                write(io_unit_err,*) "ncol, size = ", ncol, size(arr)
                stop "Program stopped."
            end if
        end if
        return
    end subroutine check_total_arg

    subroutine flux_rate(cum,prev,islot,icol,dt,out)
        ! Difference one packed flux slot into a rate, propagating the MV
        ! sentinel: a field the model does not resolve stays MV rather than
        ! becoming a spurious zero.
        implicit none
        real(wp_acc), intent(IN)  :: cum(:,:), prev(:,:)
        integer,      intent(IN)  :: islot, icol
        real(wp_acc), intent(IN)  :: dt
        real(wp),     intent(OUT) :: out
        if (cum(islot,icol) .eq. real(MV,wp_acc)) then
            out = MV
        else
            out = real((cum(islot,icol) - prev(islot,icol))/dt, wp)
        end if
        return
    end subroutine flux_rate

    subroutine check_flux_arg(name,arr,ncol)
        ! Length guard for the optional chion_get_surface_fluxes outputs.
        implicit none
        character(len=*),   intent(IN) :: name
        real(wp), optional, intent(IN) :: arr(:)
        integer,            intent(IN) :: ncol
        if (present(arr)) then
            if (size(arr) .ne. ncol) then
                write(io_unit_err,*) "chion_get_surface_fluxes:: Error: "//trim(name)//" must have length ncol."
                write(io_unit_err,*) "ncol, size = ", ncol, size(arr)
                stop "Program stopped."
            end if
        end if
        return
    end subroutine check_flux_arg

    subroutine chion_print_summary(chn,n_threads)
        ! One block of provenance for the run log: which model, how it is
        ! configured, how many columns, and whether OpenMP is live. Deliberately
        ! terse -- one line per topic, so a grep over a batch of logs is useful.

        implicit none

        type(chion_class), intent(IN) :: chn
        integer, optional, intent(IN) :: n_threads

        ! Local variables
        character(len=512) :: line
        integer            :: nt

        nt = 1
        if (present(n_threads)) nt = n_threads

        write(*,"(a)") ""
        write(*,"(a)") "== chion =================================================="
        write(*,"(a,a)")     " model      : ", trim(chn%par%model)
        write(*,"(a,i0)")    " ncol       : ", chn%grd%ncol
        write(*,"(a,i0)")    " n_active   : ", chn%grd%n_active
        write(*,"(a,l1,a,i0)") " openmp     : ", chn%par%use_omp, "   threads = ", nt
        write(*,"(a,a,a,a)") " constants  : ", trim(chn%par%phys_const_src), &
                             "  + chion's own group = ", trim(chn%par%nml_const)

        call chion_model_summary_line(chn%par,chn%bsi,chn%pdd,chn%itm,line)
        write(*,"(a,a)")     " ", trim(line)

        write(*,"(a)") "==========================================================="
        write(*,"(a)") ""

        return

    end subroutine chion_print_summary

    ! =====================================================================
    ! WP13 -- parameter loading
    ! =====================================================================

    subroutine chion_par_load(par,filename,group,init)
        ! The &chion block: which model, where the constants live, and what
        ! the model sub-groups are called.

        implicit none

        type(chion_param_class),    intent(INOUT) :: par
        character(len=*),           intent(IN)    :: filename
        character(len=*),           intent(IN)    :: group      ! usually "chion"
        logical,          optional, intent(IN)    :: init

        ! Local variables
        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        call nml_validate(filename,def_file,group,defaults_group=def_chion)

        ! Order matches the declaration order of chion_param_class and the
        ! order of the &chion group in input/chion_defaults.nml. Keep all
        ! three in step -- it is the only thing that makes a missing
        ! parameter obvious by eye.
        call nml_read(filename,group,"model",          par%model,          init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"nml_bessi",      par%nml_bessi,      init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"nml_pdd",        par%nml_pdd,        init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"nml_itm",        par%nml_itm,        init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"nml_const",      par%nml_const,      init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"phys_const_file",par%phys_const_file,init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"phys_const",     par%phys_const,     init=init_pars,defaults_file=def_file,defaults_group=def_chion)
        call nml_read(filename,group,"restart",        par%restart,        init=init_pars,defaults_file=def_file,defaults_group=def_chion)

        par%nml_chion = trim(group)

        call chion_check_enum(group,"model",par%model,CHION_MODEL_CHOICES)

        return

    end subroutine chion_par_load

    subroutine chion_const_load(c,filename,group,init)
        ! chion's own constants: the &chion_const group, read through the
        ! schema defaults (input/chion_defaults.nml) like every other chion
        ! group, so a run's par file lists only what it overrides.
        !
        ! The SHARED fields (rho_i, rho_w, ci, cw, Lm, grav, T0) are not read
        ! here: they come from the program's phys_const_class through
        ! chion_const_from_phys (docs/porting_notes.md D28).
        !
        ! The scheme flags are stored in this type rather than in
        ! chion_param_class, mirroring Chion.jl's SnowpackPhysicalConstants
        ! (docs/porting_notes.md D3), so they are read here too. They arrive
        ! as strings and are converted to integer flags (D4).

        implicit none

        type(chion_const_class),    intent(INOUT) :: c
        character(len=*),           intent(IN)    :: filename
        character(len=*),           intent(IN)    :: group        ! usually "chion_const"
        logical,          optional, intent(IN)    :: init

        ! Local variables
        logical           :: init_pars
        character(len=56) :: albedo_scheme
        character(len=56) :: semix_snow_albedo
        character(len=56) :: seb_scheme
        character(len=56) :: longwave_scheme
        character(len=56) :: turbulent_flux_scheme
        character(len=56) :: climberx_qsat
        character(len=56) :: fresh_snow_density_scheme
        character(len=56) :: low_density_densification

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        call nml_validate(filename,def_file,group,defaults_group=def_const)

        call nml_read(filename,group,"rho_s",                    c%rho_s,                   init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"rho_s_a",                  c%rho_s_a,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"rho_s_b",                  c%rho_s_b,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"rho_s_c",                  c%rho_s_c,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"fresh_snow_density_scheme",fresh_snow_density_scheme, init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"Lv",                       c%Lv,                      init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"cp_air",                   c%cp_air,                  init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"latent_heat_flux_ratio",   c%latent_heat_flux_ratio,  init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"D_sh",                     c%D_sh,                    init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"seb_scheme",               seb_scheme,                init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"turbulent_flux_scheme",    turbulent_flux_scheme,     init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"z0m_snow",                 c%z0m_snow,                init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"z0m_ice",                  c%z0m_ice,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"zm_to_zh",                 c%zm_to_zh,                init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"z_sfl",                    c%z_sfl,                   init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"karman",                   c%karman,                  init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"R_dry",                    c%R_dry,                   init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"l_neutral",                c%l_neutral,               init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"l_dew",                    c%l_dew,                   init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"climberx_qsat",            climberx_qsat,             init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"semix_karman",             c%semix_karman,            init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_surface_height",     c%semix_surface_height,    init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_z0m_snow",           c%semix_z0m_snow,          init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_z0m_ice",            c%semix_z0m_ice,           init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_zm_to_zh",           c%semix_zm_to_zh,          init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_sensible_exchange_factor", c%semix_sensible_exchange_factor, init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_stable_coefficient", c%semix_stable_coefficient, init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_latent_exchange_factor",   c%semix_latent_exchange_factor,   init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"alpha_dry",                c%alpha_dry,               init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"alpha_wet",                c%alpha_wet,               init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"alpha_ice",                c%alpha_ice,               init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"max_lwc_albedo",           c%max_lwc_albedo,          init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"albedo_scheme",            albedo_scheme,             init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"aging_cold_timescale_days",   c%aging_cold_timescale_days,   init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"aging_melting_timescale_days",c%aging_melting_timescale_days,init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"aging_snowfall_ref",       c%aging_snowfall_ref,      init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"frac_vu",                  c%frac_vu,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"alb_snow_vis_new",         c%alb_snow_vis_new,        init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"alb_snow_nir_new",         c%alb_snow_nir_new,        init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"snow_grain_fresh",         c%snow_grain_fresh,        init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"snow_grain_old",           c%snow_grain_old,          init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"d_alb_age_vis",            c%d_alb_age_vis,           init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"d_alb_age_nir",            c%d_alb_age_nir,           init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"f_age_t",                  c%f_age_t,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"dT_age",                   c%dT_age,                  init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"snow_0",                   c%snow_0,                  init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"snow_1",                   c%snow_1,                  init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"w_snow_dust",              c%w_snow_dust,             init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"dust_con_scale",           c%dust_con_scale,          init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"semix_snow_albedo",        semix_snow_albedo,         init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"dalb_snow_vis",            c%dalb_snow_vis,           init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"dalb_snow_nir",            c%dalb_snow_nir,           init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"k_sigma_orog",             c%k_sigma_orog,            init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"sigma_orog_crit",          c%sigma_orog_crit,         init=init_pars,defaults_file=def_file,defaults_group=def_const)

        call nml_read(filename,group,"eps_air",                  c%eps_air,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"longwave_scheme",          longwave_scheme,           init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"lw_emissivity_base",       c%lw_emissivity_base,      init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"lw_emissivity_temperature_slope",    c%lw_emissivity_temperature_slope,    init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"lw_emissivity_cloud_slope",          c%lw_emissivity_cloud_slope,          init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"lw_clear_sky_transmissivity",        c%lw_clear_sky_transmissivity,        init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"lw_clear_sky_transmissivity_per_km", c%lw_clear_sky_transmissivity_per_km, init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"lw_night_cloud_fraction",            c%lw_night_cloud_fraction,            init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"eps_snow",                 c%eps_snow,                init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"eps_ice",                  c%eps_ice,                 init=init_pars,defaults_file=def_file,defaults_group=def_const)
        call nml_read(filename,group,"sigma_sb",                 c%sigma_sb,                init=init_pars,defaults_file=def_file,defaults_group=def_const)


        call nml_read(filename,group,"low_density_densification",low_density_densification, init=init_pars,defaults_file=def_file,defaults_group=def_const)

        ! Validate before converting, so the error names the group and the
        ! full allowed set rather than only the offending value.
        call chion_check_enum(group,"albedo_scheme",            albedo_scheme,            CHION_ALBEDO_CHOICES)
        call chion_check_enum(group,"semix_snow_albedo",        semix_snow_albedo,        CHION_SEMIX_SNOW_ALB_CHOICES)
        call chion_check_enum(group,"seb_scheme",               seb_scheme,               CHION_SEB_CHOICES)
        call chion_check_enum(group,"longwave_scheme",          longwave_scheme,          CHION_LONGWAVE_CHOICES)
        call chion_check_enum(group,"turbulent_flux_scheme",    turbulent_flux_scheme,    CHION_TURB_CHOICES)
        call chion_check_enum(group,"climberx_qsat",            climberx_qsat,            CHION_CLIMBERX_QSAT_CHOICES)
        call chion_check_enum(group,"fresh_snow_density_scheme",fresh_snow_density_scheme,CHION_RHOS_CHOICES)
        call chion_check_enum(group,"low_density_densification",low_density_densification,CHION_DENSIFY_CHOICES)

        c%albedo_scheme             = chion_albedo_scheme_flag(albedo_scheme)
        c%semix_snow_albedo         = chion_semix_snow_albedo_flag(semix_snow_albedo)
        c%seb_scheme                = chion_seb_scheme_flag(seb_scheme)
        c%longwave_scheme           = chion_longwave_scheme_flag(longwave_scheme)
        c%turbulent_flux_scheme     = chion_turbulent_flux_scheme_flag(turbulent_flux_scheme)
        c%climberx_qsat             = chion_climberx_qsat_flag(climberx_qsat)
        c%fresh_snow_density_scheme = chion_fresh_snow_density_scheme_flag(fresh_snow_density_scheme)
        c%low_density_densification = chion_densify_scheme_flag(low_density_densification)

        call chion_const_validate(c)

        return

    end subroutine chion_const_load

    subroutine chion_const_from_phys(c,cnst)
        ! Fill the SHARED fields of chion's constants from the program's
        ! fesm-utils record, narrowing to wp. This is the only place they are
        ! set outside chion_const_init (the Chion.jl defaults, used by tests).
        !
        !   chion           phys_const_class
        !   rho_i           rho_ice
        !   rho_w           rho_w
        !   ci              cp_ice
        !   cw              cp_w
        !   Lm              L_ice
        !   grav            g            (CLIMBER-X turbulence only)
        !   T0              T0
        !
        ! Densification's DENSIFY_GRAVITY is deliberately NOT taken from g: it
        ! is part of the reference-reproduction switch (D25).

        implicit none

        type(chion_const_class), intent(INOUT) :: c
        type(phys_const_class),  intent(IN)    :: cnst

        call phys_const_require(cnst,"chion_const_from_phys")

        call phys_const_get(cnst,"rho_ice",c%rho_i)
        call phys_const_get(cnst,"rho_w",  c%rho_w)
        call phys_const_get(cnst,"cp_ice", c%ci)
        call phys_const_get(cnst,"cp_w",   c%cw)
        call phys_const_get(cnst,"L_ice",  c%Lm)
        call phys_const_get(cnst,"g",      c%grav)
        call phys_const_get(cnst,"T0",     c%T0)

        return

    end subroutine chion_const_from_phys

    subroutine bessi_par_load(par,filename,group,init)
        ! The &bessi block. One nml_read per parameter, in the declaration
        ! order of bessi_par_class (snow_bessi.f90:96-114).
        !
        ! Belongs in snow_bessi.f90 next to bessi_par_init; it is here only
        ! because WP11/WP13 may not modify the physics modules. Move it when
        ! that constraint lifts -- the signature is already the yelmo one.

        implicit none

        type(bessi_par_class),      intent(INOUT) :: par
        character(len=*),           intent(IN)    :: filename
        character(len=*),           intent(IN)    :: group      ! usually "bessi"
        logical,          optional, intent(IN)    :: init

        ! Local variables
        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        call nml_validate(filename,def_file,group,defaults_group=def_bessi)

        call nml_read(filename,group,"Ntot",            par%Ntot,            init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"mass_max",        par%mass_max,        init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"mass_split",      par%mass_split,      init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"mass_min",        par%mass_min,        init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"density_init",    par%density_init,    init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"temperature_init",par%temperature_init,init=init_pars,defaults_file=def_file,defaults_group=def_bessi)

        call nml_read(filename,group,"ice_substrate_layers",       par%ice_substrate_layers,       init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"ice_substrate_top_thickness",par%ice_substrate_top_thickness,init=init_pars,defaults_file=def_file,defaults_group=def_bessi)

        call nml_read(filename,group,"near_surface_layer_max_thicknesses",par%near_surface_layer_max_thicknesses,init=init_pars,defaults_file=def_file,defaults_group=def_bessi)

        call nml_read(filename,group,"diurnal_shortwave_substeps",            par%diurnal_shortwave_substeps,            init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_shortwave_threshold",           par%diurnal_shortwave_threshold,           init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_shortwave_max_substeps",        par%diurnal_shortwave_max_substeps,        init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_shortwave_min_air_temperature", par%diurnal_shortwave_min_air_temperature, init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_temperature_cycle",             par%diurnal_temperature_cycle,             init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_temperature_amplitude",         par%diurnal_temperature_amplitude,         init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_temperature_amplitude_gradient",         par%diurnal_temperature_amplitude_gradient,         init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_temperature_amplitude_reference_height", par%diurnal_temperature_amplitude_reference_height, init=init_pars,defaults_file=def_file,defaults_group=def_bessi)
        call nml_read(filename,group,"diurnal_temperature_amplitude_max",              par%diurnal_temperature_amplitude_max,              init=init_pars,defaults_file=def_file,defaults_group=def_bessi)

        call bessi_par_validate(par)

        return

    end subroutine bessi_par_load

    subroutine pdd_par_load(par,filename,group,init)
        ! The &pdd block. One nml_read per parameter, in the declaration order
        ! of pdd_par_class (snow_pdd.f90:56-66).
        !
        ! pdd_method arrives as a string and replaces Chion.jl's former
        ! implicit "27 <= dt_days <= 32" monthly trigger. The default is
        ! "simple", as in Chion.jl; see the defaults file.
        !
        ! Belongs in snow_pdd.f90; see the note on bessi_par_load.

        implicit none

        type(pdd_par_class),        intent(INOUT) :: par
        character(len=*),           intent(IN)    :: filename
        character(len=*),           intent(IN)    :: group      ! usually "pdd"
        logical,          optional, intent(IN)    :: init

        ! Local variables
        logical           :: init_pars
        character(len=56) :: pdd_method

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        call nml_validate(filename,def_file,group,defaults_group=def_pdd)

        call nml_read(filename,group,"ddf_snow",           par%ddf_snow,           init=init_pars,defaults_file=def_file,defaults_group=def_pdd)
        call nml_read(filename,group,"ddf_ice",            par%ddf_ice,            init=init_pars,defaults_file=def_file,defaults_group=def_pdd)
        call nml_read(filename,group,"refreezing_fraction",par%refreezing_fraction,init=init_pars,defaults_file=def_file,defaults_group=def_pdd)
        call nml_read(filename,group,"temperature_sigma",  par%temperature_sigma,  init=init_pars,defaults_file=def_file,defaults_group=def_pdd)
        call nml_read(filename,group,"H_snow_max",         par%H_snow_max,         init=init_pars,defaults_file=def_file,defaults_group=def_pdd)
        call nml_read(filename,group,"pdd_method",         pdd_method,             init=init_pars,defaults_file=def_file,defaults_group=def_pdd)

        call chion_check_enum(group,"pdd_method",pdd_method,CHION_PDD_CHOICES)
        par%pdd_method = pdd_method_flag(pdd_method)

        call pdd_par_validate(par)

        return

    end subroutine pdd_par_load

    subroutine chion_itm_par_load(par,filename,group)
        ! The &itm block, via snow_itm.f90's own itm_par_load (WP12 owns it),
        ! with the schema as defaults: the user's group may be sparse or
        ! absent, exactly as for bessi_par_load and pdd_par_load.

        implicit none

        type(itm_par_class), intent(INOUT) :: par
        character(len=*),    intent(IN)    :: filename
        character(len=*),    intent(IN)    :: group        ! usually "itm"

        call nml_validate(filename,def_file,group,defaults_group=def_itm)

        call itm_par_load(par,filename,init=.TRUE.,group=group, &
                          defaults_file=def_file,defaults_group=def_itm)

        return

    end subroutine chion_itm_par_load

end module chion_api
