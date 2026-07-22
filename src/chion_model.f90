module chion_model
    ! WP11 -- model dispatch and the column loop.
    !
    ! This module owns three things and nothing else:
    !   1. chion_pack_step_forcing -- (ncol) host arrays -> one column's
    !      chion_step_forcing_class, the shared per-column contract.
    !   2. the OpenMP loop over active columns, and the select case that
    !      routes it to one of the three snowpack models.
    !   3. the model-agnostic lifecycle helpers (alloc / init_state / reset /
    !      dealloc / cumulative-SMB extraction) that chion_api calls.
    !
    ! It contains NO physics. Every case in every select block delegates
    ! immediately.
    !
    ! -----------------------------------------------------------------------
    ! Why this module does not see chion_class
    ! -----------------------------------------------------------------------
    ! chion_class is defined in chion_api (docs/porting_notes.md D7). If this
    ! module used chion_api, and chion_api called this module, the two would
    ! be circularly dependent. So every routine here takes the components
    ! (par, c, grd, forc, bsi, pdd, itm) explicitly, and chion_api's public
    ! routines are thin forwarders. That also makes the dispatcher directly
    ! testable without building a full chion_class.
    !
    ! -----------------------------------------------------------------------
    ! The three models have DIFFERENT call shapes. Do not unify them.
    ! -----------------------------------------------------------------------
    !   bessi   call bessi_column_step(bsi,icol,fc,c)
    !               -- takes the whole state object plus a column index,
    !                  because it mutates (Ntot,ncol) layer arrays.
    !   pdd     call pdd_column_step(swe,smb,runoff,pdd_sum,fc,par,c)
    !               -- takes four scalars by reference; it has no layers and
    !                  no state object worth threading.
    !   itm     call itm_step(itm,icol,fc,z_srf,H_ice,PDDs)
    !               -- takes the state object plus three extra arguments that
    !                  are deliberately NOT in chion_step_forcing_class,
    !                  because they are ice-sheet state, not atmosphere.
    ! Forcing a common signature onto these would mean either polluting the
    ! shared forcing type with ice-sheet fields or wrapping PDD in a state
    ! object it does not need. Each case is written out instead.
    !
    ! -----------------------------------------------------------------------
    ! OpenMP correctness
    ! -----------------------------------------------------------------------
    ! The columns are independent by construction: every model writes only to
    ! index icol of every array it touches, and icol is distinct per iteration
    ! because active_idx is a strictly increasing list of distinct column
    ! indices (chion_grid_set_active builds it that way). There is therefore
    ! NO reduction anywhere in this loop, and no critical section is needed.
    !
    ! What must be private is every per-iteration temporary:
    !     i, icol, fc
    ! `fc` is a derived type of scalars with no allocatable components, so
    ! `private(fc)` gives each thread its own complete copy. Everything else
    ! (bsi, pdd, itm, forc, par, c) is shared and read-or-indexed-by-icol only.
    !
    ! default(shared) is stated explicitly rather than relied upon, so that a
    ! newly added local cannot silently become shared.

    use chion_defs, only : wp, wp_acc, io_unit_err, &
                           chion_const_class, chion_param_class, &
                           chion_grid_class, chion_forcing_class, &
                           chion_step_forcing_class

    use snow_bessi, only : bessi_class, bessi_alloc, bessi_dealloc, &
                           bessi_init_state, bessi_reset_columns, &
                           bessi_column_step
    use snow_pdd,   only : pdd_class, pdd_alloc, pdd_dealloc, &
                           pdd_init_state, pdd_reset_column, &
                           pdd_column_step
    use snow_itm,   only : itm_class, itm_alloc, itm_dealloc, &
                           itm_init_state, itm_step

    implicit none

    private

    ! The allowed values of par%model, in one place. Used by the dispatcher's
    ! error messages and by chion_api's enum validation, so the two can never
    ! disagree.
    character(len=*), parameter, public :: CHION_MODEL_CHOICES = "bessi|pdd|itm"

    public :: chion_pack_step_forcing
    public :: chion_model_alloc
    public :: chion_model_dealloc
    public :: chion_model_init_state
    public :: chion_model_reset_columns
    public :: chion_model_step
    public :: chion_model_smb_cum
    public :: chion_model_summary_line

contains

    ! =====================================================================
    ! Forcing pack
    ! =====================================================================

    subroutine chion_pack_step_forcing(forc,icol,dt_days,fc)
        ! Fill the per-column, per-step forcing contract from the host-facing
        ! (ncol) arrays. Chion.jl _step_forcing (src/forcing.jl:259-300).
        !
        ! Fields are assigned in declaration order, which is also the order of
        ! Chion.jl's SnowpackStepForcing, so the two can be diffed by eye.
        !
        ! NOTE what is NOT here: surface_height, H_ice and PDDs. The first is
        ! used by the host to derive air_pressure; the last two are ITM's, and
        ! are passed to itm_step directly (see the module header).

        implicit none

        type(chion_forcing_class),      intent(IN)  :: forc
        integer,                        intent(IN)  :: icol
        real(wp),                       intent(IN)  :: dt_days
        type(chion_step_forcing_class), intent(OUT) :: fc

        fc%air_temperature = forc%air_temperature(icol)
        fc%dt_days         = dt_days
        fc%snowfall_rate   = forc%snowfall_rate(icol)
        fc%rainfall_rate   = forc%rainfall_rate(icol)
        fc%shortwave_down  = forc%shortwave_down(icol)
        fc%wind_speed      = forc%wind_speed(icol)

        fc%q_sw_net  = forc%q_sw_net(icol)
        fc%q_lw_down = forc%q_lw_down(icol)
        fc%q_sh      = forc%q_sh(icol)
        fc%q_lh      = forc%q_lh(icol)

        fc%has_q_sw_net  = forc%has_q_sw_net(icol)
        fc%has_q_lw_down = forc%has_q_lw_down(icol)
        fc%has_q_sh      = forc%has_q_sh(icol)
        fc%has_q_lh      = forc%has_q_lh(icol)

        fc%relative_humidity     = forc%relative_humidity(icol)
        fc%has_relative_humidity = forc%has_relative_humidity(icol)

        fc%air_pressure          = forc%air_pressure(icol)
        fc%prescribed_albedo     = forc%prescribed_albedo(icol)
        fc%has_prescribed_albedo = forc%has_prescribed_albedo(icol)

        fc%latitude_deg        = forc%latitude_deg(icol)
        fc%day_of_year         = forc%day_of_year
        fc%solar_longitude_deg = forc%solar_longitude_deg

        return

    end subroutine chion_pack_step_forcing

    ! =====================================================================
    ! Lifecycle dispatch
    ! =====================================================================

    subroutine chion_model_alloc(par,bsi,pdd,itm,ncol)
        ! Allocate the state of the SELECTED model only. The other two state
        ! objects are left unallocated, which is what makes
        ! `allocated(chn%bsi%now%mass)` a meaningful test of which model is
        ! running, and what keeps a 20,000-column BESSI domain from also
        ! carrying PDD and ITM arrays it never touches.

        implicit none

        type(chion_param_class), intent(IN)    :: par
        type(bessi_class),       intent(INOUT) :: bsi
        type(pdd_class),         intent(INOUT) :: pdd
        type(itm_class),         intent(INOUT) :: itm
        integer,                 intent(IN)    :: ncol

        select case(trim(par%model))
            case("bessi")
                call bessi_alloc(bsi,ncol)
            case("pdd")
                call pdd_alloc(pdd,ncol)
            case("itm")
                call itm_alloc(itm,ncol)
            case DEFAULT
                call chion_model_error("chion_model_alloc",par%model)
        end select

        return

    end subroutine chion_model_alloc

    subroutine chion_model_dealloc(bsi,pdd,itm)
        ! Deallocate all three unconditionally. Each *_dealloc guards every
        ! array with allocated(), so this is safe whichever model ran, and it
        ! cannot leak if par%model was changed between init and end.

        implicit none

        type(bessi_class), intent(INOUT) :: bsi
        type(pdd_class),   intent(INOUT) :: pdd
        type(itm_class),   intent(INOUT) :: itm

        call bessi_dealloc(bsi)
        call pdd_dealloc(pdd)
        call itm_dealloc(itm)

        return

    end subroutine chion_model_dealloc

    subroutine chion_model_init_state(par,c,bsi,pdd,itm)
        ! Cold start. Restart reading is WP14.

        implicit none

        type(chion_param_class), intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        type(bessi_class),       intent(INOUT) :: bsi
        type(pdd_class),         intent(INOUT) :: pdd
        type(itm_class),         intent(INOUT) :: itm

        select case(trim(par%model))
            case("bessi")
                call bessi_init_state(bsi,c)
            case("pdd")
                call pdd_init_state(pdd)
            case("itm")
                ! No H_snow argument: itm_init_state then seeds the snowpack
                ! at H_snow_max, matching a fresh smbpal run (smbpal.f90:117).
                call itm_init_state(itm)
            case DEFAULT
                call chion_model_error("chion_model_init_state",par%model)
        end select

        return

    end subroutine chion_model_init_state

    subroutine chion_model_reset_columns(par,c,bsi,pdd,itm,idx)
        ! Reset the listed columns to their cold-start state. Called by
        ! chion_set_active_mask for columns that have just been switched OFF,
        ! mirroring Chion.jl _reset_model_columns! (src/runtime.jl:98-175).
        !
        ! Resetting on DEACTIVATION rather than on reactivation is Chion.jl's
        ! choice and is preserved: it means a column that is switched off
        ! immediately stops reporting stale accumulators to the output, and a
        ! column that is switched back on necessarily starts clean.

        implicit none

        type(chion_param_class), intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        type(bessi_class),       intent(INOUT) :: bsi
        type(pdd_class),         intent(INOUT) :: pdd
        type(itm_class),         intent(INOUT) :: itm
        integer,                 intent(IN)    :: idx(:)

        ! Local variables
        integer :: i

        if (size(idx) .eq. 0) return

        select case(trim(par%model))

            case("bessi")
                call bessi_reset_columns(bsi,c,idx)

            case("pdd")
                do i = 1, size(idx)
                    call pdd_reset_column(pdd,idx(i))
                end do

            case("itm")
                ! snow_itm.f90 has no per-column reset routine (its smbpal
                ! ancestor has no active-mask concept at all), so the
                ! cold-start values of itm_init_state are applied here, per
                ! column. Kept identical to itm_init_state field for field --
                ! if that routine's cold start changes, change this too.
                do i = 1, size(idx)
                    call itm_reset_column(itm,idx(i))
                end do

            case DEFAULT
                call chion_model_error("chion_model_reset_columns",par%model)

        end select

        return

    end subroutine chion_model_reset_columns

    subroutine itm_reset_column(itm,icol)
        ! Cold-start values for one ITM column. Mirrors itm_init_state
        ! (snow_itm.f90:301-346) exactly, including the H_snow_max seed and
        ! the hard-coded 273.15 for tsrf.

        implicit none

        type(itm_class), intent(INOUT) :: itm
        integer,         intent(IN)    :: icol

        if (icol .lt. 1 .or. icol .gt. itm%now%ncol) then
            write(io_unit_err,*) "itm_reset_column:: Error: column index out of range."
            write(io_unit_err,*) "icol, ncol = ", icol, itm%now%ncol
            stop "Program stopped."
        end if

        itm%now%H_snow(icol)   = itm%par%H_snow_max
        itm%now%alb_s(icol)    = itm%par%alb_snow_dry
        itm%now%smb(icol)      = 0.0_wp
        itm%now%smbi(icol)     = 0.0_wp
        itm%now%melt(icol)     = 0.0_wp
        itm%now%runoff(icol)   = 0.0_wp
        itm%now%refrz(icol)    = 0.0_wp
        itm%now%melt_net(icol) = 0.0_wp
        itm%now%tsrf(icol)     = 273.15_wp

        itm%now%smb_cum(icol)    = 0.0_wp_acc
        itm%now%smbi_cum(icol)   = 0.0_wp_acc
        itm%now%melt_cum(icol)   = 0.0_wp_acc
        itm%now%runoff_cum(icol) = 0.0_wp_acc
        itm%now%refrz_cum(icol)  = 0.0_wp_acc

        return

    end subroutine itm_reset_column

    ! =====================================================================
    ! The column loop
    ! =====================================================================

    subroutine chion_model_step(par,c,grd,forc,bsi,pdd,itm,dt_days)
        ! Advance every ACTIVE column by exactly one step of dt_days.
        !
        ! Deviation from yelmo, deliberate (docs/PLAN.md section 2.1): this
        ! takes dt_days and advances one step, rather than taking a target
        ! time and sub-cycling. The host owns the time loop. The only internal
        ! sub-cycling is BESSI's diurnal substepping, inside bessi_column_step.

        implicit none

        type(chion_param_class),   intent(IN)    :: par
        type(chion_const_class),   intent(IN)    :: c
        type(chion_grid_class),    intent(IN)    :: grd
        type(chion_forcing_class), intent(IN)    :: forc
        type(bessi_class),         intent(INOUT) :: bsi
        type(pdd_class),           intent(INOUT) :: pdd
        type(itm_class),           intent(INOUT) :: itm
        real(wp),                  intent(IN)    :: dt_days

        ! Local variables
        integer :: i, icol
        type(chion_step_forcing_class) :: fc

        if (forc%ncol .ne. grd%ncol) then
            write(io_unit_err,*) "chion_model_step:: Error: forcing and grid column counts differ."
            write(io_unit_err,*) "forc%ncol, grd%ncol = ", forc%ncol, grd%ncol
            stop "Program stopped."
        end if

        if (dt_days .le. 0.0_wp) then
            write(io_unit_err,*) "chion_model_step:: Error: dt_days must be positive."
            write(io_unit_err,*) "dt_days = ", dt_days
            stop "Program stopped."
        end if

        select case(trim(par%model))

            case("bessi")

                !$omp parallel do default(shared) private(i,icol,fc)
                do i = 1, grd%n_active
                    icol = grd%active_idx(i)
                    call chion_pack_step_forcing(forc,icol,dt_days,fc)
                    call bessi_column_step(bsi,icol,fc,c)
                end do
                !$omp end parallel do

            case("pdd")

                !$omp parallel do default(shared) private(i,icol,fc)
                do i = 1, grd%n_active
                    icol = grd%active_idx(i)
                    call chion_pack_step_forcing(forc,icol,dt_days,fc)
                    call pdd_column_step(pdd%now%snowpack_swe(icol), &
                                         pdd%now%smb_ice(icol),      &
                                         pdd%now%runoff(icol),       &
                                         pdd%now%pdd_sum(icol),      &
                                         fc,pdd%par,c)
                end do
                !$omp end parallel do

            case("itm")

                !$omp parallel do default(shared) private(i,icol,fc)
                do i = 1, grd%n_active
                    icol = grd%active_idx(i)
                    call chion_pack_step_forcing(forc,icol,dt_days,fc)
                    call itm_step(itm,icol,fc,forc%surface_height(icol), &
                                  forc%H_ice(icol),forc%PDDs(icol))
                end do
                !$omp end parallel do

            case DEFAULT

                call chion_model_error("chion_model_step",par%model)

        end select

        return

    end subroutine chion_model_step

    ! =====================================================================
    ! Model-agnostic SMB
    ! =====================================================================

    subroutine chion_model_smb_cum(par,bsi,pdd,itm,smb_cum)
        ! Cumulative NET MASS FLUX TO THE ICE SHEET, per column, since the
        ! last cold start or column reset. Units [kg m-2], positive = the ice
        ! sheet gains mass. This is the raw material for chion_get_smb; see
        ! that routine's header for the full definition and the reasoning.
        !
        ! Each model's native output is mapped as follows.
        !
        ! BESSI -- smb_ice is used unchanged.
        !   BESSI credits smb_ice with exactly the terms that cross the
        !   snowpack/ice interface: bottom export (mass_base), bare-ice
        !   deposition/sublimation, bare-surface melt and ice melt. Snow that
        !   merely sits in the pack is not in it. That is already the
        !   ice-sheet-facing quantity.
        !
        ! PDD -- smb_ice is used unchanged.
        !   chion's PDD is ice-facing by construction: it accumulates
        !   snow_to_ice - ice_melt (docs/porting_notes.md D23). No reservoir
        !   term to remove.
        !
        !   This REPLACES the earlier `smb_ice - snowpack_swe` accessor, which
        !   existed only to undo Chion.jl's whole-column convention at the API
        !   boundary. That reconciliation had two costs, both now gone: it
        !   differenced an sp-stored reservoir, so the recovered flux carried
        !   round-off bounded by eps_sp*snowpack_swe per step; and because PDD
        !   had no snow-to-ice pathway, the result was identically -ice_melt
        !   and therefore structurally never positive. With the capped
        !   reservoir, excess snow above H_snow_max converts to ice, so PDD
        !   reports a genuine positive flux in the accumulation zone.
        !
        ! ITM -- smbi_cum is used unchanged.
        !   smbpal's smbi is "the mass balance seen by the ice sheet":
        !   snow_to_ice + refrz - melted_ice (smb_itm.f90:187). Note this is
        !   NOT smbpal's smb, which is sf + rf - runoff, the whole-column
        !   balance -- the same distinction PDD gets wrong. Units are
        !   mm w.e., which is kg m-2 by definition.
        !
        ! The three are therefore the same physical quantity in the same
        ! units, and differ only in which processes each model resolves.

        implicit none

        type(chion_param_class), intent(IN)  :: par
        type(bessi_class),       intent(IN)  :: bsi
        type(pdd_class),         intent(IN)  :: pdd
        type(itm_class),         intent(IN)  :: itm
        real(wp_acc),            intent(OUT) :: smb_cum(:)

        select case(trim(par%model))

            case("bessi")
                if (size(smb_cum) .ne. bsi%now%ncol) call chion_size_error(size(smb_cum),bsi%now%ncol)
                smb_cum = bsi%now%smb_ice

            case("pdd")
                if (size(smb_cum) .ne. pdd%now%ncol) call chion_size_error(size(smb_cum),pdd%now%ncol)
                smb_cum = pdd%now%smb_ice

            case("itm")
                if (size(smb_cum) .ne. itm%now%ncol) call chion_size_error(size(smb_cum),itm%now%ncol)
                smb_cum = itm%now%smbi_cum

            case DEFAULT
                call chion_model_error("chion_model_smb_cum",par%model)

        end select

        return

    end subroutine chion_model_smb_cum

    ! =====================================================================
    ! Provenance
    ! =====================================================================

    subroutine chion_model_summary_line(par,bsi,pdd,itm,line)
        ! One line describing the selected model's configuration, for run
        ! logs. Kept here rather than in chion_api so that adding a model
        ! means editing exactly one file.

        implicit none

        type(chion_param_class), intent(IN)  :: par
        type(bessi_class),       intent(IN)  :: bsi
        type(pdd_class),         intent(IN)  :: pdd
        type(itm_class),         intent(IN)  :: itm
        character(len=*),        intent(OUT) :: line

        select case(trim(par%model))

            case("bessi")
                write(line,"(a,i0,a,f7.1,a,f7.1,a,f7.1,a,l1,a,i0)")            &
                    "bessi  : Ntot=", bsi%par%Ntot,                            &
                    "  mass_max=", bsi%par%mass_max,                           &
                    "  mass_split=", bsi%par%mass_split,                       &
                    "  mass_min=", bsi%par%mass_min,                           &
                    "  diurnal=", bsi%par%diurnal_shortwave_substeps,          &
                    "  max_substeps=", bsi%par%diurnal_shortwave_max_substeps

            case("pdd")
                write(line,"(a,f7.3,a,f7.3,a,f6.3,a,f6.3,a,i0)")               &
                    "pdd    : ddf_snow=", pdd%par%ddf_snow,                    &
                    "  ddf_ice=", pdd%par%ddf_ice,                             &
                    "  refreezing_fraction=", pdd%par%refreezing_fraction,     &
                    "  sigma=", pdd%par%temperature_sigma,                     &
                    "  pdd_method=", pdd%par%pdd_method

            case("itm")
                write(line,"(a,f8.2,a,f7.3,a,f8.1,a,f6.3)")                    &
                    "itm    : itm_c=", itm%par%itm_c,                          &
                    "  itm_t=", itm%par%itm_t,                                 &
                    "  H_snow_max=", itm%par%H_snow_max,                       &
                    "  Pmaxfrac=", itm%par%Pmaxfrac

            case DEFAULT
                call chion_model_error("chion_model_summary_line",par%model)

        end select

        return

    end subroutine chion_model_summary_line

    ! =====================================================================
    ! Errors
    ! =====================================================================

    subroutine chion_model_error(routine,model)
        ! The mandatory case DEFAULT, per the yelmo idiom: name the routine,
        ! name the offending value, and list the allowed set. Every select
        ! case in this module routes its default here so the message can
        ! never drift between them.

        implicit none

        character(len=*), intent(IN) :: routine
        character(len=*), intent(IN) :: model

        write(io_unit_err,*) trim(routine)//":: Error: model not recognized."
        write(io_unit_err,*) "model should be one of: ['bessi','pdd','itm']"
        write(io_unit_err,*) "model = ", trim(model)
        stop "Program stopped."

        return

    end subroutine chion_model_error

    subroutine chion_size_error(n_given,n_expect)

        implicit none

        integer, intent(IN) :: n_given
        integer, intent(IN) :: n_expect

        write(io_unit_err,*) "chion_model_smb_cum:: Error: smb array has the wrong length."
        write(io_unit_err,*) "size(smb), ncol = ", n_given, n_expect
        stop "Program stopped."

        return

    end subroutine chion_size_error

end module chion_model
