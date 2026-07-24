module snow_bessi
    ! WP8 -- BESSI model assembly.
    !
    ! Fortran port of Chion.jl (branch main):
    !   src/models.jl:11-78      BESSIModel        -> bessi_par_class
    !   src/state.jl:11-134      BESSIState        -> bessi_state_class
    !   src/runtime.jl:98-175    column reset      -> bessi_reset_columns
    !   src/step.jl:159-407      column_step_core! -> bessi_column_step_core
    !   src/step.jl:75-157       column_step!      -> bessi_column_step
    !   src/domain.jl:59-64      _validate_mass_partition -> bessi_par_validate
    !
    ! This module contains NO physics of its own. Every physical operation is
    ! delegated to the Level-1 kernels (WP4-WP7). What it owns is the ORDER of
    ! operations, which is the part of BESSI that is easy to get subtly wrong.
    !
    ! ---------------------------------------------------------------------
    ! ORDER OF OPERATIONS (docs/PLAN.md section 1) -- do not reorder
    ! ---------------------------------------------------------------------
    !   1  dt_seconds; snapshot started_without_surface_snow BEFORE anything
    !      else; resolve use_prescribed_albedo
    !   2  accumulation (snow + rain, split/merge, depth cap)
    !   3  fresh snow onto a column that started bare -> temperature(1) = T_air
    !   4  prescribed albedo, if any, applied BEFORE the bare-surface test
    !   5  no surface snow -> bare-ice ablation, accumulate diagnostics, RETURN
    !      (percolation and refreezing are skipped ENTIRELY -- trap, see below)
    !   6  albedo update (prescribed | dynamic-or-constant)
    !   7  HTESSEL only: snapshot mass_w(1:n) BEFORE the energy solve
    !   8  accumulation_rate -> densification
    !   9  latent-heat coefficients -> implicit energy solve
    !  10  post-solve surface vapor mass flux; accumulate its three outputs
    !  11  melt, if the solve asked for it; ice melt if the snow ran out
    !  12  percolation -> runoff
    !  13  HTESSEL only: liquid-water compaction, using the step-7 snapshot
    !  14  refreezing
    !  15  final albedo fixup
    !
    ! ---------------------------------------------------------------------
    ! Traps honoured here (docs/PLAN.md section 5)
    ! ---------------------------------------------------------------------
    ! * Step 5 is an EARLY RETURN. A column that is bare after accumulation
    !   does not percolate and does not refreeze, even if it holds liquid
    !   water (it can: rain falls into mass_w(1) whenever mass(1) > 0, which
    !   is a weaker test than surface_has_snow). Water simply sits there until
    !   the column has snow again.
    ! * Step 11 accumulates the FULL REQUESTED melt mass into the melt
    !   diagnostic, not the amount the snowpack could actually supply. The
    !   shortfall is separately charged to smb_ice and runoff as ice melt, so
    !   melt is "melt demand", not "snow melted".
    ! * Step 10 accumulates latent_heat_flux_sum in units of W m-2 DAYS:
    !   the multiplier is forc%dt_days, NOT dt_seconds. Same in step 5.
    ! * The HTESSEL snapshot is taken before the energy solve and consumed
    !   after percolation -- four kernels later.
    ! * Fortran does not short-circuit .and./.or. (docs/porting_notes.md D10),
    !   so every compound guard whose second operand is only valid under the
    !   first is written as a nested if.
    !
    ! Accumulators are wp_acc (dp) and are incremented IN PLACE. Never route
    ! an increment through a wp temporary -- that is the one mistake in this
    ! module that would compile, run, and only show up as slow mass drift.

    use chion_defs, only : wp, wp_acc, io_unit_err, &
                           TOL_TINY, TOL_EMPTY_LAYER, &
                           DEF_NTOT, DEF_MASS_MAX, DEF_MASS_SPLIT, DEF_MASS_MIN, &
                           DEF_DENSITY_INIT, DEF_TEMPERATURE_INIT, &
                           CHION_ALBEDO_PRESCRIBED, CHION_ALBEDO_SEMIX, &
                           CHION_DENSIFY_HTESSEL, &
                           chion_const_class, chion_step_forcing_class

    use snow_column_utils,  only : surface_has_snow, column_has_liquid_water

    use snow_accumulation,  only : apply_accumulation
    use snow_albedo,        only : albedo_update
    use snow_albedo_semix,  only : semix_surface_albedo, semix_dust_concentration, &
                                   semix_daily_coszm
    use snow_densify,       only : densify_column, apply_htessel_liquid_water_compaction
    use snow_energy,        only : snow_energy_flux, snow_energy_result_class
    use snow_surface_fluxes,only : bare_ice_ablation_class, bare_ice_ablation_mass, &
                                   latent_heat_coeff_class, &
                                   diagnose_latent_heat_flux_coefficients, &
                                   surface_vapor_flux_class, &
                                   apply_snow_surface_vapor_mass_flux
    use snow_melt,          only : apply_melt
    use snow_percolation,   only : apply_percolation
    use snow_refreezing,    only : apply_refreezing
    use snow_diurnal,       only : diurnal_substep_count, diurnal_substep_bounds, &
                                   diurnal_shortwave_interval_average, &
                                   diurnal_temperature_interval_average

    implicit none

    private

    real(wp_acc), parameter :: PI_ACC = 3.141592653589793238462643383279502884_wp_acc

    ! === Types ===============================================================

    type bessi_par_class
        ! Chion.jl BESSIModel (src/models.jl:11-26), minus the grid and the
        ! constants object, which chion keeps separately (chion_grid_class and
        ! chion_const_class). The six diurnal fields are Chion.jl's step
        ! `config` NamedTuple (src/step.jl:45-52), which is built directly from
        ! the model (src/simulation.jl:172-180), so they belong here.

        integer  :: Ntot                  ! [1] layer capacity
        real(wp) :: mass_max              ! [kg m-2] surface-layer split trigger
        real(wp) :: mass_split            ! [kg m-2] mass left below after a split
        real(wp) :: mass_min              ! [kg m-2] surface-layer merge trigger
        real(wp) :: density_init          ! [kg m-3] density of an inactive layer
        real(wp) :: temperature_init      ! [K] temperature of an inactive layer

        logical  :: diurnal_shortwave_substeps            ! [1] enable substepping
        real(wp) :: diurnal_shortwave_threshold           ! [W m-2] peak-minus-mean excess
        integer  :: diurnal_shortwave_max_substeps        ! [1] 1..24
        real(wp) :: diurnal_shortwave_min_air_temperature ! [K]
        logical  :: diurnal_temperature_cycle             ! [1] impose a diurnal T cycle
        real(wp) :: diurnal_temperature_amplitude         ! [K] half-amplitude
    end type bessi_par_class

    type bessi_state_class
        ! docs/PLAN.md section 2.2. Layer arrays are (Ntot,ncol), layer-major,
        ! layer 1 = surface. Byte-identical to Chion.jl's Matrix{NF} layout.

        integer :: ncol
        integer :: Ntot

        integer,  allocatable :: n_lay(:)          ! (ncol) active layer count
        real(wp), allocatable :: mass(:,:)         ! (Ntot,ncol) [kg m-2] solid
        real(wp), allocatable :: mass_w(:,:)       ! (Ntot,ncol) [kg m-2] liquid
        real(wp), allocatable :: density(:,:)      ! (Ntot,ncol) [kg m-3]
        real(wp), allocatable :: temperature(:,:)  ! (Ntot,ncol) [K]

        ! Cumulative accumulators -- MUST be wp_acc. docs/PLAN.md section 3.1.
        real(wp_acc), allocatable :: mass_base(:)             ! [kg m-2]
        real(wp_acc), allocatable :: smb_ice(:)               ! [kg m-2]
        real(wp_acc), allocatable :: runoff(:)                ! [kg m-2]
        real(wp_acc), allocatable :: melt(:)                  ! [kg m-2]
        real(wp_acc), allocatable :: refreezing(:)            ! [kg m-2]
        real(wp_acc), allocatable :: vapor_mass(:)            ! [kg m-2] + = deposition
        real(wp_acc), allocatable :: sublimation(:)           ! [kg m-2] >= 0
        real(wp_acc), allocatable :: latent_heat_flux_sum(:)  ! [W m-2 d]

        ! Instantaneous per-column scalars
        real(wp), allocatable :: t_srf(:)          ! (ncol) [K]
        real(wp), allocatable :: albedo(:)         ! (ncol) [1]

        ! SEMIX albedo bookkeeping: column SWE last step and its seasonal peak.
        ! The drawdown (w_snow_max - w_snow) drives dust concentration.
        real(wp), allocatable :: w_snow_old(:)     ! (ncol) [kg m-2]
        real(wp), allocatable :: w_snow_max(:)     ! (ncol) [kg m-2]

        ! Diagnostics, filled by snow_diagnostics (WP10)
        real(wp), allocatable :: thickness(:)      ! (ncol) [m]
        real(wp), allocatable :: wet_mass(:)       ! (ncol) [kg m-2]
        real(wp), allocatable :: bulk_density(:)   ! (ncol) [kg m-3]
        real(wp), allocatable :: liquid_water(:)   ! (ncol) [1]
    end type bessi_state_class

    type bessi_class
        type(bessi_par_class)   :: par
        type(bessi_state_class) :: now
    end type bessi_class

    public :: bessi_par_class
    public :: bessi_state_class
    public :: bessi_class

    public :: bessi_par_init
    public :: bessi_par_validate
    public :: bessi_alloc
    public :: bessi_dealloc
    public :: bessi_init_state
    public :: bessi_reset_columns
    public :: bessi_column_step_core
    public :: bessi_column_step

contains

    ! =====================================================================
    ! Parameters
    ! =====================================================================

    subroutine bessi_par_init(par)
        ! Chion.jl defaults: src/models.jl:28-78 (BESSIModel keyword defaults)
        ! and src/constants.jl (DEFAULT_NTOT etc., re-exported by chion_defs).
        !
        ! NOTE the two unit conversions Julia performs in the constructor:
        !   diurnal_shortwave_min_air_temperature_c = -8.0  ->  265.15 K
        !   diurnal_temperature_amplitude_c         =  5.0  ->    5.0 K
        ! The amplitude is a temperature DIFFERENCE, so only the minimum air
        ! temperature gains the 273.15 offset.
        !
        ! WP13 replaces this with bessi_par_load reading a namelist; these stay
        ! as the defaults.

        implicit none

        type(bessi_par_class), intent(OUT) :: par

        par%Ntot             = DEF_NTOT
        par%mass_max         = DEF_MASS_MAX
        par%mass_split       = DEF_MASS_SPLIT
        par%mass_min         = DEF_MASS_MIN
        par%density_init     = DEF_DENSITY_INIT
        par%temperature_init = DEF_TEMPERATURE_INIT

        par%diurnal_shortwave_substeps            = .FALSE.
        par%diurnal_shortwave_threshold           = 0.0_wp
        par%diurnal_shortwave_max_substeps        = 3
        par%diurnal_shortwave_min_air_temperature = 265.15_wp
        par%diurnal_temperature_cycle             = .FALSE.
        par%diurnal_temperature_amplitude         = 5.0_wp

        call bessi_par_validate(par)

        return

    end subroutine bessi_par_init

    subroutine bessi_par_validate(par)
        ! Chion.jl _validate_mass_partition (src/domain.jl:59-64) plus the four
        ! diurnal guards in the BESSIModel constructor (src/models.jl:48-51).
        !
        ! The mass partition is not merely cosmetic:
        !   mass_split < mass_max   -- otherwise the split loop cannot converge
        !   mass_min   < mass_split -- otherwise a split immediately re-merges
        !   mass_split/mass_max >= 0.5 -- so the layer left behind by a split
        !                                 holds at least half of mass_max, and
        !                                 merge_surface_layer's unguarded
        !                                 divisor (2*mass_split - mass_min)
        !                                 stays positive. See upstream defect 8
        !                                 in docs/porting_notes.md.

        implicit none

        type(bessi_par_class), intent(IN) :: par

        if (par%Ntot .lt. 1) then
            write(io_unit_err,*) "bessi_par_validate:: Error: Ntot must be at least 1."
            write(io_unit_err,*) "Ntot = ", par%Ntot
            stop "Program stopped."
        end if

        if (par%mass_max .le. 0.0_wp) then
            write(io_unit_err,*) "bessi_par_validate:: Error: mass_max must be positive."
            write(io_unit_err,*) "mass_max = ", par%mass_max
            stop "Program stopped."
        end if

        if (.not. (par%mass_split .lt. par%mass_max)) then
            write(io_unit_err,*) "bessi_par_validate:: Error: mass_split must be smaller than mass_max."
            write(io_unit_err,*) "mass_split, mass_max = ", par%mass_split, par%mass_max
            stop "Program stopped."
        end if

        if (.not. (par%mass_min .lt. par%mass_split)) then
            write(io_unit_err,*) "bessi_par_validate:: Error: mass_min must be smaller than mass_split."
            write(io_unit_err,*) "mass_min, mass_split = ", par%mass_min, par%mass_split
            stop "Program stopped."
        end if

        if (.not. (par%mass_split/par%mass_max .ge. 0.5_wp)) then
            write(io_unit_err,*) "bessi_par_validate:: Error: mass_split/mass_max must be at least 0.5."
            write(io_unit_err,*) "mass_split, mass_max, ratio = ", &
                                 par%mass_split, par%mass_max, par%mass_split/par%mass_max
            stop "Program stopped."
        end if

        if (par%diurnal_shortwave_threshold .lt. 0.0_wp) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &diurnal_shortwave_threshold must be non-negative."
            write(io_unit_err,*) "diurnal_shortwave_threshold = ", par%diurnal_shortwave_threshold
            stop "Program stopped."
        end if

        if (par%diurnal_shortwave_max_substeps .lt. 1 .or. &
            par%diurnal_shortwave_max_substeps .gt. 24) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &diurnal_shortwave_max_substeps must be between 1 and 24."
            write(io_unit_err,*) "diurnal_shortwave_max_substeps = ", &
                                 par%diurnal_shortwave_max_substeps
            stop "Program stopped."
        end if

        if (par%diurnal_temperature_amplitude .lt. 0.0_wp) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &diurnal_temperature_amplitude must be non-negative."
            write(io_unit_err,*) "diurnal_temperature_amplitude = ", &
                                 par%diurnal_temperature_amplitude
            stop "Program stopped."
        end if

        return

    end subroutine bessi_par_validate

    ! =====================================================================
    ! Allocation
    ! =====================================================================

    subroutine bessi_alloc(bsi,ncol)
        ! Allocate the state arrays. Computes nothing: bessi_init_state does
        ! that, following the yelmo _init / _init_state split.

        implicit none

        type(bessi_class), intent(INOUT) :: bsi
        integer,           intent(IN)    :: ncol

        call bessi_dealloc(bsi)

        call bessi_par_validate(bsi%par)

        if (ncol .lt. 1) then
            write(io_unit_err,*) "bessi_alloc:: Error: ncol must be positive."
            write(io_unit_err,*) "ncol = ", ncol
            stop "Program stopped."
        end if

        bsi%now%ncol = ncol
        bsi%now%Ntot = bsi%par%Ntot

        allocate(bsi%now%n_lay(ncol))
        allocate(bsi%now%mass(bsi%par%Ntot,ncol))
        allocate(bsi%now%mass_w(bsi%par%Ntot,ncol))
        allocate(bsi%now%density(bsi%par%Ntot,ncol))
        allocate(bsi%now%temperature(bsi%par%Ntot,ncol))

        allocate(bsi%now%mass_base(ncol))
        allocate(bsi%now%smb_ice(ncol))
        allocate(bsi%now%runoff(ncol))
        allocate(bsi%now%melt(ncol))
        allocate(bsi%now%refreezing(ncol))
        allocate(bsi%now%vapor_mass(ncol))
        allocate(bsi%now%sublimation(ncol))
        allocate(bsi%now%latent_heat_flux_sum(ncol))

        allocate(bsi%now%t_srf(ncol))
        allocate(bsi%now%albedo(ncol))
        allocate(bsi%now%w_snow_old(ncol))
        allocate(bsi%now%w_snow_max(ncol))

        allocate(bsi%now%thickness(ncol))
        allocate(bsi%now%wet_mass(ncol))
        allocate(bsi%now%bulk_density(ncol))
        allocate(bsi%now%liquid_water(ncol))

        return

    end subroutine bessi_alloc

    subroutine bessi_dealloc(bsi)

        implicit none

        type(bessi_class), intent(INOUT) :: bsi

        if (allocated(bsi%now%n_lay))       deallocate(bsi%now%n_lay)
        if (allocated(bsi%now%mass))        deallocate(bsi%now%mass)
        if (allocated(bsi%now%mass_w))      deallocate(bsi%now%mass_w)
        if (allocated(bsi%now%density))     deallocate(bsi%now%density)
        if (allocated(bsi%now%temperature)) deallocate(bsi%now%temperature)

        if (allocated(bsi%now%mass_base))            deallocate(bsi%now%mass_base)
        if (allocated(bsi%now%smb_ice))              deallocate(bsi%now%smb_ice)
        if (allocated(bsi%now%runoff))               deallocate(bsi%now%runoff)
        if (allocated(bsi%now%melt))                 deallocate(bsi%now%melt)
        if (allocated(bsi%now%refreezing))           deallocate(bsi%now%refreezing)
        if (allocated(bsi%now%vapor_mass))           deallocate(bsi%now%vapor_mass)
        if (allocated(bsi%now%sublimation))          deallocate(bsi%now%sublimation)
        if (allocated(bsi%now%latent_heat_flux_sum)) deallocate(bsi%now%latent_heat_flux_sum)

        if (allocated(bsi%now%t_srf))  deallocate(bsi%now%t_srf)
        if (allocated(bsi%now%albedo)) deallocate(bsi%now%albedo)
        if (allocated(bsi%now%w_snow_old)) deallocate(bsi%now%w_snow_old)
        if (allocated(bsi%now%w_snow_max)) deallocate(bsi%now%w_snow_max)

        if (allocated(bsi%now%thickness))    deallocate(bsi%now%thickness)
        if (allocated(bsi%now%wet_mass))     deallocate(bsi%now%wet_mass)
        if (allocated(bsi%now%bulk_density)) deallocate(bsi%now%bulk_density)
        if (allocated(bsi%now%liquid_water)) deallocate(bsi%now%liquid_water)

        bsi%now%ncol = 0
        bsi%now%Ntot = 0

        return

    end subroutine bessi_dealloc

    ! =====================================================================
    ! State initialization and reset
    ! =====================================================================

    subroutine bessi_init_state(bsi,c)
        ! Chion.jl _initialize_bessi_state_kernel! (src/state.jl:42-90).
        !
        ! A cold start has NO active layers (n = 0). The layer arrays are still
        ! filled with density_init / temperature_init rather than zero, because
        ! reset_layer_at_index restores exactly those values when a layer is
        ! retired, so an inactive slot always looks the same however it got
        ! there.
        !
        ! t_srf starts at T0 and the albedo at alpha_dry -- note NOT alpha_ice,
        ! even though the column is bare. The first step overwrites both.

        implicit none

        type(bessi_class),       intent(INOUT) :: bsi
        type(chion_const_class), intent(IN)    :: c

        bsi%now%n_lay       = 0
        bsi%now%mass        = 0.0_wp
        bsi%now%mass_w      = 0.0_wp
        bsi%now%density     = bsi%par%density_init
        bsi%now%temperature = bsi%par%temperature_init

        bsi%now%mass_base            = 0.0_wp_acc
        bsi%now%smb_ice              = 0.0_wp_acc
        bsi%now%runoff               = 0.0_wp_acc
        bsi%now%melt                 = 0.0_wp_acc
        bsi%now%refreezing           = 0.0_wp_acc
        bsi%now%vapor_mass           = 0.0_wp_acc
        bsi%now%sublimation          = 0.0_wp_acc
        bsi%now%latent_heat_flux_sum = 0.0_wp_acc

        bsi%now%t_srf  = c%T0
        bsi%now%albedo = c%alpha_dry
        bsi%now%w_snow_old = 0.0_wp
        bsi%now%w_snow_max = 0.0_wp

        bsi%now%thickness    = 0.0_wp
        bsi%now%wet_mass     = 0.0_wp
        bsi%now%bulk_density = 0.0_wp
        bsi%now%liquid_water = 0.0_wp

        return

    end subroutine bessi_init_state

    subroutine bessi_reset_columns(bsi,c,idx)
        ! Chion.jl _reset_bessi_columns_kernel! (src/runtime.jl:98-142), called
        ! when columns are switched off by set_active_mask!.
        !
        ! Identical to bessi_init_state per column EXCEPT that the four
        ! diagnostics (thickness, wet_mass, bulk_density, liquid_water) are
        ! deliberately NOT reset -- Julia's reset kernel does not take them.
        ! They are pure diagnostics recomputed by summarize_domain_state, so a
        ! deactivated column keeps its last values until something recomputes
        ! them. Preserved rather than "fixed".

        implicit none

        type(bessi_class),       intent(INOUT) :: bsi
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(IN)    :: idx(:)   ! columns to reset

        ! Local variables
        integer :: i, icol

        do i = 1, size(idx)

            icol = idx(i)

            if (icol .lt. 1 .or. icol .gt. bsi%now%ncol) then
                write(io_unit_err,*) "bessi_reset_columns:: Error: column index out of range."
                write(io_unit_err,*) "icol, ncol = ", icol, bsi%now%ncol
                stop "Program stopped."
            end if

            bsi%now%n_lay(icol)         = 0
            bsi%now%mass(:,icol)        = 0.0_wp
            bsi%now%mass_w(:,icol)      = 0.0_wp
            bsi%now%density(:,icol)     = bsi%par%density_init
            bsi%now%temperature(:,icol) = bsi%par%temperature_init

            bsi%now%mass_base(icol)            = 0.0_wp_acc
            bsi%now%smb_ice(icol)              = 0.0_wp_acc
            bsi%now%runoff(icol)               = 0.0_wp_acc
            bsi%now%melt(icol)                 = 0.0_wp_acc
            bsi%now%refreezing(icol)           = 0.0_wp_acc
            bsi%now%vapor_mass(icol)           = 0.0_wp_acc
            bsi%now%sublimation(icol)          = 0.0_wp_acc
            bsi%now%latent_heat_flux_sum(icol) = 0.0_wp_acc

            bsi%now%t_srf(icol)  = c%T0
            bsi%now%albedo(icol) = c%alpha_dry
            bsi%now%w_snow_old(icol) = 0.0_wp
            bsi%now%w_snow_max(icol) = 0.0_wp

        end do

        return

    end subroutine bessi_reset_columns

    ! =====================================================================
    ! The column kernel
    ! =====================================================================

    subroutine bessi_column_step_core(bsi,icol,forc,c)
        ! Chion.jl column_step_core! (src/step.jl:159-407).
        !
        ! Advances ONE column by ONE (sub)step. The grid loop belongs to the
        ! dispatcher (WP11); this routine knows nothing about neighbours.
        !
        ! Deviation from Julia, structural only (docs/porting_notes.md D8):
        ! the kernels take contiguous column slices mass(:,icol) plus the
        ! active layer count, instead of the full (Ntot,ncol) array plus an
        ! index. Fortran is column-major, so a slice is contiguous and is
        ! passed by reference -- no copy, and OpenMP privacy is obvious by
        ! inspection.

        implicit none

        type(bessi_class),              intent(INOUT) :: bsi
        integer,                        intent(IN)    :: icol
        type(chion_step_forcing_class), intent(IN)    :: forc
        type(chion_const_class),        intent(IN)    :: c

        ! Local variables
        real(wp) :: dt_seconds, accumulation_rate, melt_mass
        logical  :: started_without_surface_snow
        logical  :: use_prescribed_albedo
        logical  :: has_surface_snow, has_liquid_water
        logical  :: uses_htessel
        integer  :: n_liquid_water_before_energy
        real(wp_acc) :: melted, ice_melt, routed_runoff, refrozen_mass

        type(bare_ice_ablation_class)  :: bare_ice_fluxes
        type(latent_heat_coeff_class)  :: lh_coef
        type(snow_energy_result_class) :: energy
        type(surface_vapor_flux_class) :: snow_vapor_fluxes

        ! HTESSEL snapshot of the liquid water held BEFORE the energy solve.
        ! Automatic array, stack-local, OpenMP-private by construction. This
        ! replaces Chion.jl's persistent workspace.liquid_water_before_energy.
        real(wp) :: mass_w_before_energy(bsi%par%Ntot)

        ! SEMIX albedo scratch: column SWE, dust concentration, and the
        ! forcing scalars the albedo module takes (resolved here, since the
        ! module stays decoupled from the forcing type).
        real(wp) :: w_snow, dust_con, coszm, cloud, z_sur_std

        ! Bare-ice albedo actually used: chion's own constant unless the host
        ! supplies one per column (e.g. CLIMBER-X's slow firn-aging ice albedo).
        real(wp) :: alb_ice_use

        if (icol .lt. 1 .or. icol .gt. bsi%now%ncol) then
            write(io_unit_err,*) "bessi_column_step_core:: Error: column index out of range."
            write(io_unit_err,*) "icol, ncol = ", icol, bsi%now%ncol
            stop "Program stopped."
        end if

        associate(n           => bsi%now%n_lay(icol),               &
                  mass        => bsi%now%mass(:,icol),              &
                  mass_w      => bsi%now%mass_w(:,icol),            &
                  density     => bsi%now%density(:,icol),           &
                  temperature => bsi%now%temperature(:,icol),       &
                  mass_base   => bsi%now%mass_base(icol),           &
                  smb_ice     => bsi%now%smb_ice(icol),             &
                  runoff      => bsi%now%runoff(icol),              &
                  melt        => bsi%now%melt(icol),                &
                  refreezing  => bsi%now%refreezing(icol),          &
                  vapor_mass  => bsi%now%vapor_mass(icol),          &
                  sublimation => bsi%now%sublimation(icol),         &
                  lhf_sum     => bsi%now%latent_heat_flux_sum(icol),&
                  t_srf       => bsi%now%t_srf(icol),               &
                  albedo      => bsi%now%albedo(icol),              &
                  w_snow_old  => bsi%now%w_snow_old(icol),          &
                  w_snow_max  => bsi%now%w_snow_max(icol),          &
                  par         => bsi%par)

        ! === Step 1: setup ===================================================
        ! started_without_surface_snow MUST be sampled before accumulation:
        ! it is what tells step 3 that the snow now sitting on layer 1 fell
        ! this step onto bare ground, and so should carry the air temperature
        ! rather than whatever temperature_init the slot was reset to.

        dt_seconds = forc%dt_days*c%seconds_per_day

        started_without_surface_snow = .not. surface_has_snow(mass,n)

        ! Trap 9: PRESCRIBED is the only scheme the caller overrides, and only
        ! when the forcing actually supplies an albedo. Without one it silently
        ! behaves as dynamic. Nested if -- .and. does not short-circuit.
        use_prescribed_albedo = .FALSE.
        if (c%albedo_scheme .eq. CHION_ALBEDO_PRESCRIBED) then
            if (forc%has_prescribed_albedo) use_prescribed_albedo = .TRUE.
        end if

        alb_ice_use = c%alpha_ice
        if (forc%has_alb_ice_host) &
            alb_ice_use = min(max(forc%alb_ice_host,0.0_wp),1.0_wp)

        uses_htessel = (c%low_density_densification .eq. CHION_DENSIFY_HTESSEL)

        ! === Step 2: accumulation ============================================

        call apply_accumulation(mass,mass_w,density,temperature,n, &
                                mass_base,smb_ice,runoff,t_srf,albedo, &
                                c,par%Ntot,par%mass_max,par%mass_split,par%mass_min, &
                                forc%snowfall_rate,forc%rainfall_rate,dt_seconds, &
                                forc%air_temperature,forc%wind_speed)

        ! === Step 3: fresh snow onto a column that started bare ==============
        ! Nested ifs: n > 0 must be established before layer 1 is meaningful.

        if (forc%snowfall_rate .gt. 0.0_wp) then
            if (started_without_surface_snow) then
                if (n .gt. 0) temperature(1) = forc%air_temperature
            end if
        end if

        ! === Step 4: prescribed albedo, applied before the bare test =========

        if (use_prescribed_albedo) then
            albedo = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
        end if

        ! === Step 5: bare surface -> ablate the ice underneath and RETURN ====
        ! This is the single most consequential branch in the model. Note it
        ! is decided AFTER accumulation, so a snowfall large enough to clear
        ! TOL_EMPTY_LAYER rescues the column within the same step.
        !
        ! The return is unconditional: percolation, refreezing and the HTESSEL
        ! compaction never run on a bare column, so liquid water already in
        ! mass_w stays exactly where it is.

        has_surface_snow = surface_has_snow(mass,n)

        if (.not. has_surface_snow) then

            if (.not. use_prescribed_albedo) albedo = alb_ice_use

            bare_ice_fluxes = bare_ice_ablation_mass(c,forc,dt_seconds)

            ! Accumulate directly into the wp_acc accumulators.
            smb_ice     = smb_ice     + real(bare_ice_fluxes%net_mass_change,wp_acc)
            melt        = melt        + real(bare_ice_fluxes%melt_mass,wp_acc)
            runoff      = runoff      + real(bare_ice_fluxes%melt_mass,wp_acc)
            vapor_mass  = vapor_mass  + real(bare_ice_fluxes%vapor_mass,wp_acc)
            sublimation = sublimation + real(bare_ice_fluxes%sublimation_mass,wp_acc)

            ! W m-2 DAYS, not W m-2 seconds. dt_days, not dt_seconds.
            lhf_sum = lhf_sum + real(bare_ice_fluxes%latent_heat_flux,wp_acc) &
                                *real(forc%dt_days,wp_acc)

            return

        end if

        ! === Step 6: albedo update ===========================================
        ! Trap 5: the aging law carries no dt, so albedo_update must be called
        ! exactly once per step. That is here.

        if (use_prescribed_albedo) then
            albedo = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
        else if (c%albedo_scheme .eq. CHION_ALBEDO_SEMIX) then
            ! Column SWE and its seasonal peak drive the dust melt-amplification:
            ! meltwater scavenges little dust, so what remains concentrates as
            ! the pack draws down from its peak. Tracked here, where the SEMIX
            ! albedo is the only consumer, so the other schemes pay nothing.
            w_snow = sum(mass(1:n)) + sum(mass_w(1:n))
            if (w_snow .gt. w_snow_old) w_snow_max = w_snow

            dust_con = 0.0_wp
            if (forc%has_dust_dep) &
                dust_con = semix_dust_concentration(forc%dust_dep,forc%snowfall_rate, &
                                                    w_snow,w_snow_max,c)

            if (forc%has_coszm) then
                coszm = forc%coszm
            else
                coszm = semix_daily_coszm(forc%latitude_deg,forc%solar_longitude_deg)
            end if

            cloud = 0.0_wp
            if (forc%has_cloud) cloud = forc%cloud

            z_sur_std = 0.0_wp
            if (forc%has_z_sur_std) z_sur_std = forc%z_sur_std

            call semix_surface_albedo(c,temperature(1),forc%snowfall_rate, &
                                      coszm,cloud,z_sur_std,dust_con,albedo)

            w_snow_old = w_snow
        else
            call albedo_update(mass,mass_w,density,temperature,n,c,albedo)
        end if

        ! === Step 7: HTESSEL snapshot ========================================
        ! Taken BEFORE densification and the energy solve, consumed in step 13.
        ! n_liquid_water_before_energy doubles as the "snapshot was taken" flag,
        ! exactly as in Julia.
        !
        ! The whole array is zeroed first, unlike Julia's persistent workspace
        ! which retains values from the previous column. That difference is
        ! unreachable: between here and step 13 the layer count can only fall
        ! (melt and sublimation remove layers; splits happen in accumulation,
        ! which is already done), so step 13 never reads a slot this loop did
        ! not write.

        n_liquid_water_before_energy = 0

        if (uses_htessel) then
            mass_w_before_energy = 0.0_wp
            n_liquid_water_before_energy = n
            mass_w_before_energy(1:n) = mass_w(1:n)
        end if

        ! === Step 8: densification ===========================================
        ! Rain counts toward the accumulation rate only when the surface has
        ! snow -- which it does here, since the bare branch already returned.
        ! The conditional is kept because it is in the Julia source.

        accumulation_rate = max(forc%snowfall_rate,0.0_wp)
        if (has_surface_snow) accumulation_rate = accumulation_rate + forc%rainfall_rate

        call densify_column(mass,density,temperature,n,c,accumulation_rate,dt_seconds)

        ! === Step 9: latent-heat coefficients, then the energy solve =========

        lh_coef = diagnose_latent_heat_flux_coefficients(has_surface_snow,c, &
                                                         forc%air_temperature, &
                                                         forc%snowfall_rate, &
                                                         forc%rainfall_rate)

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,albedo, &
                              lh_coef%linear,lh_coef%constant,dt_seconds,energy)

        ! === Step 10: post-solve surface vapor mass flux =====================
        ! Evaluated exactly at the NEW surface temperature, deliberately
        ! inconsistent with the linearization used inside the solve (trap 2).

        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n, &
                                                runoff,t_srf,albedo,c,forc,dt_seconds, &
                                                par%mass_split,par%mass_min, &
                                                snow_vapor_fluxes)

        vapor_mass  = vapor_mass  + real(snow_vapor_fluxes%vapor_mass,wp_acc)
        sublimation = sublimation + real(snow_vapor_fluxes%sublimation_mass,wp_acc)
        lhf_sum     = lhf_sum     + real(snow_vapor_fluxes%latent_heat_flux,wp_acc) &
                                    *real(forc%dt_days,wp_acc)

        ! === Step 11: melt ===================================================
        ! melt_energy_available is wp_acc; melt_mass is the wp mass handed to
        ! apply_melt, matching the kernel's interface.
        !
        ! TRAP: the melt diagnostic receives the FULL REQUESTED melt_mass, not
        ! `melted`. When the snowpack cannot supply the demand AND the column
        ! is now empty, the shortfall is charged to the ice: smb_ice loses it
        ! and runoff gains it. Both halves of that guard matter -- a shortfall
        ! with layers still present (which apply_melt's general path can
        ! produce) is silently dropped.

        if (energy%needs_melt) then

            melt_mass = real(energy%melt_energy_available/real(c%Lm,wp_acc),wp)

            call apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo, &
                            par%mass_split,par%mass_min,melt_mass,c,melted)

            if (melted .lt. real(melt_mass,wp_acc)) then
                if (n .eq. 0) then
                    ice_melt = real(melt_mass,wp_acc) - melted
                    smb_ice  = smb_ice - ice_melt
                    runoff   = runoff  + ice_melt
                end if
            end if

            melt = melt + real(melt_mass,wp_acc)

        end if

        ! === Step 12: percolation ============================================

        has_liquid_water = column_has_liquid_water(mass_w,n)

        if (has_liquid_water) then

            call apply_percolation(mass,mass_w,density,n,c%rho_i,c%rho_w,routed_runoff)

            runoff = runoff + routed_runoff

            has_liquid_water = column_has_liquid_water(mass_w,n)

        end if

        ! === Step 13: HTESSEL liquid-water compaction ========================
        ! Three-way guard, all three parts required. Written nested because
        ! Fortran evaluates .and. operands in unspecified order.

        if (uses_htessel) then
            if (n_liquid_water_before_energy .gt. 0) then
                if (has_liquid_water) then
                    call apply_htessel_liquid_water_compaction(mass,mass_w,density,n, &
                                                               mass_w_before_energy,c)
                end if
            end if
        end if

        ! === Step 14: refreezing =============================================

        if (has_liquid_water) then

            call apply_refreezing(mass,mass_w,density,temperature,n, &
                                  c%T0,c%ci,c%Lm,c%rho_i,refrozen_mass)

            refreezing = refreezing + refrozen_mass

            ! Julia recomputes has_liquid_water here and never reads it again
            ! (step.jl:383). Dropped as dead code.

        end if

        ! === Step 15: final albedo fixup =====================================
        ! The column may have gone bare during melt or sublimation, in which
        ! case the albedo diagnosed in step 6 no longer describes the surface.

        if (use_prescribed_albedo) then
            albedo = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
        else if (.not. surface_has_snow(mass,n)) then
            albedo = alb_ice_use
        end if

        end associate

        return

    end subroutine bessi_column_step_core

    subroutine bessi_column_step(bsi,icol,forc,c)
        ! Chion.jl column_step! (src/step.jl:113-157): the diurnal-shortwave
        ! substep wrapper around the core kernel.
        !
        ! When enabled and warranted, the day is tiled uniformly in hour angle
        ! over [-pi, pi] and the core is run once per interval with:
        !   dt_days         scaled by the interval's fraction of the day,
        !   shortwave_down  the interval mean (nocturnal intervals get ~0),
        !   q_sw_net        likewise, if prescribed,
        !   air_temperature the interval mean, if the temperature cycle is on.
        !
        ! Everything else, INCLUDING PRECIPITATION, is passed through
        ! unchanged. Precipitation is a RATE, so shrinking dt_days shrinks the
        ! mass added per substep and the daily total is conserved exactly. The
        ! acceptance test asserts this.
        !
        ! Melt, by contrast, is NOT conserved and is not meant to be: melting
        ! is a rectified function of the surface energy balance, so resolving
        ! the daytime peak produces more melt than the daily mean does. That is
        ! the entire purpose of the option.

        implicit none

        type(bessi_class),              intent(INOUT) :: bsi
        integer,                        intent(IN)    :: icol
        type(chion_step_forcing_class), intent(IN)    :: forc
        type(chion_const_class),        intent(IN)    :: c

        ! Local variables
        integer  :: n_substeps, k
        real(wp) :: shortwave_for_criterion
        real(wp) :: hour_angle_start, hour_angle_end, fraction

        type(chion_step_forcing_class) :: subforc

        if (bsi%par%diurnal_shortwave_substeps) then

            ! The criterion uses the net shortwave when it is prescribed,
            ! because that is what actually drives the surface.
            if (forc%has_q_sw_net) then
                shortwave_for_criterion = forc%q_sw_net
            else
                shortwave_for_criterion = forc%shortwave_down
            end if

            n_substeps = diurnal_substep_count(forc%dt_days,shortwave_for_criterion, &
                                               forc%air_temperature, &
                                               bsi%par%diurnal_shortwave_min_air_temperature, &
                                               forc%latitude_deg,forc%solar_longitude_deg, &
                                               bsi%par%diurnal_shortwave_threshold, &
                                               bsi%par%diurnal_shortwave_max_substeps)

            if (n_substeps .gt. 1) then

                do k = 1, n_substeps

                    call diurnal_substep_bounds(k,n_substeps,hour_angle_start,hour_angle_end)

                    fraction = real((real(hour_angle_end,wp_acc) - real(hour_angle_start,wp_acc)) &
                                    /(2.0_wp_acc*PI_ACC),wp)

                    ! Julia skips a non-positive interval rather than erroring.
                    if (fraction .le. 0.0_wp) cycle

                    subforc = forc

                    subforc%dt_days = forc%dt_days*fraction

                    subforc%shortwave_down = &
                        diurnal_shortwave_interval_average(forc%shortwave_down, &
                                                           forc%latitude_deg, &
                                                           forc%solar_longitude_deg, &
                                                           hour_angle_start,hour_angle_end)

                    if (forc%has_q_sw_net) then
                        subforc%q_sw_net = &
                            diurnal_shortwave_interval_average(forc%q_sw_net, &
                                                               forc%latitude_deg, &
                                                               forc%solar_longitude_deg, &
                                                               hour_angle_start,hour_angle_end)
                    end if

                    if (bsi%par%diurnal_temperature_cycle) then
                        subforc%air_temperature = &
                            diurnal_temperature_interval_average(forc%air_temperature, &
                                                bsi%par%diurnal_temperature_amplitude, &
                                                hour_angle_start,hour_angle_end)
                    end if

                    call bessi_column_step_core(bsi,icol,subforc,c)

                end do

                return

            end if

        end if

        call bessi_column_step_core(bsi,icol,forc,c)

        return

    end subroutine bessi_column_step

end module snow_bessi
