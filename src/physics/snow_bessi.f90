module snow_bessi
    ! WP8 -- BESSI model assembly.
    !
    ! Fortran port of Chion.jl (branch main):
    !   src/models.jl:11-78      BESSIModel        -> bessi_par_class
    !   src/state.jl:11-134      BESSIState        -> bessi_state_class
    !   src/runtime.jl:98-175    column reset      -> bessi_reset_columns
    !   src/step.jl:159-407      column_step_core! -> bessi_column_step_core
    !   src/step.jl:75-157       column_step!      -> bessi_column_step
    !   src/step.jl _step_bare_ice_substrate! (9ec6cc7) -> bessi_bare_ice_substrate_step
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
    !   3  fresh snow onto a column that started bare -> temperature(1:n) = T_air,
    !      then the near-surface remesh (fine layers only)
    !   4  prescribed albedo, if any, applied BEFORE the bare-surface test
    !   5  no surface snow -> bare-ice ablation, accumulate diagnostics, RETURN
    !      (percolation and refreezing are skipped ENTIRELY -- trap, see below);
    !      with an ice substrate, the substrate energy solve instead of the
    !      surface held at T0 (bessi_bare_ice_substrate_step); on a land column
    !      (H_ice = 0) nothing at all (D41)
    !   6  snow albedo update (prescribed | semix | aging | dynamic-or-constant),
    !      then the albedo the energy balance sees: the snow albedo blended
    !      with the background by the snow-cover fraction (D40)
    !   7  HTESSEL only: snapshot mass_w(1:n) BEFORE the energy solve
    !   8  accumulation_rate -> densification
    !   9  latent-heat coefficients -> implicit energy solve
    !  10  post-solve surface vapor mass flux; accumulate its three outputs
    !  11  melt, if the solve asked for it; ice melt if the snow ran out
    !  12  percolation -> runoff
    !  13  HTESSEL only: liquid-water compaction, using the step-7 snapshot
    !  14  refreezing
    !  15  near-surface remesh (fine layers only)
    !  16  final albedo fixup (and snow age), the blend on the final column
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

    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    use chion_defs, only : wp, wp_acc, io_unit_err, &
                           TOL_TINY, TOL_EMPTY_LAYER, &
                           ALBEDO_THIN_SNOW_BLEND, LAND_COLUMNS_WITHOUT_ICE, &
                           DEF_NTOT, DEF_MASS_MAX, DEF_MASS_SPLIT, DEF_MASS_MIN, &
                           DEF_DENSITY_INIT, DEF_TEMPERATURE_INIT, &
                           DEF_ICE_SUBSTRATE_LAYERS, DEF_ICE_SUBSTRATE_TOP_THICKNESS, &
                           NEAR_SURFACE_LAYERS, DEF_NEAR_SURFACE_LAYER_MAX_THICKNESSES, &
                           CHION_ALBEDO_PRESCRIBED, CHION_ALBEDO_SEMIX, &
                           CHION_ALBEDO_AGING, &
                           CHION_DENSIFY_HTESSEL, DIURNAL_SINGLE_INTERVAL_AVERAGED, &
                           chion_const_class, chion_step_forcing_class

    use phys_constants,     only : sec_day
    use snow_column_utils,  only : surface_has_snow, column_has_liquid_water

    use snow_accumulation,  only : apply_accumulation
    use snow_layers,        only : remesh_near_surface_layers
    use snow_albedo,        only : albedo_update, albedo_update_aging, &
                                   snow_cover_fraction, thin_snow_albedo
    use snow_albedo_semix,  only : semix_surface_albedo, semix_dust_concentration, &
                                   semix_daily_coszm, semix_snow_cover_fraction
    use snow_densify,       only : densify_column, apply_htessel_liquid_water_compaction
    use snow_energy,        only : snow_energy_flux, snow_energy_result_class
    use snow_surface_fluxes,only : bare_ice_ablation_class, bare_ice_ablation_mass, &
                                   resolved_turbulent_latent_heat_flux, &
                                   latent_heat_coeff_class, &
                                   diagnose_latent_heat_flux_coefficients, &
                                   surface_vapor_flux_class, &
                                   apply_snow_surface_vapor_mass_flux, &
                                   with_parameterized_longwave
    use snow_melt,          only : apply_melt
    use snow_percolation,   only : apply_percolation
    use snow_refreezing,    only : apply_refreezing
    use snow_diurnal,       only : diurnal_substep_count, diurnal_substep_bounds, &
                                   diurnal_shortwave_interval_average, &
                                   diurnal_temperature_amplitude, &
                                   diurnal_temperature_interval_average

    implicit none

    private

    real(wp_acc), parameter :: PI_ACC = 3.141592653589793238462643383279502884_wp_acc

    ! === Types ===============================================================

    type bessi_par_class
        ! Chion.jl BESSIModel (src/models.jl:11-26), minus the grid and the
        ! constants object, which chion keeps separately (chion_grid_class and
        ! chion_const_class). The nine diurnal fields are Chion.jl's step
        ! `config` NamedTuple (src/step.jl:45-52), which is built directly from
        ! the model (src/simulation.jl:172-180), so they belong here.

        integer  :: Ntot                  ! [1] layer capacity
        real(wp) :: mass_max              ! [kg m-2] surface-layer split trigger
        real(wp) :: mass_split            ! [kg m-2] mass left below after a split
        real(wp) :: mass_min              ! [kg m-2] surface-layer merge trigger
        real(wp) :: density_init          ! [kg m-3] density of an inactive layer
        real(wp) :: temperature_init      ! [K] temperature of an inactive layer

        ! Thermal ice substrate (Chion.jl 03bb445): layers below the snow/firn
        ! column, top thickness doubling downward; 0 = none (bare ice held at
        ! T0). Only under ice: a column with H_ice = 0 has none (D34).
        integer  :: ice_substrate_layers          ! [1] n_ice >= 0
        real(wp) :: ice_substrate_top_thickness   ! [m] top layer, > 0

        ! Fine near-surface layers (Chion.jl 03bb445): maximum thicknesses
        ! of the top layers, held by a conservative remesh; 0 = no limit
        ! (Julia Inf), the limited ones a leading block. With layer 1
        ! limited, the 100 kg m-2 surface merge (mass_min) is off.
        real(wp) :: near_surface_layer_max_thicknesses(NEAR_SURFACE_LAYERS)   ! [m] >= 0

        logical  :: diurnal_shortwave_substeps            ! [1] enable substepping
        real(wp) :: diurnal_shortwave_threshold           ! [W m-2] peak-minus-mean excess
        integer  :: diurnal_shortwave_max_substeps        ! [1] 1..24
        real(wp) :: diurnal_shortwave_min_air_temperature ! [K]
        logical  :: diurnal_temperature_cycle             ! [1] impose a diurnal T cycle
        real(wp) :: diurnal_temperature_amplitude         ! [K] half-amplitude at z_ref
        real(wp) :: diurnal_temperature_amplitude_gradient          ! [K km-1] above z_ref
        real(wp) :: diurnal_temperature_amplitude_reference_height  ! [m] z_ref
        real(wp) :: diurnal_temperature_amplitude_max               ! [K] upper clamp
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

        ! Thermal ice substrate, top layer first. Extent 0 without one. A land
        ! column (H_ice = 0) never touches its slice.
        integer               :: n_ice
        real(wp), allocatable :: ice_temperature(:,:)  ! (n_ice,ncol) [K]

        ! Cumulative accumulators -- MUST be wp_acc. docs/PLAN.md section 3.1.
        real(wp_acc), allocatable :: mass_base(:)             ! [kg m-2]
        real(wp_acc), allocatable :: smb_ice(:)               ! [kg m-2]
        real(wp_acc), allocatable :: runoff(:)                ! [kg m-2]
        real(wp_acc), allocatable :: melt(:)                  ! [kg m-2]
        real(wp_acc), allocatable :: refreezing(:)            ! [kg m-2]
        real(wp_acc), allocatable :: vapor_mass(:)            ! [kg m-2] + = deposition
        real(wp_acc), allocatable :: sublimation(:)           ! [kg m-2] >= 0
        real(wp_acc), allocatable :: latent_heat_flux_sum(:)  ! [W m-2 d]

        ! Instantaneous per-column scalars. albedo is what the surface energy
        ! balance sees (and what is written out): the snow albedo blended with
        ! the background by the snow-cover fraction (D40). albedo_snow is the
        ! snow's own, the prognostic one the schemes age and refresh; on a
        ! bare column it holds the background, as Chion.jl's single albedo
        ! does. Without the blend (legacy_chion) the two agree under snow.
        real(wp), allocatable :: t_srf(:)          ! (ncol) [K] snow-air interface temperature
        real(wp), allocatable :: albedo(:)         ! (ncol) [1] effective
        real(wp), allocatable :: albedo_snow(:)    ! (ncol) [1] snow

        ! Time since the latest snowfall (albedo_scheme = "aging"); 0 on a
        ! bare column and under every other scheme. Chion.jl snow_age_days.
        real(wp), allocatable :: snow_age_days(:)  ! (ncol) [d]

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
        ! Chion.jl defaults: src/models.jl:74-100 at 9ec6cc7 (BESSIModel
        ! keyword defaults, the calibrated GrIS set: 8 diurnal substeps, a
        ! 1 K temperature cycle capped at 1 K) and src/constants.jl
        ! (DEFAULT_NTOT etc., re-exported by chion_defs).
        !
        ! NOTE the two unit conversions Julia performs in the constructor:
        !   diurnal_shortwave_min_air_temperature_c = -8.0  ->  265.15 K
        !   diurnal_temperature_amplitude_c         =  1.0  ->    1.0 K
        !   diurnal_temperature_amplitude_gradient_c_per_km -> K km-1 here; the
        !     /1000 to K m-1 happens in diurnal_temperature_amplitude
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

        par%ice_substrate_layers        = DEF_ICE_SUBSTRATE_LAYERS
        par%ice_substrate_top_thickness = DEF_ICE_SUBSTRATE_TOP_THICKNESS

        par%near_surface_layer_max_thicknesses = DEF_NEAR_SURFACE_LAYER_MAX_THICKNESSES

        par%diurnal_shortwave_substeps            = .TRUE.
        par%diurnal_shortwave_threshold           = 0.0_wp
        par%diurnal_shortwave_max_substeps        = 8
        par%diurnal_shortwave_min_air_temperature = 265.15_wp
        par%diurnal_temperature_cycle             = .TRUE.
        par%diurnal_temperature_amplitude         = 1.0_wp
        par%diurnal_temperature_amplitude_gradient         = 0.0_wp
        par%diurnal_temperature_amplitude_reference_height = 0.0_wp
        par%diurnal_temperature_amplitude_max              = 1.0_wp

        call bessi_par_validate(par)

        return

    end subroutine bessi_par_init

    subroutine bessi_par_validate(par)
        ! Chion.jl _validate_mass_partition (src/domain.jl:59-64) plus the six
        ! diurnal guards in the BESSIModel constructor (src/models.jl:48-51, 97-101 at 27113b6)
        ! and the two ice-substrate guards (src/models.jl:114-115 at 9ec6cc7)
        ! and the near-surface thickness guard (src/models.jl:112-113).
        !
        ! The mass partition is not merely cosmetic:
        !   mass_split < mass_max   -- otherwise the split loop cannot converge
        !   mass_min   < mass_split -- otherwise a split immediately re-merges
        !   mass_split/mass_max >= 0.5 -- so the layer left behind by a split
        !                                 holds at least half of mass_max, and
        !                                 merge_layer's unguarded
        !                                 divisor (2*mass_split - mass_min)
        !                                 stays positive. See upstream defect 8
        !                                 in docs/porting_notes.md.

        implicit none

        type(bessi_par_class), intent(IN) :: par

        ! Local variables
        integer :: k

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

        if (par%ice_substrate_layers .lt. 0) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &ice_substrate_layers must be non-negative."
            write(io_unit_err,*) "ice_substrate_layers = ", par%ice_substrate_layers
            stop "Program stopped."
        end if

        if (.not. (par%ice_substrate_top_thickness .gt. 0.0_wp)) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &ice_substrate_top_thickness must be positive."
            write(io_unit_err,*) "ice_substrate_top_thickness = ", par%ice_substrate_top_thickness
            stop "Program stopped."
        end if

        ! Julia requires positive values, Inf meaning no limit; chion writes
        ! no limit as 0. chion also requires the limited layers to be the
        ! top ones (D36): the first unlimited layer is then the first mass
        ! layer, which the split/merge below the fine layers acts on.
        do k = 1, NEAR_SURFACE_LAYERS
            if (.not. ieee_is_finite(par%near_surface_layer_max_thicknesses(k)) .or. &
                par%near_surface_layer_max_thicknesses(k) .lt. 0.0_wp) then
                write(io_unit_err,*) "bessi_par_validate:: Error: &
                                     &near_surface_layer_max_thicknesses must be finite and &
                                     &non-negative (0 = no limit)."
                write(io_unit_err,*) "near_surface_layer_max_thicknesses = ", &
                                     par%near_surface_layer_max_thicknesses
                stop "Program stopped."
            end if
        end do

        do k = 2, NEAR_SURFACE_LAYERS
            if (par%near_surface_layer_max_thicknesses(k)   .gt. 0.0_wp .and. &
                par%near_surface_layer_max_thicknesses(k-1) .le. 0.0_wp) then
                write(io_unit_err,*) "bessi_par_validate:: Error: &
                                     &near_surface_layer_max_thicknesses: only the top layers &
                                     &can be limited (no limit, 0, above a limited layer)."
                write(io_unit_err,*) "near_surface_layer_max_thicknesses = ", &
                                     par%near_surface_layer_max_thicknesses
                stop "Program stopped."
            end if
        end do

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

        if (.not. ieee_is_finite(par%diurnal_temperature_amplitude_gradient)) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &diurnal_temperature_amplitude_gradient must be finite."
            write(io_unit_err,*) "diurnal_temperature_amplitude_gradient = ", &
                                 par%diurnal_temperature_amplitude_gradient
            stop "Program stopped."
        end if

        if (.not. (par%diurnal_temperature_amplitude_max .ge. &
                   par%diurnal_temperature_amplitude)) then
            write(io_unit_err,*) "bessi_par_validate:: Error: &
                                 &diurnal_temperature_amplitude_max must be at least &
                                 &diurnal_temperature_amplitude."
            write(io_unit_err,*) "diurnal_temperature_amplitude, _max = ", &
                                 par%diurnal_temperature_amplitude, &
                                 par%diurnal_temperature_amplitude_max
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

        bsi%now%ncol  = ncol
        bsi%now%Ntot  = bsi%par%Ntot
        bsi%now%n_ice = bsi%par%ice_substrate_layers

        allocate(bsi%now%n_lay(ncol))
        allocate(bsi%now%mass(bsi%par%Ntot,ncol))
        allocate(bsi%now%mass_w(bsi%par%Ntot,ncol))
        allocate(bsi%now%density(bsi%par%Ntot,ncol))
        allocate(bsi%now%temperature(bsi%par%Ntot,ncol))
        allocate(bsi%now%ice_temperature(bsi%now%n_ice,ncol))

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
        allocate(bsi%now%albedo_snow(ncol))
        allocate(bsi%now%snow_age_days(ncol))
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
        if (allocated(bsi%now%ice_temperature)) deallocate(bsi%now%ice_temperature)

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
        if (allocated(bsi%now%albedo_snow)) deallocate(bsi%now%albedo_snow)
        if (allocated(bsi%now%snow_age_days)) deallocate(bsi%now%snow_age_days)
        if (allocated(bsi%now%w_snow_old)) deallocate(bsi%now%w_snow_old)
        if (allocated(bsi%now%w_snow_max)) deallocate(bsi%now%w_snow_max)

        if (allocated(bsi%now%thickness))    deallocate(bsi%now%thickness)
        if (allocated(bsi%now%wet_mass))     deallocate(bsi%now%wet_mass)
        if (allocated(bsi%now%bulk_density)) deallocate(bsi%now%bulk_density)
        if (allocated(bsi%now%liquid_water)) deallocate(bsi%now%liquid_water)

        bsi%now%ncol  = 0
        bsi%now%Ntot  = 0
        bsi%now%n_ice = 0

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
        ! even though the column is bare. The first step overwrites both. The
        ! ice substrate starts at temperature_init (Chion.jl state.jl).

        implicit none

        type(bessi_class),       intent(INOUT) :: bsi
        type(chion_const_class), intent(IN)    :: c

        bsi%now%n_lay       = 0
        bsi%now%mass        = 0.0_wp
        bsi%now%mass_w      = 0.0_wp
        bsi%now%density     = bsi%par%density_init
        bsi%now%temperature = bsi%par%temperature_init

        bsi%now%ice_temperature = bsi%par%temperature_init

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
        bsi%now%albedo_snow = c%alpha_dry
        bsi%now%snow_age_days = 0.0_wp
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
        !
        ! The ice substrate IS reset, to temperature_init as at a cold start
        ! (chion deviation, docs/porting_notes.md D34): Chion.jl's reset kernel
        ! leaves ice_temperature stale, so a re-activated column would start
        ! on the ice of its previous life (reported upstream, N8).

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

            bsi%now%ice_temperature(:,icol) = bsi%par%temperature_init

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
            bsi%now%albedo_snow(icol) = c%alpha_dry
            bsi%now%snow_age_days(icol) = 0.0_wp
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
        real(wp) :: surface_mass_min
        logical  :: started_without_surface_snow
        logical  :: use_prescribed_albedo
        logical  :: has_surface_snow, has_liquid_water
        logical  :: uses_htessel
        integer  :: n_liquid_water_before_energy
        integer  :: n_ice
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
        ! The background under thin or no snow: that under ice, alpha_land on
        ! a land column (D40, D41).
        real(wp) :: alb_ice_use, alb_bg
        logical  :: is_land

        if (icol .lt. 1 .or. icol .gt. bsi%now%ncol) then
            write(io_unit_err,*) "bessi_column_step_core:: Error: column index out of range."
            write(io_unit_err,*) "icol, ncol = ", icol, bsi%now%ncol
            stop "Program stopped."
        end if

        ! Ice substrate rows of this column, set before the associate below
        ! takes its slice: none on land (H_ice = 0), a chion deviation (D34);
        ! Chion.jl has no ice thickness and puts it under every column.
        n_ice = 0
        if (forc%H_ice .gt. 0.0_wp) n_ice = bsi%now%n_ice

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
                  albedo_snow => bsi%now%albedo_snow(icol),         &
                  snow_age    => bsi%now%snow_age_days(icol),       &
                  w_snow_old  => bsi%now%w_snow_old(icol),          &
                  w_snow_max  => bsi%now%w_snow_max(icol),          &
                  ice_temperature => bsi%now%ice_temperature(1:n_ice,icol), &
                  par         => bsi%par)

        ! === Step 1: setup ===================================================
        ! started_without_surface_snow MUST be sampled before accumulation:
        ! it is what tells step 3 that the snow now sitting on layer 1 fell
        ! this step onto bare ground, and so should carry the air temperature
        ! rather than whatever temperature_init the slot was reset to.

        dt_seconds = forc%dt_days*real(sec_day,wp)

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

        ! A land column (H_ice = 0) has no ice under its snow: a land
        ! background albedo, no ice ablation (D41; Chion.jl, and legacy_chion,
        ! put bare ice under every column).
        is_land = .FALSE.
        if (LAND_COLUMNS_WITHOUT_ICE) is_land = .not. (forc%H_ice .gt. 0.0_wp)

        alb_bg = alb_ice_use
        if (is_land) alb_bg = c%alpha_land

        uses_htessel = (c%low_density_densification .eq. CHION_DENSIFY_HTESSEL)

        ! The fine layers own the surface geometry: with layer 1 limited, the
        ! 100 kg m-2 surface merge in accumulation, vapour flux and melt would
        ! merge it away, so it is off there (Chion.jl 03bb445 step.jl).
        surface_mass_min = par%mass_min
        if (par%near_surface_layer_max_thicknesses(1) .gt. 0.0_wp) surface_mass_min = 0.0_wp

        ! === Step 2: accumulation ============================================

        call apply_accumulation(mass,mass_w,density,temperature,n, &
                                mass_base,smb_ice,runoff,t_srf,albedo_snow,snow_age, &
                                c,par%Ntot,par%mass_max,par%mass_split,surface_mass_min, &
                                forc%snowfall_rate,forc%rainfall_rate,dt_seconds, &
                                forc%air_temperature,forc%wind_speed)

        ! === Step 3: fresh snow onto a column that started bare ==============
        ! Every active layer takes the air temperature (Chion.jl 03bb445), not
        ! only layer 1: a snowfall above mass_max splits in accumulation and
        ! copies the reset slot temperature into layer 2, and a column with
        ! 0 < mass(1) <= TOL_EMPTY_LAYER already had layers below.

        if (forc%snowfall_rate .gt. 0.0_wp) then
            if (started_without_surface_snow) then
                temperature(1:n) = forc%air_temperature
            end if
        end if

        ! Fine layers: cap the fresh snow down into the column (after the air
        ! temperature is set, so every layer it reaches carries it).
        call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                        mass_base,smb_ice,runoff,t_srf,albedo_snow, &
                                        par%Ntot,par%mass_max,par%mass_split,par%mass_min, &
                                        par%near_surface_layer_max_thicknesses,c)

        ! === Step 4: prescribed albedo, applied before the bare test =========
        ! It replaces the snow albedo too, as it does Chion.jl's one albedo,
        ! and is never blended (D40).

        if (use_prescribed_albedo) then
            albedo_snow = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
            albedo      = albedo_snow
        end if

        ! === Step 5: bare surface -> ablate the ice underneath and RETURN ====
        ! This is the single most consequential branch in the model. Note it
        ! is decided AFTER accumulation, so a snowfall large enough to clear
        ! TOL_EMPTY_LAYER rescues the column within the same step.
        !
        ! The return is unconditional: percolation, refreezing and the HTESSEL
        ! compaction never run on a bare column, so liquid water already in
        ! mass_w stays exactly where it is.
        !
        ! With an ice substrate the bare surface is solved for, not held at T0
        ! (Chion.jl 03bb445 _step_bare_ice_substrate!).

        has_surface_snow = surface_has_snow(mass,n)

        if (.not. has_surface_snow) then

            if (.not. use_prescribed_albedo) then
                albedo_snow = alb_bg
                albedo      = alb_bg
            end if
            snow_age = 0.0_wp

            ! Snow-free land (D41): BESSI has no ground model, so the column
            ! exchanges neither mass nor energy with the atmosphere; its
            ! surface temperature is the air temperature, which is also what
            ! the next snowfall starts at (step 3). Rain has already run off
            ! (D29).
            if (is_land) then
                t_srf = forc%air_temperature
                return
            end if

            if (n_ice .gt. 0) then
                call bessi_bare_ice_substrate_step(mass,density,temperature,n, &
                                                   ice_temperature,t_srf,albedo, &
                                                   smb_ice,runoff,melt,vapor_mass, &
                                                   sublimation,lhf_sum,c,forc, &
                                                   par%ice_substrate_top_thickness, &
                                                   dt_seconds)
                return
            end if

            bare_ice_fluxes = bare_ice_ablation_mass(c,forc,dt_seconds,albedo)

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
        ! The schemes update the SNOW albedo; the energy balance then sees it
        ! blended with the background by the snow-cover fraction (D40).

        if (use_prescribed_albedo) then
            albedo_snow = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
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
                                      coszm,cloud,z_sur_std,dust_con,albedo_snow)

            w_snow_old = w_snow
        else if (c%albedo_scheme .eq. CHION_ALBEDO_AGING) then
            call albedo_update_aging(mass,temperature,n,c,forc%snowfall_rate,forc%dt_days, &
                                     albedo_snow,snow_age)
        else
            call albedo_update(mass,mass_w,density,temperature,n,c,forc%dt_days,albedo_snow)
        end if

        albedo = bessi_surface_albedo(mass,mass_w,density,n,c,forc,use_prescribed_albedo, &
                                      albedo_snow,alb_bg)

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
        ! The ice substrate, if any, is solved below the snow in the same
        ! system: the base of the snow conducts into it.

        lh_coef = diagnose_latent_heat_flux_coefficients(has_surface_snow,c, &
                                                         forc%air_temperature, &
                                                         forc%snowfall_rate, &
                                                         forc%rainfall_rate)

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,albedo, &
                              lh_coef%linear,lh_coef%constant,dt_seconds,energy, &
                              ice_temperature,par%ice_substrate_top_thickness)

        ! === Step 10: post-solve surface vapor mass flux =====================
        ! Evaluated exactly at the NEW interface temperature t_srf (Chion.jl
        ! 03bb445), deliberately inconsistent with the linearization used
        ! inside the solve (trap 2).

        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n, &
                                                runoff,t_srf,albedo_snow,c,forc,dt_seconds, &
                                                par%mass_split,surface_mass_min, &
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
        !
        ! A land column has no ice to charge (D41): its melt is the snow
        ! actually melted, and the shortfall energy is dropped, as on
        ! snow-free land.

        if (energy%needs_melt) then

            melt_mass = real(energy%melt_energy_available/real(c%Lm,wp_acc),wp)

            call apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo_snow, &
                            par%mass_split,surface_mass_min,melt_mass,c,melted)

            if (is_land) then
                melt = melt + melted
            else
                if (melted .lt. real(melt_mass,wp_acc)) then
                    if (n .eq. 0) then
                        ice_melt = real(melt_mass,wp_acc) - melted
                        smb_ice  = smb_ice - ice_melt
                        runoff   = runoff  + ice_melt
                    end if
                end if

                melt = melt + real(melt_mass,wp_acc)
            end if

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

        ! === Step 15: near-surface remesh ====================================
        ! Melt, sublimation, densification and refreezing can leave a thin or
        ! thick top cell: restore the fine-layer geometry before the next
        ! energy solve, conserving every column reservoir. Unconditional, as
        ! in Julia (steps 12-14 above are guarded on liquid water).

        call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                        mass_base,smb_ice,runoff,t_srf,albedo_snow, &
                                        par%Ntot,par%mass_max,par%mass_split,par%mass_min, &
                                        par%near_surface_layer_max_thicknesses,c)

        ! === Step 16: final albedo fixup =====================================
        ! The column may have gone bare during melt or sublimation, in which
        ! case the albedo diagnosed in step 6 no longer describes the surface.
        ! A bare column carries zero snow age; Chion.jl resets it only under
        ! the aging scheme, the only one that advances it, so this is the same.
        ! The written albedo is the blend on the final column (D40).

        if (use_prescribed_albedo) then
            albedo_snow = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
        else if (.not. surface_has_snow(mass,n)) then
            albedo_snow = alb_bg
        end if

        albedo = bessi_surface_albedo(mass,mass_w,density,n,c,forc,use_prescribed_albedo, &
                                      albedo_snow,alb_bg)

        if (.not. surface_has_snow(mass,n)) snow_age = 0.0_wp

        end associate

        return

    end subroutine bessi_column_step_core

    pure function bessi_surface_albedo(mass,mass_w,density,n,c,forc,use_prescribed_albedo, &
                                       albedo_snow,alb_bg) result(alb)
        ! The albedo the surface energy balance sees (chion, D40):
        !
        !   prescribed        -> the prescribed albedo (albedo_snow holds it)
        !   no surface snow   -> alb_bg
        !   otherwise         -> f*albedo_snow + (1 - f)*alb_bg
        !
        ! with the snow-cover fraction f
        !   albedo_scheme = semix: CLIMBER-X's tanh(h_snow/(c_fsnow*z0m_ice)),
        !       times its orography factor when the host gives z_sur_std;
        !       h_snow = solid thickness of the snow layers;
        !   dynamic, aging, constant: min(1, SWE/swe_crit_albedo), SWE = solid
        !       plus liquid mass of the snow layers. The column's, not the
        !       surface layer's: fine near-surface layers hold mass(1) near
        !       6 kg m-2 on any column (review Q15).
        ! Without ALBEDO_THIN_SNOW_BLEND (legacy_chion) f = 1: Chion.jl's
        ! switch from the snow albedo to the bare one.

        implicit none

        real(wp),                       intent(IN) :: mass(:)      ! (Ntot) [kg m-2]
        real(wp),                       intent(IN) :: mass_w(:)    ! (Ntot) [kg m-2]
        real(wp),                       intent(IN) :: density(:)   ! (Ntot) [kg m-3]
        integer,                        intent(IN) :: n
        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        logical,                        intent(IN) :: use_prescribed_albedo
        real(wp),                       intent(IN) :: albedo_snow  ! [1]
        real(wp),                       intent(IN) :: alb_bg       ! [1]
        real(wp) :: alb

        ! Local variables
        real(wp) :: f

        if (use_prescribed_albedo) then
            alb = albedo_snow
            return
        end if

        if (.not. surface_has_snow(mass,n)) then
            alb = alb_bg
            return
        end if

        if (.not. ALBEDO_THIN_SNOW_BLEND) then
            f = 1.0_wp
        else if (c%albedo_scheme .eq. CHION_ALBEDO_SEMIX) then
            f = semix_snow_cover_fraction(sum(mass(1:n)/density(1:n)),forc%z_sur_std, &
                                          forc%has_z_sur_std,c)
        else
            f = snow_cover_fraction(sum(mass(1:n)) + sum(mass_w(1:n)),c%swe_crit_albedo)
        end if

        alb = thin_snow_albedo(f,albedo_snow,alb_bg)

        return

    end function bessi_surface_albedo

    subroutine bessi_bare_ice_substrate_step(mass,density,temperature,n,ice_temperature, &
                                             t_srf,albedo,smb_ice,runoff,melt,vapor_mass, &
                                             sublimation,lhf_sum,c,forc, &
                                             ice_top_thickness,dt_seconds)
        ! Chion.jl _step_bare_ice_substrate! (src/step.jl:309-366 at 9ec6cc7).
        !
        ! Bare ice over a thermal substrate. The top substrate layer is the
        ! surface: it cools when the surface energy balance is negative and
        ! must be re-warmed to T0 before it melts -- the cold content of bare
        ! ice -- instead of being held at T0 with the negative energy discarded
        ! (bare_ice_ablation_mass). The full Robin solve of snow_energy_flux
        ! runs over the substrate rows; rain enthalpy enters the surface
        ! balance as on snow (has_surface_snow = .TRUE. in the precipitation
        ! coefficients); the latent flux is re-evaluated at the resolved Ts.
        ! The albedo has been set by the caller (bare-ice or prescribed).
        !
        ! Mass: melt = melt energy/Lm; vapor = Q_lh*dt/(Lv+Lm), bare ice being
        ! solid; smb_ice += vapor - melt; runoff += melt. No percolation, no
        ! refreezing. The latent heat in the energy flux is the phase's at Ts
        ! (Lv at T0) under the BESSI turbulence, as in Julia; under the semix
        ! turbulence it is Lv+Lm, consistent with the mass (D35; Julia's Lv
        ! under legacy_chion).
        !
        ! RAIN IS NOT ADDED TO RUNOFF HERE (D29). Julia adds the step's rain
        ! (runoff += rain + melt), but apply_accumulation has already routed
        ! it: to runoff when no layer can hold it, to mass_w(1) under a sliver
        ! surface layer (0 < mass(1) <= TOL_EMPTY_LAYER), where Julia counts it
        ! twice.
        !
        ! A sliver surface layer keeps its rows in the solve above the
        ! substrate, as in Julia: snow_energy_flux counts layers with
        ! mass(1) > 0, not TOL_EMPTY_LAYER.

        implicit none

        real(wp),     intent(IN)    :: mass(:)            ! (Ntot) [kg m-2]
        real(wp),     intent(IN)    :: density(:)         ! (Ntot) [kg m-3]
        real(wp),     intent(INOUT) :: temperature(:)     ! (Ntot) [K]
        integer,      intent(IN)    :: n                  ! active layer count
        real(wp),     intent(INOUT) :: ice_temperature(:) ! (n_ice) [K]
        real(wp),     intent(INOUT) :: t_srf              ! [K]
        real(wp),     intent(IN)    :: albedo             ! [1]
        real(wp_acc), intent(INOUT) :: smb_ice            ! [kg m-2]
        real(wp_acc), intent(INOUT) :: runoff             ! [kg m-2]
        real(wp_acc), intent(INOUT) :: melt               ! [kg m-2]
        real(wp_acc), intent(INOUT) :: vapor_mass         ! [kg m-2]
        real(wp_acc), intent(INOUT) :: sublimation        ! [kg m-2]
        real(wp_acc), intent(INOUT) :: lhf_sum            ! [W m-2 d]

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc

        real(wp),     intent(IN)    :: ice_top_thickness  ! [m]
        real(wp),     intent(IN)    :: dt_seconds         ! [s]

        ! Local variables
        real(wp) :: q_lh, vapor, melt_mass

        type(latent_heat_coeff_class)  :: lh_coef
        type(snow_energy_result_class) :: energy

        lh_coef = diagnose_latent_heat_flux_coefficients(.TRUE.,c,forc%air_temperature, &
                                                         forc%snowfall_rate, &
                                                         forc%rainfall_rate)

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,albedo, &
                              lh_coef%linear,lh_coef%constant,dt_seconds,energy, &
                              ice_temperature,ice_top_thickness)

        ! Julia re-evaluates all non-shortwave components at the resolved Ts
        ! (_resolved_nonshortwave_surface_flux_components); only the latent
        ! one feeds a budget here. A bare-ice surface: no snow depth for the
        ! CLIMBER-X roughness blend, the ice roughness and (D35) latent heat
        ! for the semix turbulence.
        q_lh = resolved_turbulent_latent_heat_flux(c,forc,t_srf,0.0_wp,.FALSE.)

        vapor = real(real(q_lh,wp_acc)*real(dt_seconds,wp_acc) &
                     /real(c%Lv + c%Lm,wp_acc),wp)

        melt_mass = 0.0_wp
        if (energy%needs_melt) &
            melt_mass = real(energy%melt_energy_available/real(c%Lm,wp_acc),wp)

        ! Accumulate directly into the wp_acc accumulators, in Julia's order.
        smb_ice     = smb_ice + real(vapor,wp_acc) - real(melt_mass,wp_acc)
        melt        = melt   + real(melt_mass,wp_acc)
        runoff      = runoff + real(melt_mass,wp_acc)
        vapor_mass  = vapor_mass  + real(vapor,wp_acc)
        sublimation = sublimation + real(max(-vapor,0.0_wp),wp_acc)

        ! W m-2 DAYS, not W m-2 seconds. dt_days, not dt_seconds.
        lhf_sum = lhf_sum + real(q_lh,wp_acc)*real(forc%dt_days,wp_acc)

        return

    end subroutine bessi_bare_ice_substrate_step

    subroutine bessi_column_step(bsi,icol,forc_in,c)
        ! Chion.jl column_step! (src/step.jl:113-157): the diurnal-shortwave
        ! substep wrapper around the core kernel.
        !
        ! The cloud-proxy longwave (longwave_scheme = "cloud_proxy", Chion.jl
        ! 03bb445) is resolved first, once, from the DAILY forcing, and enters
        ! the substeps as a prescribed downwelling flux.
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
        !
        ! A day the criterion does not split (n_substeps = 1) is stepped with
        ! its forcing as given. Chion.jl runs it as one [-pi, pi] interval,
        ! whose shortwave average is the daily mean again -- except where the
        ! geometry has no daylight (polar night, no latitude), where it is
        ! zero whatever the forcing says. legacy_chion does the same
        ! (DIURNAL_SINGLE_INTERVAL_AVERAGED, docs/porting_notes.md D39).

        implicit none

        type(bessi_class),              intent(INOUT) :: bsi
        integer,                        intent(IN)    :: icol
        type(chion_step_forcing_class), intent(IN)    :: forc_in
        type(chion_const_class),        intent(IN)    :: c

        ! Local variables
        integer  :: n_substeps, k
        real(wp) :: shortwave_for_criterion
        real(wp) :: hour_angle_start, hour_angle_end, fraction
        real(wp) :: amplitude

        type(chion_step_forcing_class) :: forc, subforc

        forc = with_parameterized_longwave(c,forc_in)

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

            if (n_substeps .gt. 1 .or. DIURNAL_SINGLE_INTERVAL_AVERAGED) then

                ! Elevation-dependent half-amplitude (Chion.jl d0146e1); the
                ! same for every interval of the day.
                amplitude = diurnal_temperature_amplitude( &
                                bsi%par%diurnal_temperature_amplitude, &
                                bsi%par%diurnal_temperature_amplitude_gradient, &
                                bsi%par%diurnal_temperature_amplitude_reference_height, &
                                bsi%par%diurnal_temperature_amplitude_max, &
                                forc%surface_height)

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
                                                amplitude,hour_angle_start,hour_angle_end)
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
