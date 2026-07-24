module chion_defs
    ! Core definitions for chion: precision, tolerances, physical constants,
    ! scheme flags, and the parameter/forcing/grid derived types.
    !
    ! Fortran port of Chion.jl (https://github.com/fesmc/Chion.jl), branch main.
    ! Source of truth for this module:
    !   Chion.jl/src/constants.jl   - tolerances, defaults, SnowpackPhysicalConstants
    !   Chion.jl/src/forcing.jl     - SnowpackForcing, SnowpackStepForcing
    !   Chion.jl/src/domain.jl      - SnowpackGrid
    !   Chion.jl/src/models.jl      - BESSIModel configuration
    !
    ! This module contains no physics. See docs/PLAN.md (WP1).

    use, intrinsic :: iso_fortran_env, only : error_unit
    use precision, only : sp, dp
    use nml, only : nml_replace

    implicit none

    private

    ! === Precision ===========================================================
    !
    ! PRECISION POLICY  (measured, not assumed -- see docs/porting_notes.md D1)
    !
    ! wp = sp for all state, forcing and interfaces. This matches yelmo and
    ! fesm-utils, so no conversion is needed at the yelmox boundary. It was
    ! verified adequate for: the layered state (mass, mass_w, density,
    ! temperature), the tridiagonal conduction solve (strictly diagonally
    ! dominant; max error 3e-5 K), the densification density gap, and the
    ! melt-energy residual.
    !
    ! wp_acc = dp is MANDATORY for cumulative accumulators. These are summed
    ! every step and never reset, so in sp small increments onto a large total
    ! are lost outright: 0.01 kg m-2 increments onto 1e5 kg m-2 drift by
    ! 80 kg m-2 (0.08%) over a 100-year daily run. Applies to smb_ice, runoff,
    ! melt, refreezing, vapor_mass, sublimation, latent_heat_flux_sum,
    ! mass_base and pdd_sum.
    !
    ! Three expressions must additionally be evaluated in dp LOCALLY, even
    ! though their inputs and outputs are wp. Each is a difference of nearly
    ! equal numbers feeding a division or a tolerance test:
    !   1. pore volume   phi = m/rho - m/rho_i     (percolation, albedo)
    !      In sp, phi is quantized at ~1e-5 m, so the TOL_TINY guard can never
    !      fire and the division lwc = m_w/rho_w/phi is unstable.
    !   2. cold content  (T0 - T)*ci*m_s           (refreezing)
    !      sp resolution at 273 K is 3e-5 K; smaller offsets are pure noise.
    !   3. surface energy accumulation over diurnal substeps.
    ! Use real(wp_acc) locals, convert back on store. This is a few operations
    ! per column per step -- the cost is not measurable.

    ! wp is selectable at COMPILE TIME: `make ... precision=dp` defines
    ! CHION_DP. sp is the production setting and the yelmo/yelmox interface
    ! kind; dp exists so that the port can be compared against Chion.jl, which
    ! is Float64 throughout, without sp round-off being mistaken for a porting
    ! error -- and so that the sp-vs-dp difference is measured rather than
    ! assumed. See docs/porting_notes.md D19 and validation/.
    !
    ! Note this moves ONLY wp. wp_acc stays dp in both builds: the accumulator
    ! argument below is about summation over 10^4-10^5 steps, not about the
    ! precision of the state, and it holds regardless of what wp is.
#ifdef CHION_DP
    integer, parameter, public :: wp     = dp
#else
    integer, parameter, public :: wp     = sp
#endif
    integer, parameter, public :: wp_acc = dp      ! cumulative accumulators
    integer, parameter, public :: wp_chion = wp    ! exposed to external models

    public :: sp, dp

    ! === IO units and sentinels ==============================================

    integer,  parameter, public :: io_unit_err = error_unit
    real(wp), parameter, public :: MV          = -9999.0_wp

    ! === Numeric tolerances ==================================================
    ! Chion.jl/src/constants.jl:9-10.
    !
    ! WARNING 1: these two are NOT interchangeable. Different routines gate on
    ! different thresholds, deliberately. See docs/PLAN.md section 5, item 1.
    !
    ! WARNING 2: both are declared dp so that comparisons promote correctly.
    ! TOL_TINY is below sp resolution for any quantity of order 1 or larger,
    ! so a guard of the form "x <= TOL_TINY" is only meaningful when x itself
    ! was computed in dp. See the precision policy above.

    real(wp_acc), parameter, public :: TOL_TINY        = 1.0e-12_wp_acc  ! EPS_TINY
    real(wp_acc), parameter, public :: TOL_EMPTY_LAYER = 1.0e-10_wp_acc  ! EPS_EMPTY_LAYER

    ! === Scheme flags ========================================================
    ! Chion.jl/src/constants.jl:43-49. Integer-valued so that physics routines
    ! can branch on them without string comparisons in inner loops. Namelist
    ! input uses the string names; see chion_*_scheme_flag below.

    integer, parameter, public :: CHION_FRESH_SNOW_DENSITY_CONSTANT      = 1
    integer, parameter, public :: CHION_FRESH_SNOW_DENSITY_PARAMETERIZED = 2

    integer, parameter, public :: CHION_ALBEDO_CONSTANT   = 1
    integer, parameter, public :: CHION_ALBEDO_DYNAMIC    = 2
    integer, parameter, public :: CHION_ALBEDO_PRESCRIBED = 3
    integer, parameter, public :: CHION_ALBEDO_SEMIX      = 4

    ! Which spectral snow-albedo parameterization the SEMIX scheme uses.
    ! CLIMBER-X defaults to Dang (its isnow_albedo = 2).
    integer, parameter, public :: SEMIX_SNOW_ALBEDO_WW   = 1
    integer, parameter, public :: SEMIX_SNOW_ALBEDO_DANG = 2

    ! Surface energy balance family. "bessi" is the bulk-coefficient scheme
    ! chion was ported with; "semix" is CLIMBER-X SEMIX's aerodynamic scheme.
    ! Orthogonal to albedo_scheme and to Ntot -- see docs/semix_port_scope.md.
    integer, parameter, public :: CHION_SEB_BESSI = 1
    integer, parameter, public :: CHION_SEB_SEMIX = 2

    ! Which saturation-specific-humidity parameterization the SEMIX surface
    ! scheme uses for the turbulent latent flux. SEMIX is CLIMBER-X's own
    ! q_sat_i/dqsat_dT_i; BESSI routes chion's ice saturation vapour pressure
    ! through the same 0.622/p conversion.
    integer, parameter, public :: SEMIX_QSAT_SEMIX = 1
    integer, parameter, public :: SEMIX_QSAT_BESSI = 2

    integer, parameter, public :: CHION_DENSIFY_BESSI   = 1
    integer, parameter, public :: CHION_DENSIFY_HTESSEL = 2

    ! === Defaults ============================================================
    ! Chion.jl/src/constants.jl:16-40.

    real(wp), parameter, public :: DEF_SECONDS_PER_DAY = 86400.0_wp

    real(wp), parameter, public :: DEF_SEA_LEVEL_AIR_PRESSURE = 101325.0_wp
    real(wp), parameter, public :: DEF_GRAVITY                = 9.80665_wp
    real(wp), parameter, public :: DEF_MOLAR_MASS_DRY_AIR     = 0.0289644_wp
    real(wp), parameter, public :: DEF_UNIVERSAL_GAS_CONSTANT = 8.31446261815324_wp

    ! === Reference-reproduction mode =========================================
    !
    ! `make ... legacy_chion=1` defines CHION_LEGACY, which reverts the
    ! DELIBERATE PHYSICS CORRECTIONS chion has made against Chion.jl, so that
    ! the validation harness can still prove the port is faithful.
    !
    ! It exists because "is the port correct?" and "is the reference correct?"
    ! are different questions and must not be conflated. Without it, every
    ! upstream bug chion fixes would show up as a WP16 gate failure, and the
    ! only way to keep the gate green would be to stop testing those fields --
    ! so the harness would get weaker exactly as the port got better.
    !
    ! THIS IS NOT A PRODUCTION SETTING. It selects physics believed to be
    ! wrong. Nothing but validation/ should ever build with it.
    !
    ! Currently reverted under CHION_LEGACY:
    !   * DENSIFY_R_GAS  -> 8.13, Chion.jl's typo for the gas constant
    !     (Chion.jl issue #18, docs/porting_notes.md D22).
    !   * DENSIFY_GRAVITY -> 9.81, Chion.jl's second gravity constant
    !     (docs/porting_notes.md D25).
    !
    ! Deliberately NOT covered: the PDD smb_ice convention (Chion.jl issue #19,
    ! D23). Chion.jl's PDD is not authoritative (docs/PLAN.md section 3.2), so
    ! reproducing its convention would mean maintaining a second PDD core --
    ! which is upstream defect 13 (three diverged copies) reintroduced on
    ! purpose. PDD is compared to Chion.jl as a REPORTED diagnostic instead,
    ! and gated on its own mass-closure identity.
#ifdef CHION_LEGACY
    real(wp_acc), parameter, public :: DENSIFY_R_GAS   = 8.13_wp_acc
    real(wp_acc), parameter, public :: DENSIFY_GRAVITY = 9.81_wp_acc
    logical,      parameter, public :: CHION_LEGACY_MODE = .TRUE.
#else
    real(wp_acc), parameter, public :: DENSIFY_R_GAS   = real(DEF_UNIVERSAL_GAS_CONSTANT,wp_acc)
    real(wp_acc), parameter, public :: DENSIFY_GRAVITY = real(DEF_GRAVITY,wp_acc)
    logical,      parameter, public :: CHION_LEGACY_MODE = .FALSE.
#endif

    integer,  parameter, public :: DEF_NTOT             = 15
    real(wp), parameter, public :: DEF_MASS_MAX         = 500.0_wp
    real(wp), parameter, public :: DEF_MASS_SPLIT       = 300.0_wp
    real(wp), parameter, public :: DEF_MASS_MIN         = 100.0_wp
    real(wp), parameter, public :: DEF_DENSITY_INIT     = 300.0_wp
    real(wp), parameter, public :: DEF_TEMPERATURE_INIT = 273.0_wp

    ! Depth cap reference values. NOTE: the cap uses this hard-coded layer
    ! count, NOT the configured Ntot. Chion.jl/src/constants.jl:37-38 and
    ! processes/layer_structure.jl. See docs/PLAN.md section 5, item 11.
    integer,  parameter, public :: BESSI_REFERENCE_LAYER_COUNT   = 15
    real(wp), parameter, public :: BESSI_REFERENCE_DEPTH_DENSITY = 300.0_wp

    ! === Physical constants ==================================================

    type chion_const_class
        ! Mirrors Chion.jl SnowpackPhysicalConstants (src/constants.jl:61-88),
        ! field for field and in the same order. Defaults at src/constants.jl:173-200.
        !
        ! Porting note: Chion.jl keeps the three scheme flags inside this struct.
        ! chion does the same, so that every physics routine (which already
        ! receives c) can branch without threading an extra argument.

        ! Densities
        real(wp) :: rho_s              ! [kg m-3] fresh snow density (constant scheme)
        real(wp) :: rho_i              ! [kg m-3] ice density
        real(wp) :: rho_w              ! [kg m-3] water density

        ! Parameterized fresh-snow density: rho = a + b*(T-T0) + c*sqrt(wind)
        real(wp) :: rho_s_a            ! [kg m-3]
        real(wp) :: rho_s_b            ! [kg m-3 K-1]
        real(wp) :: rho_s_c            ! [kg m-3 (m s-1)^-1/2]
        integer  :: fresh_snow_density_scheme   ! CHION_FRESH_SNOW_DENSITY_*

        ! Thermal properties
        real(wp) :: Ki                 ! [W m-1 K-1] thermal conductivity of ice
        real(wp) :: ci                 ! [J kg-1 K-1] heat capacity of ice
        real(wp) :: cw                 ! [J kg-1 K-1] heat capacity of water
        real(wp) :: Lm                 ! [J kg-1] latent heat of melting
        real(wp) :: Lv                 ! [J kg-1] latent heat of vaporization
        real(wp) :: cp_air             ! [J kg-1 K-1] heat capacity of air
        real(wp) :: latent_heat_flux_ratio  ! [1] scaling of turbulent latent flux

        ! Turbulent exchange
        real(wp) :: D_sh               ! [W m-2 K-1] sensible heat exchange coefficient

        ! Surface energy balance scheme, and the aerodynamic exchange it needs
        ! (CHION_SEB_SEMIX only). Roughness lengths and the surface-layer height
        ! are CLIMBER-X smb_par / constants values; karman, grav and R_dry are
        ! the universal constants SEMIX pulls from its constants module. The
        ! heat capacity of air is chion's existing cp_air (1003 vs SEMIX's
        ! 1000 J kg-1 K-1, 0.3% on f_sh) rather than a second constant for the
        ! same quantity.
        integer  :: seb_scheme         ! CHION_SEB_*
        real(wp) :: z0m_snow           ! [m] momentum roughness length, snow
        real(wp) :: z0m_ice            ! [m] momentum roughness length, ice
        real(wp) :: zm_to_zh           ! [1] heat/momentum roughness ratio
        real(wp) :: z_sfl              ! [m] surface layer height
        real(wp) :: karman             ! [1] von Karman constant
        real(wp) :: grav               ! [m s-2] gravitational acceleration
        real(wp) :: R_dry              ! [J kg-1 K-1] gas constant of dry air
        logical  :: l_neutral          ! [1] force neutral stratification
        logical  :: l_dew              ! [1] allow dew/frost deposition
        integer  :: semix_qsat         ! SEMIX_QSAT_*

        ! Albedo
        real(wp) :: alpha_dry          ! [1] dry snow albedo (upper bound)
        real(wp) :: alpha_wet          ! [1] wet snow albedo (lower bound)
        real(wp) :: alpha_ice          ! [1] bare ice albedo
        real(wp) :: max_lwc_albedo     ! [1] LWC at which albedo reaches alpha_wet
        integer  :: albedo_scheme      ! CHION_ALBEDO_*

        ! SEMIX spectral albedo (CHION_ALBEDO_SEMIX). Warren & Wiscombe 1980
        ! bands, collapsed to broadband by the incoming-SW spectral weights.
        real(wp) :: frac_vu            ! [1]  visible+UV fraction of incoming solar
        real(wp) :: alb_snow_vis_new   ! [1]  fresh-snow visible diffuse albedo
        real(wp) :: alb_snow_nir_new   ! [1]  fresh-snow near-IR diffuse albedo
        real(wp) :: snow_grain_fresh   ! [um] fresh snow grain size
        real(wp) :: snow_grain_old     ! [um] aged snow grain size
        real(wp) :: d_alb_age_vis      ! [1]  visible aging albedo reduction
        real(wp) :: d_alb_age_nir      ! [1]  near-IR aging albedo reduction
        real(wp) :: f_age_t            ! [K-1] dry-snow temperature aging factor
        real(wp) :: dT_age             ! [K]  aging temperature offset
        real(wp) :: snow_0             ! [kg m-2 d-1] critical snowfall rate for aging
        real(wp) :: snow_1             ! [1]  snowfall-rate aging exponent
        real(wp) :: w_snow_dust        ! [kg m-2] SWE melt that doubles dust concentration
        real(wp) :: dust_con_scale     ! [1]  dust concentration scaling
        integer  :: semix_snow_albedo  ! SEMIX_SNOW_ALBEDO_*
        real(wp) :: dalb_snow_vis      ! [1]  visible snow albedo offset (Dang)
        real(wp) :: dalb_snow_nir      ! [1]  near-IR snow albedo offset (Dang)
        real(wp) :: k_sigma_orog       ! [1]  orographic albedo reduction scale (Dang)
        real(wp) :: sigma_orog_crit    ! [m]  orographic roughness scale (Dang)

        ! Radiation. eps_ice is consulted ONLY by seb_scheme = semix, which
        ! carries SEMIX's snow/ice emissivity pair; the bessi scheme applies
        ! eps_snow to bare ice as well, as Chion.jl does.
        real(wp) :: eps_air            ! [1] emissivity of air
        real(wp) :: eps_snow           ! [1] emissivity of snow
        real(wp) :: eps_ice            ! [1] emissivity of bare ice (semix SEB)
        real(wp) :: sigma_sb           ! [W m-2 K-4] Stefan-Boltzmann constant

        ! Reference values
        real(wp) :: T0                 ! [K] freezing point of water
        real(wp) :: seconds_per_day    ! [s]

        ! Densification
        integer  :: low_density_densification   ! CHION_DENSIFY_*
    end type chion_const_class

    ! === Per-column, per-substep forcing =====================================

    type chion_step_forcing_class
        ! Mirrors Chion.jl SnowpackStepForcing (src/forcing.jl:234-257) field
        ! for field and in the same order. This is the contract that makes the
        ! snowpack models interchangeable: every model kernel takes exactly
        ! this type and ignores the fields it does not need.
        !
        ! The has_* flags select prescribed fluxes over internal
        ! parameterizations. See docs/PLAN.md section 1.

        real(wp) :: air_temperature     ! [K]
        real(wp) :: dt_days             ! [d]
        real(wp) :: snowfall_rate       ! [kg m-2 s-1]
        real(wp) :: rainfall_rate       ! [kg m-2 s-1]
        real(wp) :: shortwave_down      ! [W m-2]
        real(wp) :: wind_speed          ! [m s-1]

        real(wp) :: q_sw_net = 0.0_wp            ! [W m-2] net shortwave, if prescribed
        real(wp) :: q_lw_down = 0.0_wp           ! [W m-2] downward longwave, if prescribed
        real(wp) :: q_sh = 0.0_wp                ! [W m-2] sensible heat flux, if prescribed
        real(wp) :: q_lh = 0.0_wp                ! [W m-2] latent heat flux, if prescribed

        logical  :: has_q_sw_net = .FALSE.
        logical  :: has_q_lw_down = .FALSE.
        logical  :: has_q_sh = .FALSE.
        logical  :: has_q_lh = .FALSE.

        real(wp) :: relative_humidity = 0.0_wp   ! [1] or [%]; >1 is interpreted as percent
        logical  :: has_relative_humidity = .FALSE.

        real(wp) :: air_pressure        ! [Pa]
        real(wp) :: prescribed_albedo = 0.0_wp   ! [1]
        logical  :: has_prescribed_albedo = .FALSE.

        ! SEMIX albedo inputs. coszm falls back to a lat/solar-longitude daily
        ! mean when absent; cloud defaults to clear-sky (all-direct) when absent.
        real(wp) :: coszm = 0.0_wp               ! [1] daily-mean cos(solar zenith)
        logical  :: has_coszm = .FALSE.
        real(wp) :: cloud = 0.0_wp               ! [1] cloud fraction, [0,1]
        logical  :: has_cloud = .FALSE.
        real(wp) :: dust_dep = 0.0_wp            ! [kg m-2 s-1] dust deposition rate
        logical  :: has_dust_dep = .FALSE.
        real(wp) :: z_sur_std = 0.0_wp           ! [m] subgrid surface-height std deviation
        logical  :: has_z_sur_std = .FALSE.
        real(wp) :: alb_ice_host = 0.0_wp        ! [1] bare-ice albedo supplied by the host
        logical  :: has_alb_ice_host = .FALSE.

        real(wp) :: latitude_deg        ! [deg N]
        real(wp) :: day_of_year         ! [d] fractional, 1-based
        real(wp) :: solar_longitude_deg ! [deg]
    end type chion_step_forcing_class

    ! === Host-facing forcing =================================================

    type chion_forcing_class
        ! Per-column forcing arrays, written directly by the host between
        ! calls to chion_update. Corresponds to one time slice of Chion.jl's
        ! SnowpackForcing (src/forcing.jl:180-224), which stores (ncol,ntime).

        integer :: ncol

        real(wp), allocatable :: air_temperature(:)      ! [K]
        real(wp), allocatable :: snowfall_rate(:)        ! [kg m-2 s-1]
        real(wp), allocatable :: rainfall_rate(:)        ! [kg m-2 s-1]
        real(wp), allocatable :: shortwave_down(:)       ! [W m-2]
        real(wp), allocatable :: wind_speed(:)           ! [m s-1]

        real(wp), allocatable :: q_sw_net(:)             ! [W m-2]
        real(wp), allocatable :: q_lw_down(:)            ! [W m-2]
        real(wp), allocatable :: q_sh(:)                 ! [W m-2]
        real(wp), allocatable :: q_lh(:)                 ! [W m-2]

        logical,  allocatable :: has_q_sw_net(:)
        logical,  allocatable :: has_q_lw_down(:)
        logical,  allocatable :: has_q_sh(:)
        logical,  allocatable :: has_q_lh(:)

        real(wp), allocatable :: relative_humidity(:)    ! [1]
        logical,  allocatable :: has_relative_humidity(:)

        real(wp), allocatable :: surface_height(:)       ! [m] used for air pressure
        real(wp), allocatable :: air_pressure(:)         ! [Pa]
        real(wp), allocatable :: prescribed_albedo(:)    ! [1]
        logical,  allocatable :: has_prescribed_albedo(:)

        real(wp), allocatable :: coszm(:)                ! [1] daily-mean cos(zenith)
        logical,  allocatable :: has_coszm(:)
        real(wp), allocatable :: cloud(:)                ! [1] cloud fraction
        logical,  allocatable :: has_cloud(:)
        real(wp), allocatable :: dust_dep(:)             ! [kg m-2 s-1] dust deposition
        logical,  allocatable :: has_dust_dep(:)
        real(wp), allocatable :: z_sur_std(:)            ! [m] subgrid height std dev
        logical,  allocatable :: has_z_sur_std(:)
        real(wp), allocatable :: alb_ice_host(:)         ! [1] host bare-ice albedo
        logical,  allocatable :: has_alb_ice_host(:)

        real(wp), allocatable :: latitude_deg(:)         ! [deg N]

        ! --- ITM-only fields (WP11) ------------------------------------
        !
        ! These three are deliberately NOT part of chion_step_forcing_class.
        ! That type mirrors Chion.jl's SnowpackStepForcing and is the shared,
        ! model-neutral contract every kernel takes; adding ice-sheet state to
        ! it would make BESSI and PDD carry fields they can never use. ITM
        ! instead receives them as explicit arguments from the dispatcher
        ! (itm_step(itm,icol,fc,z_srf,H_ice,PDDs)).
        !
        ! ITM's z_srf is the EXISTING surface_height(:) field above -- there is
        ! no separate array for it.
        !
        ! PDDs is a WHOLE-YEAR total, not a per-step value. smbpal recomputes
        ! it once per year from the annual temperature series
        ! (smbpal.f90: calc_pdds / the annual loop) and holds it fixed for
        ! every step of that year. It is used only to interpolate the critical
        ! snow depth between the "desert" and "forest" end members in
        ! calc_albedo_surface. A host that overwrites it every step with a
        ! per-step degree-day increment will get the desert branch always, and
        ! a systematically different albedo. Set it once per year.

        real(wp), allocatable :: H_ice(:)                ! [m] ice thickness
        real(wp), allocatable :: PDDs(:)                 ! [K d] ANNUAL positive degree days

        ! Time metadata, uniform across columns
        real(wp) :: day_of_year                          ! [d] fractional, 1-based
        real(wp) :: solar_longitude_deg                  ! [deg]
    end type chion_forcing_class

    ! === Grid ================================================================

    type chion_grid_class
        ! Chion.jl SnowpackGrid (src/domain.jl:20-27): a packed list of ncol
        ! independent columns, with an optional mapping back onto a 2-D grid
        ! for IO only. The physics never uses the spatial coordinates.

        integer :: ncol                                  ! total columns
        integer :: n_active                              ! currently active columns

        logical :: has_spatial                           ! are x/y/js/is/mask set?
        real(wp), allocatable :: x(:)                    ! [m or deg] grid x axis
        real(wp), allocatable :: y(:)                    ! [m or deg] grid y axis
        integer,  allocatable :: js(:)                   ! (ncol) y index of each column
        integer,  allocatable :: is(:)                   ! (ncol) x index of each column
        real(wp), allocatable :: mask(:,:)               ! (ny,nx) domain mask

        logical,  allocatable :: active(:)               ! (ncol) column on/off
        integer,  allocatable :: active_idx(:)           ! (n_active) packed active list
    end type chion_grid_class

    ! === Model parameters ====================================================

    type chion_param_class
        ! Top-level configuration. Model-specific parameters live in the
        ! model modules (bessi_par_class, pdd_par_class, itm_par_class).

        character(len=56)  :: model          ! "bessi" | "pdd" | "itm"

        ! Namelist group names, so chion can be instantiated more than once
        character(len=56)  :: nml_chion
        character(len=56)  :: nml_bessi
        character(len=56)  :: nml_pdd
        character(len=56)  :: nml_itm
        character(len=56)  :: nml_const

        character(len=512) :: phys_const_file ! path to the physical constants file
        character(len=512) :: phys_const     ! group name within phys_const_file
        character(len=512) :: restart        ! restart file, or "none"

        logical :: use_omp                   ! set at init from OpenMP availability
    end type chion_param_class

    ! === Public interface ====================================================

    public :: chion_const_class
    public :: chion_step_forcing_class
    public :: chion_forcing_class
    public :: chion_grid_class
    public :: chion_param_class

    public :: chion_const_init
    public :: chion_const_print

    public :: chion_forcing_alloc
    public :: chion_forcing_dealloc
    public :: chion_grid_init
    public :: chion_grid_dealloc
    public :: chion_grid_set_active

    public :: chion_albedo_scheme_flag
    public :: chion_semix_snow_albedo_flag
    public :: chion_seb_scheme_flag
    public :: chion_semix_qsat_flag
    public :: chion_fresh_snow_density_scheme_flag
    public :: chion_densify_scheme_flag

    public :: chion_check_enum
    public :: chion_check_file
    public :: chion_parse_path
    public :: chion_load_command_line_args

contains

    subroutine chion_const_init(c)
        ! Populate a constants object with the Chion.jl defaults
        ! (src/constants.jl:173-200). Parameter loading from a namelist
        ! overrides these; see WP13.

        implicit none

        type(chion_const_class), intent(OUT) :: c

        c%rho_s   = 315.0_wp
        c%rho_i   = 917.0_wp
        c%rho_w   = 1000.0_wp

        c%rho_s_a = 109.0_wp
        c%rho_s_b = 6.0_wp
        c%rho_s_c = 26.0_wp
        c%fresh_snow_density_scheme = CHION_FRESH_SNOW_DENSITY_CONSTANT

        c%Ki      = 2.1_wp
        c%ci      = 2110.0_wp
        c%cw      = 4181.0_wp
        c%Lm      = 334000.0_wp
        c%Lv      = 2.501e6_wp
        c%cp_air  = 1003.0_wp
        c%latent_heat_flux_ratio = 1.0_wp

        c%D_sh    = 10.0_wp

        ! SEMIX aerodynamic exchange defaults (CLIMBER-X smb_par.nml /
        ! smb_params.f90 / constants.f90).
        c%seb_scheme  = CHION_SEB_BESSI
        c%z0m_snow    = 0.0024_wp
        c%z0m_ice     = 0.002_wp
        c%zm_to_zh    = exp(-2.0_wp)
        c%z_sfl       = 100.0_wp
        c%karman      = 0.4_wp
        c%grav        = 9.81_wp
        c%R_dry       = 287.058_wp
        c%l_neutral   = .FALSE.
        c%l_dew       = .TRUE.
        c%semix_qsat  = SEMIX_QSAT_SEMIX

        c%alpha_dry      = 0.81_wp
        c%alpha_wet      = 0.70_wp
        c%alpha_ice      = 0.30_wp
        c%max_lwc_albedo = 0.10_wp
        c%albedo_scheme  = CHION_ALBEDO_DYNAMIC

        ! SEMIX spectral albedo defaults (CLIMBER-X smb_par / constants).
        c%frac_vu          = 0.45_wp
        c%alb_snow_vis_new = 0.99_wp
        c%alb_snow_nir_new = 0.65_wp
        c%snow_grain_fresh = 50.0_wp
        c%snow_grain_old   = 1000.0_wp
        c%d_alb_age_vis    = 0.05_wp
        c%d_alb_age_nir    = 0.25_wp
        c%f_age_t          = 0.1_wp
        c%dT_age           = 0.0_wp
        c%snow_0           = 1.0_wp
        c%snow_1           = 0.5_wp
        c%w_snow_dust      = 10.0_wp
        c%dust_con_scale   = 1.0_wp
        c%semix_snow_albedo = SEMIX_SNOW_ALBEDO_DANG
        c%dalb_snow_vis     = 0.0_wp
        c%dalb_snow_nir     = 0.0_wp
        c%k_sigma_orog      = 0.0_wp
        c%sigma_orog_crit   = 1000.0_wp

        c%eps_air  = 0.80_wp
        c%eps_snow = 0.98_wp
        c%eps_ice  = 0.98_wp
        c%sigma_sb = 5.670373e-8_wp

        c%T0              = 273.15_wp
        c%seconds_per_day = DEF_SECONDS_PER_DAY

        c%low_density_densification = CHION_DENSIFY_BESSI

        return

    end subroutine chion_const_init

    subroutine chion_const_print(c)
        ! Write the full constants set to stdout, for provenance in run logs.

        implicit none

        type(chion_const_class), intent(IN) :: c

        write(*,"(a)") "chion physical constants:"
        write(*,"(a25,g14.6,a)") "rho_s   = ", c%rho_s,   "  [kg m-3]"
        write(*,"(a25,g14.6,a)") "rho_i   = ", c%rho_i,   "  [kg m-3]"
        write(*,"(a25,g14.6,a)") "rho_w   = ", c%rho_w,   "  [kg m-3]"
        write(*,"(a25,g14.6,a)") "rho_s_a = ", c%rho_s_a, "  [kg m-3]"
        write(*,"(a25,g14.6,a)") "rho_s_b = ", c%rho_s_b, "  [kg m-3 K-1]"
        write(*,"(a25,g14.6,a)") "rho_s_c = ", c%rho_s_c, "  [kg m-3 (m s-1)^-1/2]"
        write(*,"(a25,i14)")     "fresh_snow_density_scheme = ", c%fresh_snow_density_scheme
        write(*,"(a25,g14.6,a)") "Ki      = ", c%Ki,      "  [W m-1 K-1]"
        write(*,"(a25,g14.6,a)") "ci      = ", c%ci,      "  [J kg-1 K-1]"
        write(*,"(a25,g14.6,a)") "cw      = ", c%cw,      "  [J kg-1 K-1]"
        write(*,"(a25,g14.6,a)") "Lm      = ", c%Lm,      "  [J kg-1]"
        write(*,"(a25,g14.6,a)") "Lv      = ", c%Lv,      "  [J kg-1]"
        write(*,"(a25,g14.6,a)") "cp_air  = ", c%cp_air,  "  [J kg-1 K-1]"
        write(*,"(a25,g14.6,a)") "latent_heat_flux_ratio = ", c%latent_heat_flux_ratio, "  [1]"
        write(*,"(a25,g14.6,a)") "D_sh    = ", c%D_sh,    "  [W m-2 K-1]"
        write(*,"(a25,i14)")     "seb_scheme = ", c%seb_scheme
        write(*,"(a25,g14.6,a)") "z0m_snow = ", c%z0m_snow, "  [m]"
        write(*,"(a25,g14.6,a)") "z0m_ice  = ", c%z0m_ice,  "  [m]"
        write(*,"(a25,g14.6,a)") "zm_to_zh = ", c%zm_to_zh, "  [1]"
        write(*,"(a25,g14.6,a)") "z_sfl    = ", c%z_sfl,    "  [m]"
        write(*,"(a25,l14)")     "l_neutral = ", c%l_neutral
        write(*,"(a25,l14)")     "l_dew     = ", c%l_dew
        write(*,"(a25,i14)")     "semix_qsat = ", c%semix_qsat
        write(*,"(a25,g14.6,a)") "alpha_dry = ", c%alpha_dry, "  [1]"
        write(*,"(a25,g14.6,a)") "alpha_wet = ", c%alpha_wet, "  [1]"
        write(*,"(a25,g14.6,a)") "alpha_ice = ", c%alpha_ice, "  [1]"
        write(*,"(a25,g14.6,a)") "max_lwc_albedo = ", c%max_lwc_albedo, "  [1]"
        write(*,"(a25,i14)")     "albedo_scheme = ", c%albedo_scheme
        write(*,"(a25,g14.6,a)") "eps_air  = ", c%eps_air,  "  [1]"
        write(*,"(a25,g14.6,a)") "eps_snow = ", c%eps_snow, "  [1]"
        write(*,"(a25,g14.6,a)") "eps_ice  = ", c%eps_ice,  "  [1]"
        write(*,"(a25,g14.6,a)") "sigma_sb = ", c%sigma_sb, "  [W m-2 K-4]"
        write(*,"(a25,g14.6,a)") "T0       = ", c%T0,       "  [K]"
        write(*,"(a25,g14.6,a)") "seconds_per_day = ", c%seconds_per_day, "  [s]"
        write(*,"(a25,i14)")     "low_density_densification = ", c%low_density_densification

        return

    end subroutine chion_const_print

    subroutine chion_forcing_alloc(forc,ncol)
        ! Allocate all forcing arrays and set neutral defaults: no prescribed
        ! fluxes, no precipitation, sea-level pressure.

        implicit none

        type(chion_forcing_class), intent(INOUT) :: forc
        integer,                   intent(IN)    :: ncol

        call chion_forcing_dealloc(forc)

        forc%ncol = ncol

        allocate(forc%air_temperature(ncol))
        allocate(forc%snowfall_rate(ncol))
        allocate(forc%rainfall_rate(ncol))
        allocate(forc%shortwave_down(ncol))
        allocate(forc%wind_speed(ncol))

        allocate(forc%q_sw_net(ncol))
        allocate(forc%q_lw_down(ncol))
        allocate(forc%q_sh(ncol))
        allocate(forc%q_lh(ncol))

        allocate(forc%has_q_sw_net(ncol))
        allocate(forc%has_q_lw_down(ncol))
        allocate(forc%has_q_sh(ncol))
        allocate(forc%has_q_lh(ncol))

        allocate(forc%relative_humidity(ncol))
        allocate(forc%has_relative_humidity(ncol))

        allocate(forc%surface_height(ncol))
        allocate(forc%air_pressure(ncol))
        allocate(forc%prescribed_albedo(ncol))
        allocate(forc%has_prescribed_albedo(ncol))

        allocate(forc%coszm(ncol))
        allocate(forc%has_coszm(ncol))
        allocate(forc%cloud(ncol))
        allocate(forc%has_cloud(ncol))
        allocate(forc%dust_dep(ncol))
        allocate(forc%has_dust_dep(ncol))
        allocate(forc%z_sur_std(ncol))
        allocate(forc%has_z_sur_std(ncol))
        allocate(forc%alb_ice_host(ncol))
        allocate(forc%has_alb_ice_host(ncol))

        allocate(forc%latitude_deg(ncol))

        allocate(forc%H_ice(ncol))
        allocate(forc%PDDs(ncol))

        forc%air_temperature = 273.15_wp
        forc%snowfall_rate   = 0.0_wp
        forc%rainfall_rate   = 0.0_wp
        forc%shortwave_down  = 0.0_wp
        forc%wind_speed      = 0.0_wp

        forc%q_sw_net  = 0.0_wp
        forc%q_lw_down = 0.0_wp
        forc%q_sh      = 0.0_wp
        forc%q_lh      = 0.0_wp

        forc%has_q_sw_net  = .FALSE.
        forc%has_q_lw_down = .FALSE.
        forc%has_q_sh      = .FALSE.
        forc%has_q_lh      = .FALSE.

        forc%relative_humidity     = 0.0_wp
        forc%has_relative_humidity = .FALSE.

        forc%surface_height        = 0.0_wp
        forc%air_pressure          = DEF_SEA_LEVEL_AIR_PRESSURE
        forc%prescribed_albedo     = 0.0_wp
        forc%has_prescribed_albedo = .FALSE.

        forc%coszm     = 0.0_wp
        forc%has_coszm = .FALSE.
        forc%cloud     = 0.0_wp
        forc%has_cloud = .FALSE.
        forc%dust_dep     = 0.0_wp
        forc%has_dust_dep = .FALSE.
        forc%z_sur_std     = 0.0_wp
        forc%has_z_sur_std = .FALSE.
        forc%alb_ice_host     = 0.0_wp
        forc%has_alb_ice_host = .FALSE.

        forc%latitude_deg = 0.0_wp

        ! ITM-only. H_ice = 0 selects calc_albedo_surface's land branch, and
        ! PDDs = 0 selects the "desert" critical snow depth. Both are neutral
        ! starting points; a host running model="itm" must set them.
        forc%H_ice = 0.0_wp
        forc%PDDs  = 0.0_wp

        forc%day_of_year         = 1.0_wp
        forc%solar_longitude_deg = 0.0_wp

        return

    end subroutine chion_forcing_alloc

    subroutine chion_forcing_dealloc(forc)

        implicit none

        type(chion_forcing_class), intent(INOUT) :: forc

        if (allocated(forc%air_temperature))       deallocate(forc%air_temperature)
        if (allocated(forc%snowfall_rate))         deallocate(forc%snowfall_rate)
        if (allocated(forc%rainfall_rate))         deallocate(forc%rainfall_rate)
        if (allocated(forc%shortwave_down))        deallocate(forc%shortwave_down)
        if (allocated(forc%wind_speed))            deallocate(forc%wind_speed)
        if (allocated(forc%q_sw_net))              deallocate(forc%q_sw_net)
        if (allocated(forc%q_lw_down))             deallocate(forc%q_lw_down)
        if (allocated(forc%q_sh))                  deallocate(forc%q_sh)
        if (allocated(forc%q_lh))                  deallocate(forc%q_lh)
        if (allocated(forc%has_q_sw_net))          deallocate(forc%has_q_sw_net)
        if (allocated(forc%has_q_lw_down))         deallocate(forc%has_q_lw_down)
        if (allocated(forc%has_q_sh))              deallocate(forc%has_q_sh)
        if (allocated(forc%has_q_lh))              deallocate(forc%has_q_lh)
        if (allocated(forc%relative_humidity))     deallocate(forc%relative_humidity)
        if (allocated(forc%has_relative_humidity)) deallocate(forc%has_relative_humidity)
        if (allocated(forc%surface_height))        deallocate(forc%surface_height)
        if (allocated(forc%air_pressure))          deallocate(forc%air_pressure)
        if (allocated(forc%prescribed_albedo))     deallocate(forc%prescribed_albedo)
        if (allocated(forc%has_prescribed_albedo)) deallocate(forc%has_prescribed_albedo)
        if (allocated(forc%coszm))                 deallocate(forc%coszm)
        if (allocated(forc%has_coszm))             deallocate(forc%has_coszm)
        if (allocated(forc%cloud))                 deallocate(forc%cloud)
        if (allocated(forc%has_cloud))             deallocate(forc%has_cloud)
        if (allocated(forc%dust_dep))              deallocate(forc%dust_dep)
        if (allocated(forc%has_dust_dep))          deallocate(forc%has_dust_dep)
        if (allocated(forc%z_sur_std))             deallocate(forc%z_sur_std)
        if (allocated(forc%has_z_sur_std))         deallocate(forc%has_z_sur_std)
        if (allocated(forc%alb_ice_host))          deallocate(forc%alb_ice_host)
        if (allocated(forc%has_alb_ice_host))      deallocate(forc%has_alb_ice_host)
        if (allocated(forc%latitude_deg))          deallocate(forc%latitude_deg)
        if (allocated(forc%H_ice))                 deallocate(forc%H_ice)
        if (allocated(forc%PDDs))                  deallocate(forc%PDDs)

        forc%ncol = 0

        return

    end subroutine chion_forcing_dealloc

    subroutine chion_grid_init(grd,ncol,x,y,js,is,mask)
        ! Initialize a column list. Spatial coordinates are optional and are
        ! used only for NetCDF output; supply either all of x/y/js/is or none.

        implicit none

        type(chion_grid_class), intent(INOUT) :: grd
        integer,                intent(IN)    :: ncol
        real(wp), optional,     intent(IN)    :: x(:)
        real(wp), optional,     intent(IN)    :: y(:)
        integer,  optional,     intent(IN)    :: js(:)
        integer,  optional,     intent(IN)    :: is(:)
        real(wp), optional,     intent(IN)    :: mask(:,:)

        ! Local variables
        integer :: n_present

        call chion_grid_dealloc(grd)

        if (ncol .le. 0) then
            write(io_unit_err,*) "chion_grid_init:: Error: ncol must be positive."
            write(io_unit_err,*) "ncol = ", ncol
            stop "Program stopped."
        end if

        grd%ncol = ncol

        allocate(grd%active(ncol))
        allocate(grd%active_idx(ncol))

        grd%active = .TRUE.
        call chion_grid_set_active(grd,grd%active)

        n_present = 0
        if (present(x))  n_present = n_present + 1
        if (present(y))  n_present = n_present + 1
        if (present(js)) n_present = n_present + 1
        if (present(is)) n_present = n_present + 1

        if (n_present .eq. 0) then
            grd%has_spatial = .FALSE.
            return
        else if (n_present .lt. 4) then
            write(io_unit_err,*) "chion_grid_init:: Error: provide all of x, y, js, is, or none."
            stop "Program stopped."
        end if

        if (size(js) .ne. ncol .or. size(is) .ne. ncol) then
            write(io_unit_err,*) "chion_grid_init:: Error: js and is must have length ncol."
            write(io_unit_err,*) "ncol, size(js), size(is) = ", ncol, size(js), size(is)
            stop "Program stopped."
        end if

        grd%has_spatial = .TRUE.

        allocate(grd%x(size(x)))
        allocate(grd%y(size(y)))
        allocate(grd%js(ncol))
        allocate(grd%is(ncol))
        allocate(grd%mask(size(y),size(x)))

        grd%x  = x
        grd%y  = y
        grd%js = js
        grd%is = is

        if (present(mask)) then
            if (size(mask,1) .ne. size(y) .or. size(mask,2) .ne. size(x)) then
                write(io_unit_err,*) "chion_grid_init:: Error: mask must have shape (ny,nx)."
                write(io_unit_err,*) "shape(mask), ny, nx = ", shape(mask), size(y), size(x)
                stop "Program stopped."
            end if
            grd%mask = mask
        else
            grd%mask = 1.0_wp
        end if

        return

    end subroutine chion_grid_init

    subroutine chion_grid_dealloc(grd)

        implicit none

        type(chion_grid_class), intent(INOUT) :: grd

        if (allocated(grd%x))          deallocate(grd%x)
        if (allocated(grd%y))          deallocate(grd%y)
        if (allocated(grd%js))         deallocate(grd%js)
        if (allocated(grd%is))         deallocate(grd%is)
        if (allocated(grd%mask))       deallocate(grd%mask)
        if (allocated(grd%active))     deallocate(grd%active)
        if (allocated(grd%active_idx)) deallocate(grd%active_idx)

        grd%ncol        = 0
        grd%n_active    = 0
        grd%has_spatial = .FALSE.

        return

    end subroutine chion_grid_dealloc

    subroutine chion_grid_set_active(grd,active)
        ! Set the active-column mask and rebuild the packed active index list.
        ! Mirrors Chion.jl set_active_mask! (src/integrators.jl). Resetting the
        ! state of newly-deactivated columns is the caller's responsibility
        ! (see WP11), because it needs the model state.

        implicit none

        type(chion_grid_class), intent(INOUT) :: grd
        logical,                intent(IN)    :: active(:)

        ! Local variables
        integer :: i, n

        if (size(active) .ne. grd%ncol) then
            write(io_unit_err,*) "chion_grid_set_active:: Error: mask length must equal ncol."
            write(io_unit_err,*) "ncol, size(active) = ", grd%ncol, size(active)
            stop "Program stopped."
        end if

        grd%active = active

        n = 0
        do i = 1, grd%ncol
            if (grd%active(i)) then
                n = n + 1
                grd%active_idx(n) = i
            end if
        end do

        grd%n_active = n

        return

    end subroutine chion_grid_set_active

    function chion_semix_snow_albedo_flag(name) result(flag)
        ! Map a namelist string onto a SEMIX spectral snow-albedo flag.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("ww","warren","warren_wiscombe")
                flag = SEMIX_SNOW_ALBEDO_WW
            case("dang")
                flag = SEMIX_SNOW_ALBEDO_DANG
            case DEFAULT
                write(io_unit_err,*) "chion_semix_snow_albedo_flag:: Error: scheme not recognized."
                write(io_unit_err,*) "semix_snow_albedo should be one of: ['ww','dang'] &
                                     &(aliases: 'warren','warren_wiscombe' -> 'ww')"
                write(io_unit_err,*) "semix_snow_albedo = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_semix_snow_albedo_flag

    function chion_seb_scheme_flag(name) result(flag)
        ! Map a namelist string onto a surface-energy-balance scheme flag.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("bessi")
                flag = CHION_SEB_BESSI
            case("semix")
                flag = CHION_SEB_SEMIX
            case DEFAULT
                write(io_unit_err,*) "chion_seb_scheme_flag:: Error: seb scheme not recognized."
                write(io_unit_err,*) "seb_scheme should be one of: ['bessi','semix']"
                write(io_unit_err,*) "seb_scheme = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_seb_scheme_flag

    function chion_semix_qsat_flag(name) result(flag)
        ! Map a namelist string onto a saturation-humidity parameterization.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("semix","climberx")
                flag = SEMIX_QSAT_SEMIX
            case("bessi","chion")
                flag = SEMIX_QSAT_BESSI
            case DEFAULT
                write(io_unit_err,*) "chion_semix_qsat_flag:: Error: scheme not recognized."
                write(io_unit_err,*) "semix_qsat should be one of: ['semix','bessi'] &
                                     &(aliases: 'climberx' -> 'semix', 'chion' -> 'bessi')"
                write(io_unit_err,*) "semix_qsat = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_semix_qsat_flag

    function chion_albedo_scheme_flag(name) result(flag)
        ! Map a namelist string onto an albedo scheme flag.
        ! Chion.jl aliases :bessi and :legacy to :constant
        ! (src/constants.jl:151-165); those aliases are preserved.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("constant","bessi","legacy")
                flag = CHION_ALBEDO_CONSTANT
            case("dynamic")
                flag = CHION_ALBEDO_DYNAMIC
            case("prescribed")
                flag = CHION_ALBEDO_PRESCRIBED
            case("semix")
                flag = CHION_ALBEDO_SEMIX
            case DEFAULT
                write(io_unit_err,*) "chion_albedo_scheme_flag:: Error: albedo scheme not recognized."
                write(io_unit_err,*) "albedo_scheme should be one of: &
                                     &['constant','dynamic','prescribed','semix'] &
                                     &(aliases: 'bessi','legacy' -> 'constant')"
                write(io_unit_err,*) "albedo_scheme = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_albedo_scheme_flag

    function chion_fresh_snow_density_scheme_flag(name) result(flag)
        ! Chion.jl aliases :bessi -> :constant and :htessel -> :parameterized
        ! (src/constants.jl:131-147).

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("constant","bessi")
                flag = CHION_FRESH_SNOW_DENSITY_CONSTANT
            case("parameterized","htessel")
                flag = CHION_FRESH_SNOW_DENSITY_PARAMETERIZED
            case DEFAULT
                write(io_unit_err,*) "chion_fresh_snow_density_scheme_flag:: Error: &
                                     &fresh snow density scheme not recognized."
                write(io_unit_err,*) "fresh_snow_density_scheme should be one of: &
                                     &['constant','parameterized'] &
                                     &(aliases: 'bessi','htessel')"
                write(io_unit_err,*) "fresh_snow_density_scheme = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_fresh_snow_density_scheme_flag

    function chion_densify_scheme_flag(name) result(flag)
        ! NOTE: Chion.jl dispatches this with an if/else on HTESSEL only, so
        ! any unrecognized value silently falls through to BESSI
        ! (src/processes/densification.jl:262). chion validates instead.
        ! See docs/PLAN.md section 5, item 8.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("bessi")
                flag = CHION_DENSIFY_BESSI
            case("htessel")
                flag = CHION_DENSIFY_HTESSEL
            case DEFAULT
                write(io_unit_err,*) "chion_densify_scheme_flag:: Error: &
                                     &densification scheme not recognized."
                write(io_unit_err,*) "low_density_densification should be one of: ['bessi','htessel']"
                write(io_unit_err,*) "low_density_densification = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_densify_scheme_flag

    subroutine chion_check_enum(group,varname,value,allowed)
        ! Validate a character parameter against a '|'-delimited list of
        ! allowed values. Follows yelmo_check_enum (yelmo/src/yelmo_defs.f90:1130).

        implicit none

        character(len=*), intent(IN) :: group
        character(len=*), intent(IN) :: varname
        character(len=*), intent(IN) :: value
        character(len=*), intent(IN) :: allowed

        ! Local variables
        integer :: i0, i1
        logical :: found

        found = .FALSE.
        i0    = 1

        do while (i0 .le. len_trim(allowed))
            i1 = index(allowed(i0:),"|")
            if (i1 .eq. 0) then
                i1 = len_trim(allowed) + 1
            else
                i1 = i0 + i1 - 1
            end if
            if (trim(adjustl(allowed(i0:i1-1))) .eq. trim(adjustl(value))) then
                found = .TRUE.
                exit
            end if
            i0 = i1 + 1
        end do

        if (.not. found) then
            write(io_unit_err,*) "chion_check_enum:: Error: parameter value not recognized."
            write(io_unit_err,*) "group   = ", trim(group)
            write(io_unit_err,*) "name    = ", trim(varname)
            write(io_unit_err,*) "value   = ", trim(value)
            write(io_unit_err,*) "allowed = ", trim(allowed)
            stop "Program stopped."
        end if

        return

    end subroutine chion_check_enum

    subroutine chion_check_file(filename)
        ! Stop with a clear message if a required input file is missing.
        ! Follows yelmo_check_file (yelmo/src/yelmo_defs.f90:1174).

        implicit none

        character(len=*), intent(IN) :: filename

        ! Local variables
        logical :: exists

        inquire(file=trim(filename),exist=exists)

        if (.not. exists) then
            write(io_unit_err,*) "chion_check_file:: Error: file does not exist."
            write(io_unit_err,*) "filename = ", trim(filename)
            stop "Program stopped."
        end if

        return

    end subroutine chion_check_file

    subroutine chion_parse_path(path,domain,grid_name,rundir)
        ! Expand {domain}, {grid_name} and {rundir} placeholders in a path
        ! string. Follows yelmo_parse_path (yelmo/src/yelmo_defs.f90:1228).

        implicit none

        character(len=*),           intent(INOUT) :: path
        character(len=*), optional, intent(IN)    :: domain
        character(len=*), optional, intent(IN)    :: grid_name
        character(len=*), optional, intent(IN)    :: rundir

        if (present(domain))    call nml_replace(path,"{domain}",   trim(domain))
        if (present(grid_name)) call nml_replace(path,"{grid_name}",trim(grid_name))
        if (present(rundir))    call nml_replace(path,"{rundir}",   trim(rundir))

        return

    end subroutine chion_parse_path

    subroutine chion_load_command_line_args(path_par)
        ! Read the parameter file path from argv(1). Exactly one argument is
        ! required; this is the contract runme depends on
        ! (runme/.runme/info.json: par_path_as_argument).
        ! Follows yelmo_load_command_line_args (yelmo/src/yelmo_defs.f90:1267).

        implicit none

        character(len=*), intent(OUT) :: path_par

        ! Local variables
        integer :: narg

        narg = command_argument_count()

        if (narg .ne. 1) then
            write(io_unit_err,*) "chion_load_command_line_args:: Error: &
                                 &exactly one argument is required, the parameter file path."
            write(io_unit_err,*) "n arguments = ", narg
            write(io_unit_err,*) "usage: chion_<program>.x path/to/par_file.nml"
            stop "Program stopped."
        end if

        call get_command_argument(1,path_par)

        return

    end subroutine chion_load_command_line_args

end module chion_defs
