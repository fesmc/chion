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
    use phys_constants, only : sec_year_360d, sec_day

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
    integer, parameter, public :: CHION_ALBEDO_AGING      = 5   ! Chion.jl ALBEDO_AGING = 4

    ! Which spectral snow-albedo parameterization the SEMIX scheme uses.
    ! CLIMBER-X defaults to Dang (its isnow_albedo = 2).
    integer, parameter, public :: SEMIX_SNOW_ALBEDO_WW   = 1
    integer, parameter, public :: SEMIX_SNOW_ALBEDO_DANG = 2

    ! Surface energy balance: since Chion.jl d0146e1 it selects the LONGWAVE
    ! treatment only. "bessi": LWdn - eps_snow*sigma*Ts^4 (the downwelling
    ! flux absorbed in full); "semix": eps_s*(LWdn - sigma*Ts^4), with eps_s
    ! = eps_snow on snow and eps_ice on bare ice (CLIMBER-X SEMIX's ebal).
    integer, parameter, public :: CHION_SEB_BESSI = 1
    integer, parameter, public :: CHION_SEB_SEMIX = 2

    ! Turbulent sensible and latent heat, independent of seb_scheme (Chion.jl
    ! turbulent_flux_scheme). "bessi": the bulk coefficient D_sh and BESSI's
    ! vapour-pressure gradient; "semix": Chion.jl's bulk scheme of that name
    ! (log-law, Richardson stability, calibrated exchange factors;
    ! snow_turbulence); "climberx": CLIMBER-X SEMIX's aerodynamic resistance
    ! (snow_seb_semix; docs/semix_port_scope.md). Orthogonal to albedo_scheme
    ! and to Ntot.
    integer, parameter, public :: CHION_TURB_BESSI    = 1
    integer, parameter, public :: CHION_TURB_SEMIX    = 2
    integer, parameter, public :: CHION_TURB_CLIMBERX = 3

    ! Downwelling longwave when the host does not prescribe it (Chion.jl
    ! 03bb445 longwave_scheme): "graybody" eps_air*sigma*T_air^4, or
    ! "cloud_proxy", an emissivity from air temperature and a shortwave
    ! cloudiness proxy, resolved once per forcing step.
    integer, parameter, public :: CHION_LONGWAVE_GRAYBODY    = 1
    integer, parameter, public :: CHION_LONGWAVE_CLOUD_PROXY = 2

    ! Which saturation-specific-humidity parameterization the CLIMBER-X
    ! turbulence uses for the latent flux. "climberx" is CLIMBER-X's own
    ! q_sat_i/dqsat_dT_i; "bessi" routes chion's ice saturation vapour
    ! pressure through the same 0.622/p conversion.
    integer, parameter, public :: CLIMBERX_QSAT_CLIMBERX = 1
    integer, parameter, public :: CLIMBERX_QSAT_BESSI    = 2

    integer, parameter, public :: CHION_DENSIFY_BESSI   = 1
    integer, parameter, public :: CHION_DENSIFY_HTESSEL = 2

    ! === Defaults ============================================================
    ! Chion.jl/src/constants.jl:16-40.

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
    !   * DENSIFY_R_GAS  -> 8.314, Chion.jl's rounded gas constant (it had
    !     the typo 8.13 until dev_nils 408e91c; Chion.jl issue #18,
    !     docs/porting_notes.md D22). Drop once Chion.jl uses 8.31446...
    !   * DENSIFY_GRAVITY -> 9.81, Chion.jl's second gravity constant
    !     (docs/porting_notes.md D25).
    !   * ITM_FIRN_DAYS_YEAR -> 1, i.e. ITM's tsrf applies firn_fac to the
    !     daily melt_net rate, as Chion.jl's ITMModel does (docs/porting_notes.md
    !     D27). chion scales it to the annual rate firn_fac is calibrated on.
    !   * ALBEDO_AGING_BINARY_REFRESH -> .TRUE., i.e. the aging albedo scheme
    !     resets to alpha_dry (age 0) on any snowfall, as Chion.jl does.
    !     chion rejuvenates in proportion to the step's snowfall
    !     (docs/porting_notes.md D30).
    !   * NEAR_SURFACE_SPLIT_MERGE_BELOW -> .FALSE., i.e. with fine
    !     near-surface layers everything below them stays in one layer that
    !     is never split or merged, as in Chion.jl. chion splits and merges
    !     it by mass like the surface layer (docs/porting_notes.md D32).
    !   * TURB_SEMIX_ICE_SUBLIMATION -> .FALSE., i.e. the semix turbulence's
    !     latent exchange over bare ice carries the phase's latent heat at
    !     the surface temperature (Lv at T0), as in Chion.jl; chion uses
    !     Lv + Lm, the latent heat its vapour mass is converted with
    !     (docs/porting_notes.md D35).
    !   * TURB_SEMIX_R_AIR_LITERAL -> .TRUE., i.e. the semix turbulence's air
    !     density uses Chion.jl's literal 287.05 instead of c%R_dry
    !     (docs/porting_notes.md D38).
    !
    ! Not covered: the PDD budget (D23). Chion.jl adopted it (ce6a68d), so
    ! the plain build is gated against Chion.jl for PDD.
    !
    ! ITM_FIRN_DAYS_YEAR: days per year of the calendar ITM's firn_fac is
    ! calibrated on. smbpal's annual totals are on a 360-day year; a property
    ! of the calibration, not of the host's calendar.
#ifdef CHION_LEGACY
    real(wp_acc), parameter, public :: DENSIFY_R_GAS   = 8.314_wp_acc
    real(wp_acc), parameter, public :: DENSIFY_GRAVITY = 9.81_wp_acc
    real(wp),     parameter, public :: ITM_FIRN_DAYS_YEAR = 1.0_wp
    logical,      parameter, public :: ALBEDO_AGING_BINARY_REFRESH = .TRUE.
    logical,      parameter, public :: NEAR_SURFACE_SPLIT_MERGE_BELOW = .FALSE.
    logical,      parameter, public :: TURB_SEMIX_ICE_SUBLIMATION = .FALSE.
    logical,      parameter, public :: TURB_SEMIX_R_AIR_LITERAL = .TRUE.
    logical,      parameter, public :: CHION_LEGACY_MODE = .TRUE.
#else
    real(wp_acc), parameter, public :: DENSIFY_R_GAS   = real(DEF_UNIVERSAL_GAS_CONSTANT,wp_acc)
    real(wp_acc), parameter, public :: DENSIFY_GRAVITY = real(DEF_GRAVITY,wp_acc)
    real(wp),     parameter, public :: ITM_FIRN_DAYS_YEAR = real(sec_year_360d/sec_day,wp)
    logical,      parameter, public :: ALBEDO_AGING_BINARY_REFRESH = .FALSE.
    logical,      parameter, public :: NEAR_SURFACE_SPLIT_MERGE_BELOW = .TRUE.
    logical,      parameter, public :: TURB_SEMIX_ICE_SUBLIMATION = .TRUE.
    logical,      parameter, public :: TURB_SEMIX_R_AIR_LITERAL = .FALSE.
    logical,      parameter, public :: CHION_LEGACY_MODE = .FALSE.
#endif

    integer,  parameter, public :: DEF_NTOT             = 15
    real(wp), parameter, public :: DEF_MASS_MAX         = 500.0_wp
    real(wp), parameter, public :: DEF_MASS_SPLIT       = 300.0_wp
    real(wp), parameter, public :: DEF_MASS_MIN         = 100.0_wp
    real(wp), parameter, public :: DEF_DENSITY_INIT     = 300.0_wp
    real(wp), parameter, public :: DEF_TEMPERATURE_INIT = 273.0_wp

    ! Thermal ice substrate below the snow/firn column (Chion.jl 03bb445,
    ! src/models.jl:98-99). Upstream defaults to 5 layers; chion keeps 0
    ! (no substrate, the bare-ice-at-T0 treatment) until the default switch.
    integer,  parameter, public :: DEF_ICE_SUBSTRATE_LAYERS        = 0
    real(wp), parameter, public :: DEF_ICE_SUBSTRATE_TOP_THICKNESS = 0.05_wp

    ! Fine near-surface layers (Chion.jl 03bb445, src/models.jl:97): maximum
    ! thicknesses of the top NEAR_SURFACE_LAYERS layers, held by a
    ! conservative remesh. 0 = no limit (Julia Inf). Upstream defaults to
    ! (0.02, 0.05, 0.10, 0.30) m; chion keeps them off until the default
    ! switch (docs/porting_notes.md D36).
    integer,  parameter, public :: NEAR_SURFACE_LAYERS = 4
    real(wp), parameter, public :: DEF_NEAR_SURFACE_LAYER_MAX_THICKNESSES(NEAR_SURFACE_LAYERS) = 0.0_wp

    ! Depth cap: a fixed total solid depth, independent of Ntot and of
    ! mass_split (Chion.jl 03bb445, src/constants.jl:39). It replaces the
    ! former 15*mass_split*1.5/300, identical at mass_split = 300.
    real(wp), parameter, public :: BESSI_REFERENCE_SNOW_DEPTH_M = 22.5_wp

    ! === Physical constants ==================================================

    type chion_const_class
        ! Mirrors Chion.jl SnowpackPhysicalConstants (src/constants.jl:61-88),
        ! field for field and in the same order. Defaults at src/constants.jl:173-200.
        !
        ! Porting note: Chion.jl keeps the three scheme flags inside this struct.
        ! chion does the same, so that every physics routine (which already
        ! receives c) can branch without threading an extra argument.
        !
        ! SHARED fields -- rho_i, rho_w, ci, cw, Lm, grav, T0 -- are not chion's
        ! to choose: chion_const_from_phys fills them from the program's
        ! fesm-utils phys_const_class (rho_ice, rho_w, cp_ice, cp_w, L_ice, g,
        ! T0), so a host hands chion the same constants as every other
        ! component. All other fields are chion's own (&chion_const in
        ! input/chion_defaults.nml). Day length is the named convention
        ! phys_constants:sec_day, not a field. See docs/porting_notes.md D28.

        ! Densities
        real(wp) :: rho_s              ! [kg m-3] fresh snow density (constant scheme)
        real(wp) :: rho_i              ! [kg m-3] ice density              (shared: rho_ice)
        real(wp) :: rho_w              ! [kg m-3] water density            (shared: rho_w)

        ! Parameterized fresh-snow density: rho = a + b*(T-T0) + c*sqrt(wind)
        real(wp) :: rho_s_a            ! [kg m-3]
        real(wp) :: rho_s_b            ! [kg m-3 K-1]
        real(wp) :: rho_s_c            ! [kg m-3 (m s-1)^-1/2]
        integer  :: fresh_snow_density_scheme   ! CHION_FRESH_SNOW_DENSITY_*

        ! Thermal properties
        real(wp) :: ci                 ! [J kg-1 K-1] heat capacity of ice   (shared: cp_ice)
        real(wp) :: cw                 ! [J kg-1 K-1] heat capacity of water (shared: cp_w)
        real(wp) :: Lm                 ! [J kg-1] latent heat of melting     (shared: L_ice)
        real(wp) :: Lv                 ! [J kg-1] latent heat of vaporization
        real(wp) :: cp_air             ! [J kg-1 K-1] heat capacity of air
        real(wp) :: latent_heat_flux_ratio  ! [1] scaling of turbulent latent flux

        ! Turbulent exchange
        real(wp) :: D_sh               ! [W m-2 K-1] sensible heat exchange coefficient

        ! Surface energy balance (longwave) and turbulent-flux schemes.
        integer  :: seb_scheme             ! CHION_SEB_*
        integer  :: turbulent_flux_scheme  ! CHION_TURB_*

        ! CLIMBER-X SEMIX aerodynamic exchange (CHION_TURB_CLIMBERX only).
        ! Roughness lengths and the surface-layer height are CLIMBER-X smb_par
        ! / constants values; karman, grav and R_dry are the universal
        ! constants SEMIX pulls from its constants module. The heat capacity
        ! of air is chion's existing cp_air (1003 vs SEMIX's 1000 J kg-1 K-1,
        ! 0.3% on f_sh) rather than a second constant for the same quantity.
        real(wp) :: z0m_snow           ! [m] momentum roughness length, snow
        real(wp) :: z0m_ice            ! [m] momentum roughness length, ice
        real(wp) :: zm_to_zh           ! [1] heat/momentum roughness ratio
        real(wp) :: z_sfl              ! [m] surface layer height
        real(wp) :: karman             ! [1] von Karman constant
        real(wp) :: grav               ! [m s-2] gravitational acceleration (shared: g)
        real(wp) :: R_dry              ! [J kg-1 K-1] gas constant of dry air
        logical  :: l_neutral          ! [1] force neutral stratification
        logical  :: l_dew              ! [1] allow dew/frost deposition
        integer  :: climberx_qsat      ! CLIMBERX_QSAT_*

        ! Chion.jl's bulk turbulence (CHION_TURB_SEMIX only; Chion.jl's
        ! semix_* keywords, 03bb445 defaults). z0h = z0m/semix_zm_to_zh; the
        ! sensible exchange factor and the stable coefficient b of the
        ! Richardson damping 1/(1 + b Ri) are calibrated against MAR.
        real(wp) :: semix_karman                    ! [1] von Karman constant
        real(wp) :: semix_surface_height            ! [m] reference height of T_a, wind
        real(wp) :: semix_z0m_snow                  ! [m] momentum roughness, snow
        real(wp) :: semix_z0m_ice                   ! [m] momentum roughness, bare ice
        real(wp) :: semix_zm_to_zh                  ! [1] z0m/z0h
        real(wp) :: semix_sensible_exchange_factor  ! [1]
        real(wp) :: semix_stable_coefficient        ! [1] b
        real(wp) :: semix_latent_exchange_factor    ! [1]

        ! Albedo
        real(wp) :: alpha_dry          ! [1] dry snow albedo (upper bound)
        real(wp) :: alpha_wet          ! [1] wet snow albedo (lower bound)
        real(wp) :: alpha_ice          ! [1] bare ice albedo
        real(wp) :: max_lwc_albedo     ! [1] LWC at which albedo reaches alpha_wet
        integer  :: albedo_scheme      ! CHION_ALBEDO_*

        ! Snowfall-age albedo (CHION_ALBEDO_AGING): e-folding time of the
        ! relaxation towards alpha_wet, cold vs melting surface.
        real(wp) :: aging_cold_timescale_days     ! [d]
        real(wp) :: aging_melting_timescale_days  ! [d]
        ! Step snowfall that fully rejuvenates the aging albedo (D30).
        real(wp) :: aging_snowfall_ref            ! [kg m-2]

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
        real(wp) :: eps_air            ! [1] emissivity of air (graybody longwave)

        ! Cloud-proxy downwelling longwave (CHION_LONGWAVE_CLOUD_PROXY,
        ! Chion.jl 03bb445): eps = base + temperature_slope*(T_a - T0)
        ! + cloud_slope*n, n = 1 - SWdn/(TOA*tau_clear(z)), tau_clear =
        ! clear_sky_transmissivity + clear_sky_transmissivity_per_km*z/1000;
        ! n = night_cloud_fraction without a daily-mean TOA.
        integer  :: longwave_scheme                     ! CHION_LONGWAVE_*
        real(wp) :: lw_emissivity_base                  ! [1]
        real(wp) :: lw_emissivity_temperature_slope     ! [K-1]
        real(wp) :: lw_emissivity_cloud_slope           ! [1]
        real(wp) :: lw_clear_sky_transmissivity         ! [1] at sea level
        real(wp) :: lw_clear_sky_transmissivity_per_km  ! [km-1]
        real(wp) :: lw_night_cloud_fraction             ! [1]

        real(wp) :: eps_snow           ! [1] emissivity of snow
        real(wp) :: eps_ice            ! [1] emissivity of bare ice (semix SEB)
        real(wp) :: sigma_sb           ! [W m-2 K-4] Stefan-Boltzmann constant

        ! Reference values
        real(wp) :: T0                 ! [K] freezing point of water       (shared: T0)

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
        real(wp) :: surface_height = 0.0_wp  ! [m] diurnal T amplitude gradient; non-finite = no excess
        real(wp) :: day_of_year         ! [d] fractional, 1-based
        real(wp) :: solar_longitude_deg ! [deg]

        ! chion only, not in SnowpackStepForcing: the host's ice thickness.
        ! BESSI puts its thermal ice substrate only under ice (H_ice > 0); a
        ! land column (H_ice = 0, the forcing default) has none
        ! (docs/porting_notes.md D34). ITM still takes it as an argument.
        real(wp) :: H_ice = 0.0_wp      ! [m]

        ! chion only: the host's daily-mean top-of-atmosphere shortwave. When
        ! given, the cloud-proxy longwave divides by it instead of chion's
        ! fixed-orbit TOA from latitude and season (docs/porting_notes.md D33).
        real(wp) :: toa_shortwave = 0.0_wp       ! [W m-2]
        logical  :: has_toa_shortwave = .FALSE.
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

        real(wp), allocatable :: surface_height(:)       ! [m] host air pressure; BESSI diurnal amplitude; ITM
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

        ! Optional daily-mean top-of-atmosphere shortwave for the cloud-proxy
        ! longwave (chion only, D33); unset, chion computes a fixed-orbit TOA.
        real(wp), allocatable :: toa_shortwave(:)        ! [W m-2]
        logical,  allocatable :: has_toa_shortwave(:)

        ! --- Ice-sheet fields (WP11; H_ice also BESSI since C3) ---------
        !
        ! PDDs is deliberately NOT part of chion_step_forcing_class. That type
        ! mirrors Chion.jl's SnowpackStepForcing and is the shared,
        ! model-neutral contract every kernel takes. ITM receives z_srf, H_ice
        ! and PDDs as explicit arguments from the dispatcher
        ! (itm_step(itm,icol,fc,z_srf,H_ice,PDDs)). H_ice is packed into the
        ! step forcing as well, because BESSI places its thermal ice substrate
        ! only where H_ice > 0 (docs/porting_notes.md D34).
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
        character(len=56)  :: nml_const      ! group holding chion's own constants

        ! The shared constants (fesm-utils phys_const_class), read only when the
        ! host passes none to chion_init.
        character(len=512) :: phys_const_file ! path to the phys_const file
        character(len=512) :: phys_const     ! group name within phys_const_file
        character(len=512) :: phys_const_src ! provenance of the shared constants used
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
    public :: chion_const_validate

    public :: chion_forcing_alloc
    public :: chion_forcing_dealloc
    public :: chion_grid_init
    public :: chion_grid_dealloc
    public :: chion_grid_set_active

    public :: chion_albedo_scheme_flag
    public :: chion_semix_snow_albedo_flag
    public :: chion_seb_scheme_flag
    public :: chion_longwave_scheme_flag
    public :: chion_turbulent_flux_scheme_flag
    public :: chion_climberx_qsat_flag
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

        c%ci      = 2110.0_wp
        c%cw      = 4181.0_wp
        c%Lm      = 334000.0_wp
        c%Lv      = 2.501e6_wp
        c%cp_air  = 1003.0_wp
        c%latent_heat_flux_ratio = 1.0_wp

        c%D_sh    = 10.0_wp

        c%seb_scheme            = CHION_SEB_BESSI
        c%turbulent_flux_scheme = CHION_TURB_BESSI

        ! CLIMBER-X SEMIX aerodynamic exchange defaults (CLIMBER-X
        ! smb_par.nml / smb_params.f90 / constants.f90).
        c%z0m_snow    = 0.0024_wp
        c%z0m_ice     = 0.002_wp
        c%zm_to_zh    = exp(-2.0_wp)
        c%z_sfl       = 100.0_wp
        c%karman      = 0.4_wp
        c%grav        = 9.81_wp
        c%R_dry       = 287.058_wp
        c%l_neutral   = .FALSE.
        c%l_dew       = .TRUE.
        c%climberx_qsat = CLIMBERX_QSAT_CLIMBERX

        ! Chion.jl bulk turbulence defaults (src/constants.jl, 03bb445).
        c%semix_karman                   = 0.4_wp
        c%semix_surface_height           = 10.0_wp
        c%semix_z0m_snow                 = 0.001_wp
        c%semix_z0m_ice                  = 0.01_wp
        c%semix_zm_to_zh                 = 10.0_wp
        c%semix_sensible_exchange_factor = 2.5_wp
        c%semix_stable_coefficient       = 40.0_wp
        c%semix_latent_exchange_factor   = 1.0_wp

        c%alpha_dry      = 0.81_wp
        c%alpha_wet      = 0.70_wp
        c%alpha_ice      = 0.30_wp
        c%max_lwc_albedo = 0.10_wp
        c%albedo_scheme  = CHION_ALBEDO_DYNAMIC

        ! dev_nils defaults (27113b6); 12407a3 had 5 d for the melting surface.
        c%aging_cold_timescale_days    = 20.0_wp
        c%aging_melting_timescale_days =  2.0_wp
        c%aging_snowfall_ref           = 10.0_wp

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

        ! Chion.jl 03bb445 coefficients (fitted to daily MAR longwave over
        ! Greenland). Upstream defaults to the cloud proxy; chion keeps the
        ! graybody until the Stage C default switch.
        c%longwave_scheme                    = CHION_LONGWAVE_GRAYBODY
        c%lw_emissivity_base                 = 0.624_wp
        c%lw_emissivity_temperature_slope    = 0.0032_wp
        c%lw_emissivity_cloud_slope          = 0.613_wp
        c%lw_clear_sky_transmissivity        = 0.85_wp
        c%lw_clear_sky_transmissivity_per_km = 0.075_wp
        c%lw_night_cloud_fraction            = 0.389_wp

        c%eps_snow = 0.98_wp
        c%eps_ice  = 0.98_wp
        c%sigma_sb = 5.670373e-8_wp

        c%T0      = 273.15_wp

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
        write(*,"(a25,g14.6,a)") "ci      = ", c%ci,      "  [J kg-1 K-1]"
        write(*,"(a25,g14.6,a)") "cw      = ", c%cw,      "  [J kg-1 K-1]"
        write(*,"(a25,g14.6,a)") "Lm      = ", c%Lm,      "  [J kg-1]"
        write(*,"(a25,g14.6,a)") "Lv      = ", c%Lv,      "  [J kg-1]"
        write(*,"(a25,g14.6,a)") "cp_air  = ", c%cp_air,  "  [J kg-1 K-1]"
        write(*,"(a25,g14.6,a)") "latent_heat_flux_ratio = ", c%latent_heat_flux_ratio, "  [1]"
        write(*,"(a25,g14.6,a)") "D_sh    = ", c%D_sh,    "  [W m-2 K-1]"
        write(*,"(a25,i14)")     "seb_scheme = ", c%seb_scheme
        write(*,"(a25,i14)")     "turbulent_flux_scheme = ", c%turbulent_flux_scheme
        write(*,"(a25,g14.6,a)") "z0m_snow = ", c%z0m_snow, "  [m]"
        write(*,"(a25,g14.6,a)") "z0m_ice  = ", c%z0m_ice,  "  [m]"
        write(*,"(a25,g14.6,a)") "zm_to_zh = ", c%zm_to_zh, "  [1]"
        write(*,"(a25,g14.6,a)") "z_sfl    = ", c%z_sfl,    "  [m]"
        write(*,"(a25,l14)")     "l_neutral = ", c%l_neutral
        write(*,"(a25,l14)")     "l_dew     = ", c%l_dew
        write(*,"(a25,i14)")     "climberx_qsat = ", c%climberx_qsat
        write(*,"(a37,g14.6,a)") "semix_karman = ", c%semix_karman, "  [1]"
        write(*,"(a37,g14.6,a)") "semix_surface_height = ", c%semix_surface_height, "  [m]"
        write(*,"(a37,g14.6,a)") "semix_z0m_snow = ", c%semix_z0m_snow, "  [m]"
        write(*,"(a37,g14.6,a)") "semix_z0m_ice = ", c%semix_z0m_ice, "  [m]"
        write(*,"(a37,g14.6,a)") "semix_zm_to_zh = ", c%semix_zm_to_zh, "  [1]"
        write(*,"(a37,g14.6,a)") "semix_sensible_exchange_factor = ", &
                                 c%semix_sensible_exchange_factor, "  [1]"
        write(*,"(a37,g14.6,a)") "semix_stable_coefficient = ", c%semix_stable_coefficient, "  [1]"
        write(*,"(a37,g14.6,a)") "semix_latent_exchange_factor = ", &
                                 c%semix_latent_exchange_factor, "  [1]"
        write(*,"(a25,g14.6,a)") "alpha_dry = ", c%alpha_dry, "  [1]"
        write(*,"(a25,g14.6,a)") "alpha_wet = ", c%alpha_wet, "  [1]"
        write(*,"(a25,g14.6,a)") "alpha_ice = ", c%alpha_ice, "  [1]"
        write(*,"(a25,g14.6,a)") "max_lwc_albedo = ", c%max_lwc_albedo, "  [1]"
        write(*,"(a25,i14)")     "albedo_scheme = ", c%albedo_scheme
        write(*,"(a25,g14.6,a)") "aging_cold_timescale_days = ",    c%aging_cold_timescale_days,    "  [d]"
        write(*,"(a25,g14.6,a)") "aging_melting_timescale_days = ", c%aging_melting_timescale_days, "  [d]"
        write(*,"(a25,g14.6,a)") "aging_snowfall_ref = ", c%aging_snowfall_ref, "  [kg m-2]"
        write(*,"(a25,g14.6,a)") "eps_air  = ", c%eps_air,  "  [1]"
        write(*,"(a25,i14)")     "longwave_scheme = ", c%longwave_scheme
        write(*,"(a37,g14.6,a)") "lw_emissivity_base = ", c%lw_emissivity_base, "  [1]"
        write(*,"(a37,g14.6,a)") "lw_emissivity_temperature_slope = ", &
                                 c%lw_emissivity_temperature_slope, "  [K-1]"
        write(*,"(a37,g14.6,a)") "lw_emissivity_cloud_slope = ", c%lw_emissivity_cloud_slope, "  [1]"
        write(*,"(a37,g14.6,a)") "lw_clear_sky_transmissivity = ", &
                                 c%lw_clear_sky_transmissivity, "  [1]"
        write(*,"(a37,g14.6,a)") "lw_clear_sky_transmissivity_per_km = ", &
                                 c%lw_clear_sky_transmissivity_per_km, "  [km-1]"
        write(*,"(a37,g14.6,a)") "lw_night_cloud_fraction = ", c%lw_night_cloud_fraction, "  [1]"
        write(*,"(a25,g14.6,a)") "eps_snow = ", c%eps_snow, "  [1]"
        write(*,"(a25,g14.6,a)") "eps_ice  = ", c%eps_ice,  "  [1]"
        write(*,"(a25,g14.6,a)") "sigma_sb = ", c%sigma_sb, "  [W m-2 K-4]"
        write(*,"(a25,g14.6,a)") "T0       = ", c%T0,       "  [K]"
        write(*,"(a25,i14)")     "low_density_densification = ", c%low_density_densification

        return

    end subroutine chion_const_print

    subroutine chion_const_validate(c)
        ! Constraints Chion.jl's SnowpackPhysicalConstants constructor checks
        ! (src/constants.jl, 6d06af6, 03bb445): the aging timescales are
        ! positive, under the aging scheme 0 <= alpha_wet <= alpha_dry <= 1,
        ! the clear-sky transmissivity is positive and the night cloud
        ! fraction in [0,1], the semix turbulence's karman constant, height,
        ! roughness lengths, roughness ratio and exchange factors positive and
        ! its stable coefficient non-negative. Plus chion's
        ! aging_snowfall_ref > 0 (D30).

        implicit none

        type(chion_const_class), intent(IN) :: c

        if (.not. (c%semix_karman .gt. 0.0_wp .and. c%semix_surface_height .gt. 0.0_wp &
                   .and. c%semix_z0m_snow .gt. 0.0_wp .and. c%semix_z0m_ice .gt. 0.0_wp &
                   .and. c%semix_zm_to_zh .gt. 0.0_wp &
                   .and. c%semix_sensible_exchange_factor .gt. 0.0_wp &
                   .and. c%semix_latent_exchange_factor .gt. 0.0_wp)) then
            write(io_unit_err,*) "chion_const_validate:: Error: semix_karman, &
                                 &semix_surface_height, semix_z0m_snow, semix_z0m_ice, &
                                 &semix_zm_to_zh and the semix exchange factors must be positive."
            stop "Program stopped."
        end if

        if (.not. c%semix_stable_coefficient .ge. 0.0_wp) then
            write(io_unit_err,*) "chion_const_validate:: Error: semix_stable_coefficient &
                                 &must be non-negative."
            write(io_unit_err,*) "semix_stable_coefficient = ", c%semix_stable_coefficient
            stop "Program stopped."
        end if

        if (.not. c%lw_clear_sky_transmissivity .gt. 0.0_wp) then
            write(io_unit_err,*) "chion_const_validate:: Error: lw_clear_sky_transmissivity &
                                 &must be positive."
            write(io_unit_err,*) "lw_clear_sky_transmissivity = ", c%lw_clear_sky_transmissivity
            stop "Program stopped."
        end if

        if (c%lw_night_cloud_fraction .lt. 0.0_wp .or. c%lw_night_cloud_fraction .gt. 1.0_wp) then
            write(io_unit_err,*) "chion_const_validate:: Error: lw_night_cloud_fraction &
                                 &must be in [0,1]."
            write(io_unit_err,*) "lw_night_cloud_fraction = ", c%lw_night_cloud_fraction
            stop "Program stopped."
        end if

        if (c%aging_cold_timescale_days .le. 0.0_wp .or. &
            c%aging_melting_timescale_days .le. 0.0_wp) then
            write(io_unit_err,*) "chion_const_validate:: Error: aging timescales must be positive."
            write(io_unit_err,*) "aging_cold_timescale_days    = ", c%aging_cold_timescale_days
            write(io_unit_err,*) "aging_melting_timescale_days = ", c%aging_melting_timescale_days
            stop "Program stopped."
        end if

        if (c%aging_snowfall_ref .le. 0.0_wp) then
            write(io_unit_err,*) "chion_const_validate:: Error: aging_snowfall_ref must be positive."
            write(io_unit_err,*) "aging_snowfall_ref = ", c%aging_snowfall_ref
            stop "Program stopped."
        end if

        if (c%albedo_scheme .eq. CHION_ALBEDO_AGING) then
            if (c%alpha_wet .lt. 0.0_wp .or. c%alpha_wet .gt. c%alpha_dry &
                                         .or. c%alpha_dry .gt. 1.0_wp) then
                write(io_unit_err,*) "chion_const_validate:: Error: albedo_scheme = 'aging' &
                                     &requires 0 <= alpha_wet <= alpha_dry <= 1."
                write(io_unit_err,*) "alpha_wet, alpha_dry = ", c%alpha_wet, c%alpha_dry
                stop "Program stopped."
            end if
        end if

        return

    end subroutine chion_const_validate

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

        allocate(forc%toa_shortwave(ncol))
        allocate(forc%has_toa_shortwave(ncol))

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

        forc%toa_shortwave     = 0.0_wp
        forc%has_toa_shortwave = .FALSE.

        ! H_ice = 0 selects ITM's land albedo branch and gives BESSI no ice
        ! substrate; PDDs = 0 selects ITM's "desert" critical snow depth. A
        ! host running model="itm", or BESSI with ice_substrate_layers > 0,
        ! must set H_ice.
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
        if (allocated(forc%toa_shortwave))         deallocate(forc%toa_shortwave)
        if (allocated(forc%has_toa_shortwave))     deallocate(forc%has_toa_shortwave)
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
            case("warren_wiscombe")
                flag = SEMIX_SNOW_ALBEDO_WW
            case("dang")
                flag = SEMIX_SNOW_ALBEDO_DANG
            case DEFAULT
                write(io_unit_err,*) "chion_semix_snow_albedo_flag:: Error: scheme not recognized."
                write(io_unit_err,*) "semix_snow_albedo should be one of: ['warren_wiscombe','dang']"
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

    function chion_longwave_scheme_flag(name) result(flag)
        ! Map a namelist string onto a downwelling-longwave scheme flag.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("graybody")
                flag = CHION_LONGWAVE_GRAYBODY
            case("cloud_proxy")
                flag = CHION_LONGWAVE_CLOUD_PROXY
            case DEFAULT
                write(io_unit_err,*) "chion_longwave_scheme_flag:: Error: longwave scheme not recognized."
                write(io_unit_err,*) "longwave_scheme should be one of: ['graybody','cloud_proxy']"
                write(io_unit_err,*) "longwave_scheme = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_longwave_scheme_flag

    function chion_turbulent_flux_scheme_flag(name) result(flag)
        ! Map a namelist string onto a turbulent-flux scheme flag.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("bessi")
                flag = CHION_TURB_BESSI
            case("semix")
                flag = CHION_TURB_SEMIX
            case("climberx")
                flag = CHION_TURB_CLIMBERX
            case DEFAULT
                write(io_unit_err,*) "chion_turbulent_flux_scheme_flag:: Error: &
                                     &turbulent flux scheme not recognized."
                write(io_unit_err,*) "turbulent_flux_scheme should be one of: &
                                     &['bessi','semix','climberx']"
                write(io_unit_err,*) "turbulent_flux_scheme = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_turbulent_flux_scheme_flag

    function chion_climberx_qsat_flag(name) result(flag)
        ! Map a namelist string onto a saturation-humidity parameterization.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("climberx")
                flag = CLIMBERX_QSAT_CLIMBERX
            case("bessi")
                flag = CLIMBERX_QSAT_BESSI
            case DEFAULT
                write(io_unit_err,*) "chion_climberx_qsat_flag:: Error: scheme not recognized."
                write(io_unit_err,*) "climberx_qsat should be one of: ['climberx','bessi']"
                write(io_unit_err,*) "climberx_qsat = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_climberx_qsat_flag

    function chion_albedo_scheme_flag(name) result(flag)
        ! Map a namelist string onto an albedo scheme flag. Canonical names
        ! only: Chion.jl 03bb445 dropped the aliases :bessi and :legacy.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("constant")
                flag = CHION_ALBEDO_CONSTANT
            case("dynamic")
                flag = CHION_ALBEDO_DYNAMIC
            case("prescribed")
                flag = CHION_ALBEDO_PRESCRIBED
            case("semix")
                flag = CHION_ALBEDO_SEMIX
            case("aging")
                flag = CHION_ALBEDO_AGING
            case DEFAULT
                write(io_unit_err,*) "chion_albedo_scheme_flag:: Error: albedo scheme not recognized."
                write(io_unit_err,*) "albedo_scheme should be one of: &
                                     &['constant','dynamic','prescribed','semix','aging']"
                write(io_unit_err,*) "albedo_scheme = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function chion_albedo_scheme_flag

    function chion_fresh_snow_density_scheme_flag(name) result(flag)
        ! Canonical names only: Chion.jl 03bb445 dropped the aliases
        ! :bessi -> :constant and :htessel -> :parameterized.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("constant")
                flag = CHION_FRESH_SNOW_DENSITY_CONSTANT
            case("parameterized")
                flag = CHION_FRESH_SNOW_DENSITY_PARAMETERIZED
            case DEFAULT
                write(io_unit_err,*) "chion_fresh_snow_density_scheme_flag:: Error: &
                                     &fresh snow density scheme not recognized."
                write(io_unit_err,*) "fresh_snow_density_scheme should be one of: &
                                     &['constant','parameterized']"
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
