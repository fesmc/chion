program test_surface
    ! WP5 (surface-fluxes half) acceptance test: snow_surface_fluxes.
    !
    ! Covers:
    !   * every has_* flag branch selects the prescribed value over the
    !     internal parameterization
    !   * the linearized latent flux reproduces the exact one at T = T^n
    !     (consistency at the linearization point -- the ONLY place where the
    !     two evaluations of trap 2 are required to agree)
    !   * saturation vapor pressures against hand-computed values
    !   * bare-ice energy-balance closure and the vapor_mass sign convention
    !   * diagnose_latent_heat_flux_coefficients in all three branches,
    !     including snowfall beating rainfall
    !   * the cloud-proxy downwelling longwave (Chion.jl 03bb445): the
    !     emissivity at known inputs, the night / sub-daily / no-latitude
    !     fallback, a missing or negative surface height, the clamp, the
    !     host TOA override (D33), and the pass-through of prescribed
    !     longwave and of the graybody scheme
    !   * phase-dependent latent heat and the gradient-based vapour mass
    !     (Chion.jl d0146e1): L = Lv+Lm below T0, Lv at it, Lv+Lm on bare
    !     ice; apply_snow_surface_vapor_mass_flux applies E*dt in both the
    !     solid and the liquid branch, and q_lh*dt/L(Ts) when prescribed

    use chion_defs,          only : wp, wp_acc, chion_const_class, &
                                    chion_step_forcing_class, chion_const_init, &
                                    DEF_SEA_LEVEL_AIR_PRESSURE, &
                                    CHION_SEB_BESSI, CHION_SEB_SEMIX, &
                                    CHION_TURB_BESSI, CHION_TURB_SEMIX, CHION_TURB_CLIMBERX, &
                                    CHION_LONGWAVE_GRAYBODY, CHION_LONGWAVE_CLOUD_PROXY
    use snow_surface_fluxes
    use snow_diurnal, only : daily_toa_shortwave
    use snow_vapor
    use snow_seb_semix
    use snow_turbulence

    implicit none

    type(chion_const_class)        :: c
    type(chion_step_forcing_class) :: forc
    type(nonshortwave_flux_class)  :: nsw, nsw_ice
    type(bare_ice_flux_class)      :: bif
    type(bare_ice_ablation_class)  :: abl
    type(latent_vapor_flux_lin_class) :: lin
    type(latent_heat_coeff_class)     :: coef

    ! Snow depth, needed only by turbulent_flux_scheme = "climberx" for its
    ! roughness blend. Every check below except the SEMIX section runs the
    ! BESSI scheme, which ignores it entirely.
    real(wp), parameter :: H_NONE = 0.0_wp
    real(wp), parameter :: H_DEEP = 1.0_wp

    real(wp) :: Tn, q_exact, q_lin, dt_seconds, q_net, expected, E0
    real(wp) :: f_sh_semix
    integer  :: nfail

    type(semix_exchange_class) :: sx
    type(turb_semix_lin_class) :: tx

    nfail = 0

    ! BESSI's original surface scheme, which the closed forms below are
    ! written for (the defaults are Chion.jl's calibrated semix set since
    ! C11); the blocks for the other schemes switch explicitly.
    call chion_const_init(c)
    c%seb_scheme            = CHION_SEB_BESSI
    c%turbulent_flux_scheme = CHION_TURB_BESSI
    c%longwave_scheme       = CHION_LONGWAVE_GRAYBODY

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP5 acceptance test: snow_surface_fluxes"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! === Vapor-pressure parameterizations ================================
    write(*,"(a)") "--- vapor pressures (hand-computed) ---"

    ! At T = T0 the exponent vanishes in both forms, so both must return the
    ! shared prefactor 611.2 Pa exactly.
    call check_close("es_water(273.15) = 611.2", &
                     water_saturation_vapor_pressure(273.15_wp,c%T0), 611.2_wp, 1.0e-6_wp, nfail)
    call check_close("es_ice(273.15)   = 611.2", &
                     ice_saturation_vapor_pressure(273.15_wp,c%T0), 611.2_wp, 1.0e-6_wp, nfail)

    ! 611.2*exp(17.27*(-10)/(-10+243.12)) = 291.3729493
    call check_close("es_water(263.15) = 291.37295", &
                     water_saturation_vapor_pressure(263.15_wp,c%T0), 291.3729493_wp, 1.0e-6_wp, nfail)

    ! 611.2*exp(22.46*(-10)/(-10+272.62)) = 259.8738060
    call check_close("es_ice(263.15)   = 259.87381", &
                     ice_saturation_vapor_pressure(263.15_wp,c%T0), 259.8738060_wp, 1.0e-6_wp, nfail)

    ! des/dT = es*22.46*272.62/(T-T0+272.62)^2 = 23.0714228 Pa K-1 at 263.15 K
    call check_close("des_ice/dT(263.15) = 23.071423", &
                     ice_saturation_vapor_pressure_derivative(263.15_wp,c%T0, &
                         ice_saturation_vapor_pressure(263.15_wp,c%T0)), &
                     23.0714228_wp, 1.0e-6_wp, nfail)

    ! Saturation over ice is below saturation over water at subfreezing
    ! temperatures -- the whole point of using two different fits.
    call check("es_ice < es_water below freezing", &
               ice_saturation_vapor_pressure(263.15_wp,c%T0) .lt. &
               water_saturation_vapor_pressure(263.15_wp,c%T0), nfail)

    write(*,*)
    write(*,"(a)") "--- safe_positive and relative_humidity_fraction ---"

    call check_close("safe_positive passes a normal value", &
                     safe_positive(1.0e5_wp), 1.0e5_wp, 1.0e-6_wp, nfail)
    call check_close("safe_positive floors zero at 1e-12", &
                     safe_positive(0.0_wp), 1.0e-12_wp, 1.0e-6_wp, nfail)
    call check_close("safe_positive floors negatives at 1e-12", &
                     safe_positive(-5.0_wp), 1.0e-12_wp, 1.0e-6_wp, nfail)

    call check_close("rh = 0.6 stays a fraction", &
                     relative_humidity_fraction(0.6_wp), 0.6_wp, 1.0e-6_wp, nfail)
    call check_close("rh = 60 is read as percent", &
                     relative_humidity_fraction(60.0_wp), 0.6_wp, 1.0e-6_wp, nfail)
    call check_close("rh = 150 clamps to 1", &
                     relative_humidity_fraction(150.0_wp), 1.0_wp, 1.0e-6_wp, nfail)
    call check_close("rh < 0 clamps to 0", &
                     relative_humidity_fraction(-0.3_wp), 0.0_wp, 1.0e-6_wp, nfail)

    ! === Linearization consistency at T = T^n ============================
    write(*,*)
    write(*,"(a)") "--- linearized latent flux at the linearization point ---"

    Tn = 265.0_wp

    q_exact = latent_vapor_flux(Tn,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE, &
                                surface_vapor_latent_heat(Tn,c))
    lin     = latent_vapor_flux_linearized(Tn,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE)
    q_lin   = lin%constant - lin%linear*Tn

    call check_close("Q_lin(T^n) = Q_exact(T^n)", q_lin, q_exact, 1.0e-5_wp, nfail)

    ! The linear coefficient must be the exchange-weighted des/dT, which is
    ! strictly positive: a warmer surface loses more vapor.
    call check("linear coefficient is positive", lin%linear .gt. 0.0_wp, nfail)

    ! Away from T^n the two DO differ -- this is trap 2, and it is deliberate.
    call check("linearized and exact differ away from T^n (trap 2 preserved)", &
               abs((lin%constant - lin%linear*(Tn+5.0_wp)) &
                   - latent_vapor_flux(Tn+5.0_wp,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE, &
                                       surface_vapor_latent_heat(Tn+5.0_wp,c))) &
               .gt. 1.0e-4_wp, nfail)

    ! === Phase-dependent latent heat, gradient vapour mass ===============
    write(*,*)
    write(*,"(a)") "--- phase-dependent latent heat (Chion.jl d0146e1) ---"

    call check_close("L below T0 = Lv + Lm", surface_vapor_latent_heat(c%T0 - 0.01_wp,c), &
                     c%Lv + c%Lm, 1.0e-6_wp, nfail)
    call check_close("L at T0 = Lv", surface_vapor_latent_heat(c%T0,c), c%Lv, 1.0e-6_wp, nfail)

    E0 = vapor_mass_flux(c%T0,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE)
    call check("E(T0) < 0 in subsaturated cold air (evaporation)", E0 .lt. 0.0_wp, nfail)
    call check_close("Q(T0) = Lv*E(T0)", &
                     latent_vapor_flux(c%T0,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE,c%Lv), &
                     c%Lv*E0, 1.0e-5_wp, nfail)

    lin   = latent_vapor_flux_linearized(c%T0,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE)
    q_lin = lin%constant - lin%linear*c%T0
    call check_close("linearized at T0 carries Lv: Q_lin(T0) = Lv*E(T0)", &
                     q_lin, c%Lv*E0, 1.0e-5_wp, nfail)

    call forcing_init(forc)
    forc%air_temperature       = 270.0_wp
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,c%T0,H_NONE,.TRUE.)
    call check_close("melting snow surface: Q = Lv*E(T0)", nsw%latent, c%Lv*E0, &
                     1.0e-5_wp, nfail)
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,c%T0,H_NONE,.FALSE.)
    call check_close("bare ice at T0 stays solid: Q = (Lv+Lm)*E(T0)", nsw%latent, &
                     (c%Lv + c%Lm)*E0, 1.0e-5_wp, nfail)
    abl = bare_ice_ablation_mass(c,forc,86400.0_wp,c%alpha_ice)
    call check_close("bare ice vapour mass = E(T0)*dt", abl%vapor_mass, E0*86400.0_wp, &
                     1.0e-5_wp, nfail)

    call test_snow_vapor_mass(nfail)
    call test_cloud_proxy(nfail)

    ! === has_* flag branches ============================================
    write(*,*)
    write(*,"(a)") "--- has_* flags select prescribed over parameterized ---"

    call forcing_init(forc)
    forc%air_temperature = 268.0_wp
    forc%rainfall_rate   = 0.0_wp

    ! Baseline: nothing prescribed, no humidity.
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)

    call check_close("LW parameterized = sigma*(eps_a*Ta^4 - eps_s*Ts^4)", &
                     nsw%longwave, &
                     c%sigma_sb*(c%eps_air*268.0_wp**4 - c%eps_snow*265.0_wp**4), &
                     1.0e-5_wp, nfail)
    call check_close("SH parameterized = D_sh*(Ta - Ts)", &
                     nsw%sensible, c%D_sh*3.0_wp, 1.0e-5_wp, nfail)
    call check_close("LH = 0 with neither q_lh nor humidity", &
                     nsw%latent, 0.0_wp, 1.0e-6_wp, nfail)
    call check_close("Q_rain = 0 with no rainfall", nsw%rain, 0.0_wp, 1.0e-6_wp, nfail)

    ! has_q_lw_down: the downward term is replaced, the upward term is kept.
    forc%has_q_lw_down = .TRUE.
    forc%q_lw_down     = 250.0_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)
    call check_close("has_q_lw_down -> q_lw_down - sigma*eps_s*Ts^4", &
                     nsw%longwave, 250.0_wp - c%sigma_sb*c%eps_snow*265.0_wp**4, &
                     1.0e-5_wp, nfail)

    ! has_q_sh: taken verbatim, no dependence on Ta or Ts at all.
    forc%has_q_sh = .TRUE.
    forc%q_sh     = -17.5_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)
    call check_close("has_q_sh -> prescribed value verbatim", &
                     nsw%sensible, -17.5_wp, 1.0e-6_wp, nfail)

    ! has_relative_humidity alone activates the internal vapor flux.
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)
    call check_close("has_relative_humidity -> BESSI vapor flux", &
                     nsw%latent, &
                     latent_vapor_flux(265.0_wp,c,268.0_wp,0.75_wp,forc%air_pressure, &
                                       c%Lv + c%Lm), &
                     1.0e-5_wp, nfail)

    ! has_q_lh takes precedence over has_relative_humidity.
    forc%has_q_lh = .TRUE.
    forc%q_lh     = -8.25_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)
    call check_close("has_q_lh beats has_relative_humidity", &
                     nsw%latent, -8.25_wp, 1.0e-6_wp, nfail)
    call check_close("resolved_turbulent_latent_heat_flux agrees", &
                     resolved_turbulent_latent_heat_flux(c,forc,265.0_wp,H_NONE,.TRUE.), -8.25_wp, &
                     1.0e-6_wp, nfail)

    ! Rain heat flux is always parameterized; there is no has_* flag for it.
    forc%rainfall_rate = 1.0e-4_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)
    call check_close("Q_rain = P_rain*cw*(Ta - T0)", &
                     nsw%rain, 1.0e-4_wp*c%cw*(268.0_wp - c%T0), 1.0e-5_wp, nfail)

    ! === CLIMBER-X SEMIX: seb_scheme = "semix", turbulence = "climberx" ===
    ! The CLIMBER-X aerodynamic scheme replaces the sensible and
    ! turbulent-latent terms at BOTH exact-flux sites, its longwave the
    ! longwave. Its own physics is pinned in test_seb; what is checked here is
    ! the DISPATCH: that these two entry points route to it, that rain is
    ! untouched by the switches, and that a prescribed flux still wins.
    write(*,*)
    write(*,"(a)") "--- seb_scheme = semix, turbulent_flux_scheme = climberx dispatch ---"

    call forcing_init(forc)
    forc%air_temperature       = 268.0_wp
    forc%wind_speed            = 5.0_wp
    forc%rainfall_rate         = 0.0_wp
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp

    c%seb_scheme            = CHION_SEB_SEMIX
    c%turbulent_flux_scheme = CHION_TURB_CLIMBERX

    sx  = semix_turbulent_exchange(c,H_DEEP,268.0_wp,265.0_wp,5.0_wp, &
                                   forc%air_pressure,0.75_wp,.TRUE.)
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)

    call check_close("semix SH = f_sh*(Ta - Ts)", nsw%sensible, &
                     semix_sensible_heat_flux(sx,268.0_wp,265.0_wp), 1.0e-5_wp, nfail)
    call check_close("semix LH = -f_lh*(qsat - q_air)", nsw%latent, &
                     semix_latent_heat_flux(sx), 1.0e-5_wp, nfail)
    call check_close("semix latent site agrees with the flux-components site", &
                     resolved_turbulent_latent_heat_flux(c,forc,265.0_wp,H_DEEP,.TRUE.), &
                     nsw%latent, 1.0e-6_wp, nfail)

    ! Longwave: semix absorbs the downwelling flux with the surface
    ! emissivity, BESSI absorbs it in full. Same emitted term either way.
    call check_close("semix LW = emiss*(LWdn - sigma*Ts^4)", nsw%longwave, &
                     c%eps_snow*(c%sigma_sb*c%eps_air*268.0_wp**4) &
                     - c%eps_snow*c%sigma_sb*265.0_wp**4, 1.0e-5_wp, nfail)
    call check("semix absorbs LESS longwave than bessi", &
               nsw%longwave .lt. c%sigma_sb*(c%eps_air*268.0_wp**4 &
                                             - c%eps_snow*265.0_wp**4), nfail)

    ! Snow vs ice emissivity reaches the longwave term (ebal's mask_snow).
    c%eps_ice = 0.90_wp
    nsw_ice = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.FALSE.)
    call check_close("has_snow = F selects eps_ice", nsw_ice%longwave, &
                     c%eps_ice*(c%sigma_sb*c%eps_air*268.0_wp**4) &
                     - c%eps_ice*c%sigma_sb*265.0_wp**4, 1.0e-5_wp, nfail)
    c%eps_ice = 0.98_wp

    ! A prescribed downwelling flux is absorbed with emissivity too.
    forc%has_q_lw_down = .TRUE.
    forc%q_lw_down     = 250.0_wp
    nsw_ice = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check_close("semix absorbs a prescribed LWdn with emissivity", &
                     nsw_ice%longwave, &
                     c%eps_snow*250.0_wp - c%eps_snow*c%sigma_sb*265.0_wp**4, &
                     1.0e-5_wp, nfail)
    forc%has_q_lw_down = .FALSE.

    call check_close("semix leaves the rain flux alone", nsw%rain, 0.0_wp, &
                     1.0e-6_wp, nfail)

    ! Snow depth actually reaches the roughness blend: deeper snow is rougher
    ! (z0m_snow > z0m_ice), so it exchanges more strongly.
    f_sh_semix = nsw%sensible
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE,.TRUE.)
    call check("h_snow reaches the roughness blend", &
               nsw%sensible .lt. f_sh_semix, nfail)

    ! Prescribed fluxes still beat the scheme, exactly as under BESSI.
    forc%has_q_sh = .TRUE.
    forc%q_sh     = -17.5_wp
    forc%has_q_lh = .TRUE.
    forc%q_lh     = -8.25_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check_close("has_q_sh still beats the semix scheme", nsw%sensible, &
                     -17.5_wp, 1.0e-6_wp, nfail)
    call check_close("has_q_lh still beats the semix scheme", nsw%latent, &
                     -8.25_wp, 1.0e-6_wp, nfail)

    ! Bare ice takes the same branch, at h_snow = 0 and T = T0.
    call forcing_init(forc)
    forc%air_temperature = 271.0_wp
    forc%wind_speed      = 5.0_wp
    bif = resolved_bare_ice_surface_flux_components(c,forc,0.30_wp)
    sx  = semix_turbulent_exchange(c,0.0_wp,271.0_wp,c%T0,5.0_wp, &
                                   forc%air_pressure,0.0_wp,.FALSE.)
    call check_close("bare ice takes the semix sensible flux at T0", &
                     bif%sensible, semix_sensible_heat_flux(sx,271.0_wp,c%T0), &
                     1.0e-5_wp, nfail)

    ! The two switches are independent (Chion.jl d0146e1): the semix
    ! longwave with BESSI's turbulence, and BESSI's longwave with the
    ! CLIMBER-X turbulence.
    call forcing_init(forc)
    forc%air_temperature       = 268.0_wp
    forc%wind_speed            = 5.0_wp
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp

    c%seb_scheme            = CHION_SEB_SEMIX
    c%turbulent_flux_scheme = CHION_TURB_BESSI
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check_close("seb semix + bessi turbulence: semix longwave", nsw%longwave, &
                     c%eps_snow*(c%sigma_sb*c%eps_air*268.0_wp**4) &
                     - c%eps_snow*c%sigma_sb*265.0_wp**4, 1.0e-5_wp, nfail)
    call check_close("seb semix + bessi turbulence: SH = D_sh*(Ta - Ts)", nsw%sensible, &
                     c%D_sh*(268.0_wp - 265.0_wp), 1.0e-5_wp, nfail)
    call check_close("seb semix + bessi turbulence: BESSI latent flux", nsw%latent, &
                     latent_vapor_flux(265.0_wp,c,268.0_wp,0.75_wp,forc%air_pressure, &
                                       c%Lv + c%Lm), 1.0e-5_wp, nfail)

    c%seb_scheme            = CHION_SEB_BESSI
    c%turbulent_flux_scheme = CHION_TURB_CLIMBERX
    sx  = semix_turbulent_exchange(c,H_DEEP,268.0_wp,265.0_wp,5.0_wp, &
                                   forc%air_pressure,0.75_wp,.TRUE.)
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check_close("seb bessi + climberx turbulence: BESSI longwave", nsw%longwave, &
                     c%sigma_sb*(c%eps_air*268.0_wp**4 - c%eps_snow*265.0_wp**4), &
                     1.0e-5_wp, nfail)
    call check_close("seb bessi + climberx turbulence: CLIMBER-X sensible", nsw%sensible, &
                     semix_sensible_heat_flux(sx,268.0_wp,265.0_wp), 1.0e-5_wp, nfail)

    c%seb_scheme            = CHION_SEB_BESSI
    c%turbulent_flux_scheme = CHION_TURB_BESSI

    ! === turbulent_flux_scheme = "semix" (Chion.jl bulk turbulence) =====
    ! Dispatch only; the scheme is pinned in test_seb. Exact fluxes are
    ! constant - linear*T of the linearization at T, with the snow roughness
    ! on snow and the ice roughness and (D35) latent heat on bare ice.
    write(*,*)
    write(*,"(a)") "--- turbulent_flux_scheme = semix dispatch ---"

    call forcing_init(forc)
    forc%air_temperature       = 268.0_wp
    forc%wind_speed            = 5.0_wp
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp

    c%turbulent_flux_scheme = CHION_TURB_SEMIX

    tx  = turb_semix_flux_linearized(c,265.0_wp,268.0_wp,0.75_wp,forc%air_pressure, &
                                     5.0_wp,c%semix_z0m_snow, &
                                     surface_vapor_latent_heat(265.0_wp,c))
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check_close("semix turbulence on snow: SH = const - lin*Ts", nsw%sensible, &
                     tx%sensible_constant - tx%sensible_linear*265.0_wp, 1.0e-5_wp, nfail)
    call check_close("semix turbulence on snow: LH = const - lin*Ts", nsw%latent, &
                     tx%latent_constant - tx%latent_linear*265.0_wp, 1.0e-5_wp, nfail)
    call check_close("semix turbulence: latent site agrees with the flux-components site", &
                     resolved_turbulent_latent_heat_flux(c,forc,265.0_wp,H_DEEP,.TRUE.), &
                     nsw%latent, 1.0e-6_wp, nfail)
    call check_close("semix turbulence leaves the BESSI longwave alone", nsw%longwave, &
                     c%sigma_sb*(c%eps_air*268.0_wp**4 - c%eps_snow*265.0_wp**4), &
                     1.0e-5_wp, nfail)

    ! Bare ice at T0: ice roughness and the ice latent heat.
    tx  = turb_semix_flux_linearized(c,c%T0,268.0_wp,0.75_wp,forc%air_pressure, &
                                     5.0_wp,c%semix_z0m_ice, &
                                     turb_semix_latent_heat(c,c%T0,.TRUE.))
    bif = resolved_bare_ice_surface_flux_components(c,forc,0.30_wp)
    call check_close("semix turbulence on bare ice: ice roughness (sensible)", bif%sensible, &
                     tx%sensible_constant - tx%sensible_linear*c%T0, 1.0e-5_wp, nfail)
    call check_close("semix turbulence on bare ice: ice latent heat (D35)", bif%latent, &
                     tx%latent_constant - tx%latent_linear*c%T0, 1.0e-5_wp, nfail)
    call check_close("substrate latent site = bare-ice latent", &
                     resolved_turbulent_latent_heat_flux(c,forc,c%T0,H_NONE,.FALSE.), &
                     bif%latent, 1.0e-6_wp, nfail)

    ! No humidity forcing: no latent flux, sensible unaffected.
    forc%has_relative_humidity = .FALSE.
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check("semix turbulence: no humidity forcing -> no latent flux", &
               nsw%latent .eq. 0.0_wp .and. &
               resolved_turbulent_latent_heat_flux(c,forc,265.0_wp,H_DEEP,.TRUE.) .eq. 0.0_wp, &
               nfail)
    call check("semix turbulence: sensible flux without humidity", nsw%sensible .gt. 0.0_wp, nfail)

    ! Prescribed fluxes win.
    forc%has_relative_humidity = .TRUE.
    forc%has_q_sh = .TRUE.
    forc%q_sh     = 12.0_wp
    forc%has_q_lh = .TRUE.
    forc%q_lh     = -8.0_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP,.TRUE.)
    call check("semix turbulence: prescribed SH and LH win", &
               nsw%sensible .eq. 12.0_wp .and. nsw%latent .eq. -8.0_wp, nfail)

    c%turbulent_flux_scheme = CHION_TURB_BESSI

    ! === Shortwave: has_q_sw_net and the max(SWdn,0) clamp ==============
    write(*,*)
    write(*,"(a)") "--- shortwave absorption (bare ice) ---"

    call forcing_init(forc)
    forc%air_temperature = 271.0_wp
    forc%shortwave_down  = 300.0_wp

    bif = resolved_bare_ice_surface_flux_components(c,forc,0.30_wp)
    call check_close("SW_abs = max(SWdn,0)*(1-albedo)", &
                     bif%shortwave_absorbed, 300.0_wp*0.70_wp, 1.0e-5_wp, nfail)

    ! surface_fluxes.jl:74 clamps SWdn at 0; energy_flux.jl:99 does NOT.
    forc%shortwave_down = -300.0_wp
    bif = resolved_bare_ice_surface_flux_components(c,forc,0.30_wp)
    call check_close("negative SWdn is clamped to 0 here (NOT in snow_energy)", &
                     bif%shortwave_absorbed, 0.0_wp, 1.0e-6_wp, nfail)

    forc%shortwave_down = 300.0_wp
    forc%has_q_sw_net   = .TRUE.
    forc%q_sw_net       = 42.0_wp
    bif = resolved_bare_ice_surface_flux_components(c,forc,0.30_wp)
    call check_close("has_q_sw_net -> prescribed value verbatim", &
                     bif%shortwave_absorbed, 42.0_wp, 1.0e-6_wp, nfail)

    ! Albedo is clamped to [0,1] before use.
    forc%has_q_sw_net = .FALSE.
    bif = resolved_bare_ice_surface_flux_components(c,forc,-0.5_wp)
    call check_close("albedo < 0 clamps to 0", bif%shortwave_absorbed, 300.0_wp, &
                     1.0e-5_wp, nfail)
    bif = resolved_bare_ice_surface_flux_components(c,forc,2.0_wp)
    call check_close("albedo > 1 clamps to 1", bif%shortwave_absorbed, 0.0_wp, &
                     1.0e-6_wp, nfail)

    ! Bare ice is ASSUMED at the melting point: the non-shortwave components
    ! must equal those evaluated at T = T0, not at any other temperature.
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,c%T0,H_NONE,.FALSE.)
    call check_close("bare-ice LW is evaluated at T0", bif%longwave, nsw%longwave, &
                     1.0e-5_wp, nfail)
    call check_close("bare-ice SH is evaluated at T0", bif%sensible, nsw%sensible, &
                     1.0e-5_wp, nfail)

    ! === Bare-ice ablation: energy closure and sign convention ==========
    write(*,*)
    write(*,"(a)") "--- bare_ice_ablation_mass ---"

    dt_seconds = 86400.0_wp

    call forcing_init(forc)
    forc%air_temperature       = 278.0_wp
    forc%shortwave_down        = 250.0_wp
    forc%rainfall_rate         = 2.0e-5_wp
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.60_wp

    bif = resolved_bare_ice_surface_flux_components(c,forc,c%alpha_ice)
    abl = bare_ice_ablation_mass(c,forc,dt_seconds,c%alpha_ice)

    q_net = bif%shortwave_absorbed + bif%longwave + bif%sensible + bif%latent + bif%rain

    call check("net surface flux is positive in this melting case", &
               q_net .gt. 0.0_wp, nfail)

    ! Energy closure: melt uses the FULL Q_net. It is NOT reduced by the energy
    ! that went into the vapor exchange -- melt and sublimation are computed
    ! independently from the same Q_net. Preserve this.
    call check_close("melt_mass*Lm = max(Q_net,0)*dt (melt NOT reduced by vapor)", &
                     abl%melt_mass*c%Lm, max(q_net,0.0_wp)*dt_seconds, 1.0e-5_wp, nfail)

    call check_close("latent_heat_flux is the resolved latent component", &
                     abl%latent_heat_flux, bif%latent, 1.0e-6_wp, nfail)

    ! Vapor mass on bare ice always uses (Lv + Lm), regardless of temperature.
    call check_close("vapor_mass = LH*dt/(Lv+Lm) on bare ice", &
                     abl%vapor_mass, bif%latent*dt_seconds/(c%Lv + c%Lm), 1.0e-5_wp, nfail)

    call check_close("net_mass_change = vapor_mass - melt_mass", &
                     abl%net_mass_change, abl%vapor_mass - abl%melt_mass, 1.0e-5_wp, nfail)

    ! Sign convention, dry air: ea < es -> latent flux negative -> mass is lost
    ! by sublimation, so vapor_mass < 0 and sublimation_mass = -vapor_mass > 0.
    forc%relative_humidity = 0.10_wp
    abl = bare_ice_ablation_mass(c,forc,dt_seconds,c%alpha_ice)
    call check("dry air -> negative latent flux", abl%latent_heat_flux .lt. 0.0_wp, nfail)
    call check("dry air -> vapor_mass < 0 (mass loss)", abl%vapor_mass .lt. 0.0_wp, nfail)
    call check_close("sublimation_mass = -vapor_mass", &
                     abl%sublimation_mass, -abl%vapor_mass, 1.0e-6_wp, nfail)

    ! Sign convention, saturated warm air: ea > es -> deposition -> mass gain,
    ! and sublimation_mass is clipped at zero.
    forc%relative_humidity = 1.0_wp
    abl = bare_ice_ablation_mass(c,forc,dt_seconds,c%alpha_ice)
    call check("saturated warm air -> positive latent flux", &
               abl%latent_heat_flux .gt. 0.0_wp, nfail)
    call check("deposition -> vapor_mass > 0 (mass gain)", abl%vapor_mass .gt. 0.0_wp, nfail)
    call check_close("sublimation_mass = 0 under deposition", &
                     abl%sublimation_mass, 0.0_wp, 1.0e-6_wp, nfail)

    ! Melt is floored at zero when the net flux is negative.
    call forcing_init(forc)
    forc%air_temperature = 240.0_wp
    forc%shortwave_down  = 0.0_wp
    abl = bare_ice_ablation_mass(c,forc,dt_seconds,c%alpha_ice)
    call check_close("melt_mass = 0 when Q_net < 0", abl%melt_mass, 0.0_wp, 1.0e-6_wp, nfail)

    ! The albedo is the caller's (snow_bessi: the prescribed one, or the
    ! background, D40): the absorbed shortwave follows it.
    call forcing_init(forc)
    forc%air_temperature        = 275.0_wp
    forc%shortwave_down         = 400.0_wp

    abl = bare_ice_ablation_mass(c,forc,dt_seconds,0.90_wp)
    bif = resolved_bare_ice_surface_flux_components(c,forc,0.90_wp)
    q_net = bif%shortwave_absorbed + bif%longwave + bif%sensible + bif%latent + bif%rain
    call check_close("the albedo argument sets the absorbed shortwave", &
                     abl%melt_mass*c%Lm, max(q_net,0.0_wp)*dt_seconds, 1.0e-5_wp, nfail)

    ! === diagnose_latent_heat_flux_coefficients =========================
    write(*,*)
    write(*,"(a)") "--- diagnose_latent_heat_flux_coefficients (linear,constant) ---"

    ! Branch 1: snowfall.
    coef = diagnose_latent_heat_flux_coefficients(.TRUE.,c,268.0_wp,1.0e-5_wp,0.0_wp)
    call check_close("snowfall -> linear = P_snow*ci", &
                     coef%linear, 1.0e-5_wp*c%ci, 1.0e-5_wp, nfail)
    call check_close("snowfall -> constant = P_snow*ci*Ta", &
                     coef%constant, 1.0e-5_wp*c%ci*268.0_wp, 1.0e-5_wp, nfail)

    ! Branch 1 wins even when it is also raining -- snowfall takes precedence.
    coef = diagnose_latent_heat_flux_coefficients(.TRUE.,c,268.0_wp,1.0e-5_wp,5.0e-5_wp)
    call check_close("snowfall beats rainfall: linear unchanged", &
                     coef%linear, 1.0e-5_wp*c%ci, 1.0e-5_wp, nfail)
    call check_close("snowfall beats rainfall: constant has no rain term", &
                     coef%constant, 1.0e-5_wp*c%ci*268.0_wp, 1.0e-5_wp, nfail)

    ! Snowfall branch does not consult has_surface_snow.
    coef = diagnose_latent_heat_flux_coefficients(.FALSE.,c,268.0_wp,1.0e-5_wp,0.0_wp)
    call check_close("snowfall branch ignores has_surface_snow", &
                     coef%constant, 1.0e-5_wp*c%ci*268.0_wp, 1.0e-5_wp, nfail)

    ! Branch 2: rainfall onto an existing snow surface.
    coef = diagnose_latent_heat_flux_coefficients(.TRUE.,c,278.0_wp,0.0_wp,5.0e-5_wp)
    call check_close("rain on snow -> linear = 0", coef%linear, 0.0_wp, 1.0e-6_wp, nfail)
    call check_close("rain on snow -> constant = P_rain*cw*(Ta-T0)", &
                     coef%constant, 5.0e-5_wp*c%cw*(278.0_wp - c%T0), 1.0e-5_wp, nfail)

    ! Branch 3: rainfall with no surface snow, and the fully quiescent case.
    coef = diagnose_latent_heat_flux_coefficients(.FALSE.,c,278.0_wp,0.0_wp,5.0e-5_wp)
    call check_close("rain without surface snow -> linear = 0", &
                     coef%linear, 0.0_wp, 1.0e-6_wp, nfail)
    call check_close("rain without surface snow -> constant = 0", &
                     coef%constant, 0.0_wp, 1.0e-6_wp, nfail)

    coef = diagnose_latent_heat_flux_coefficients(.TRUE.,c,268.0_wp,0.0_wp,0.0_wp)
    call check_close("no precipitation -> linear = 0", coef%linear, 0.0_wp, 1.0e-6_wp, nfail)
    call check_close("no precipitation -> constant = 0", coef%constant, 0.0_wp, 1.0e-6_wp, nfail)

    ! === Summary ========================================================
    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP5 (surface fluxes): ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP5 (surface fluxes): ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    subroutine test_snow_vapor_mass(nfail)
        ! apply_snow_surface_vapor_mass_flux with humidity on, one step on a
        ! two-layer column, in each reservoir. The applied mass equals the
        ! reported vapor_mass (nothing is clipped here), and for the
        ! parameterized BESSI turbulence it is E(Ts)*dt whatever the phase:
        ! at T0 the old Q/Lv conversion of an (Lv+Lm)-weighted flux would
        ! have moved (Lv+Lm)/Lv = 1.13 times as much water.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_step_forcing_class) :: fv
        type(surface_vapor_flux_class) :: vflux
        real(wp)     :: mass(4), mass_w(4), density(4), temperature(4)
        real(wp)     :: t_srf, albedo, dt, E
        real(wp_acc) :: runoff
        integer      :: n

        dt = 86400.0_wp

        call forcing_init(fv)
        fv%air_temperature       = 268.0_wp
        fv%has_relative_humidity = .TRUE.
        fv%relative_humidity     = 0.5_wp

        ! Solid branch: surface below T0, sublimation from mass(1).
        mass = 0.0_wp; mass_w = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        n = 2
        mass(1:2) = 200.0_wp; density(1:2) = 350.0_wp; temperature(1:2) = 265.0_wp
        t_srf = 265.0_wp; albedo = 0.8_wp; runoff = 0.0_wp_acc

        E = vapor_mass_flux(265.0_wp,c,fv%air_temperature,fv%relative_humidity, &
                            fv%air_pressure)
        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n,runoff, &
                                                t_srf,albedo,c,fv,dt,300.0_wp,100.0_wp,vflux)
        call check_close("solid: vapor_mass = E(Ts)*dt", vflux%vapor_mass, E*dt, 1.0e-5_wp, nfail)
        call check_close("solid: mass(1) changed by vapor_mass", mass(1) - 200.0_wp, &
                         vflux%vapor_mass, 1.0e-3_wp, nfail)
        call check_close("solid: latent flux = (Lv+Lm)*E(Ts)", vflux%latent_heat_flux, &
                         (c%Lv + c%Lm)*E, 1.0e-5_wp, nfail)

        ! Liquid branch: surface at T0, evaporation from mass_w(1).
        mass = 0.0_wp; mass_w = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        n = 2
        mass(1:2) = 200.0_wp; density(1:2) = 350.0_wp; temperature(1:2) = c%T0
        mass_w(1) = 5.0_wp
        t_srf = c%T0; runoff = 0.0_wp_acc

        E = vapor_mass_flux(c%T0,c,fv%air_temperature,fv%relative_humidity,fv%air_pressure)
        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n,runoff, &
                                                t_srf,albedo,c,fv,dt,300.0_wp,100.0_wp,vflux)
        call check_close("liquid: vapor_mass = E(T0)*dt (not Q*dt/Lv)", vflux%vapor_mass, &
                         E*dt, 1.0e-5_wp, nfail)
        call check_close("liquid: mass_w(1) changed by vapor_mass", mass_w(1) - 5.0_wp, &
                         vflux%vapor_mass, 1.0e-4_wp, nfail)
        call check_close("liquid: latent flux = Lv*E(T0)", vflux%latent_heat_flux, &
                         c%Lv*E, 1.0e-5_wp, nfail)

        ! Prescribed q_lh at T0: the flux controls its mass, q_lh*dt/Lv.
        fv%has_q_lh = .TRUE.
        fv%q_lh     = -20.0_wp
        mass_w(1) = 5.0_wp
        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n,runoff, &
                                                t_srf,albedo,c,fv,dt,300.0_wp,100.0_wp,vflux)
        call check_close("prescribed q_lh at T0: vapor_mass = q_lh*dt/Lv", vflux%vapor_mass, &
                         -20.0_wp*dt/c%Lv, 1.0e-5_wp, nfail)

        ! Chion.jl 03bb445: the flux and the reservoir follow the INTERFACE
        ! temperature t_srf, not the top cell's centre temperature(1), and
        ! t_srf is left as the energy solve set it. Interface at T0 over a
        ! subfreezing cell centre: liquid branch, E(T0).
        fv%has_q_lh = .FALSE.
        mass = 0.0_wp; mass_w = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        n = 2
        mass(1:2) = 200.0_wp; density(1:2) = 350.0_wp; temperature(1:2) = 265.0_wp
        mass_w(1) = 5.0_wp
        t_srf = c%T0; runoff = 0.0_wp_acc

        E = vapor_mass_flux(c%T0,c,fv%air_temperature,fv%relative_humidity,fv%air_pressure)
        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n,runoff, &
                                                t_srf,albedo,c,fv,dt,300.0_wp,100.0_wp,vflux)
        call check_close("interface at T0 over T1 < T0: vapor_mass = E(Tsrf)*dt", &
                         vflux%vapor_mass, E*dt, 1.0e-5_wp, nfail)
        call check("interface at T0 over T1 < T0: liquid reservoir, solid untouched", &
                   mass(1) .eq. 200.0_wp .and. mass_w(1) .lt. 5.0_wp, nfail)
        call check("t_srf is not reset to temperature(1)", t_srf .eq. c%T0, nfail)

        ! Chion.jl's semix turbulence converts its own flux with L(Ts):
        ! Lv+Lm from the solid surface, Lv at T0.
        c%turbulent_flux_scheme = CHION_TURB_SEMIX
        fv%wind_speed = 5.0_wp
        mass = 0.0_wp; mass_w = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        n = 2
        mass(1:2) = 200.0_wp; density(1:2) = 350.0_wp; temperature(1:2) = 265.0_wp
        t_srf = 265.0_wp; runoff = 0.0_wp_acc
        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n,runoff, &
                                                t_srf,albedo,c,fv,dt,300.0_wp,100.0_wp,vflux)
        call check_close("semix solid: vapor_mass = Q_lh*dt/(Lv+Lm)", vflux%vapor_mass, &
                         vflux%latent_heat_flux*dt/(c%Lv + c%Lm), 1.0e-5_wp, nfail)
        call check("semix solid: sublimation into dry air", vflux%vapor_mass .lt. 0.0_wp, nfail)

        mass = 0.0_wp; mass_w = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        n = 2
        mass(1:2) = 200.0_wp; density(1:2) = 350.0_wp; temperature(1:2) = c%T0
        mass_w(1) = 5.0_wp
        t_srf = c%T0; runoff = 0.0_wp_acc
        call apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n,runoff, &
                                                t_srf,albedo,c,fv,dt,300.0_wp,100.0_wp,vflux)
        call check_close("semix liquid: vapor_mass = Q_lh*dt/Lv", vflux%vapor_mass, &
                         vflux%latent_heat_flux*dt/c%Lv, 1.0e-5_wp, nfail)
        c%turbulent_flux_scheme = CHION_TURB_BESSI

        return

    end subroutine test_snow_vapor_mass

    subroutine test_cloud_proxy(nfail)
        ! Chion.jl surface_fluxes.jl _cloud_proxy_emissivity and
        ! _with_parameterized_longwave (03bb445), and chion's host TOA (D33).

        use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

        implicit none

        integer, intent(INOUT) :: nfail

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: fp, fc
        real(wp) :: toa, tau, n_cloud, eps, eps_ref, nan_wp

        write(*,*)
        write(*,"(a)") "--- cloud-proxy longwave ---"

        call chion_const_init(c)

        ! Greenland summer day: 70 N at the June solstice, 1.5 km.
        call forcing_init(fp)
        fp%air_temperature     = 268.0_wp
        fp%shortwave_down      = 250.0_wp
        fp%latitude_deg        = 70.0_wp
        fp%solar_longitude_deg = 90.0_wp
        fp%day_of_year         = 172.0_wp
        fp%surface_height      = 1500.0_wp

        toa     = daily_toa_shortwave(70.0_wp,90.0_wp,172.0_wp)
        tau     = 0.85_wp + 0.075_wp*1.5_wp
        n_cloud = 1.0_wp - 250.0_wp/(toa*tau)
        eps_ref = 0.624_wp + 0.0032_wp*(268.0_wp - c%T0) + 0.613_wp*n_cloud

        call check("70 N solstice: TOA well above the 50 W m-2 night limit", &
                   toa .gt. 400.0_wp, nfail)
        call check("70 N solstice: cloudiness strictly inside (0,1)", &
                   n_cloud .gt. 0.0_wp .and. n_cloud .lt. 1.0_wp, nfail)
        eps = cloud_proxy_emissivity(c,fp)
        call check_close("eps = 0.624 + 0.0032(Ta-T0) + 0.613(1 - SW/(TOA tau_clear(z)))", &
                         eps, eps_ref, 1.0e-6_wp, nfail)

        ! Graybody (the default) passes the forcing through untouched.
        c%longwave_scheme = CHION_LONGWAVE_GRAYBODY
        fc = with_parameterized_longwave(c,fp)
        call check("graybody: forcing passed through, no longwave prescribed", &
                   (.not. fc%has_q_lw_down) .and. fc%q_lw_down .eq. fp%q_lw_down, nfail)

        c%longwave_scheme = CHION_LONGWAVE_CLOUD_PROXY
        fc = with_parameterized_longwave(c,fp)
        call check("cloud proxy: longwave resolved as prescribed", fc%has_q_lw_down, nfail)
        call check_close("cloud proxy: q_lw_down = eps*sigma*Ta^4", fc%q_lw_down, &
                         eps_ref*c%sigma_sb*268.0_wp**4, 1.0e-6_wp, nfail)
        call check("cloud proxy: nothing else in the forcing changes", &
                   fc%air_temperature .eq. fp%air_temperature .and. &
                   fc%shortwave_down .eq. fp%shortwave_down, nfail)

        ! A prescribed longwave wins over the proxy.
        fp%has_q_lw_down = .TRUE.
        fp%q_lw_down     = 222.0_wp
        fc = with_parameterized_longwave(c,fp)
        call check("cloud proxy: a prescribed longwave passes through", &
                   fc%has_q_lw_down .and. fc%q_lw_down .eq. 222.0_wp, nfail)
        fp%has_q_lw_down = .FALSE.
        fp%q_lw_down     = 0.0_wp

        ! Night fallbacks: polar night, a sub-daily step, no latitude.
        eps_ref = 0.624_wp + 0.0032_wp*(268.0_wp - c%T0) + 0.613_wp*0.389_wp

        fc = fp
        fc%solar_longitude_deg = 270.0_wp
        fc%day_of_year         = 355.0_wp
        fc%latitude_deg        = 80.0_wp
        call check("80 N midwinter: TOA below the night limit", &
                   daily_toa_shortwave(80.0_wp,270.0_wp,355.0_wp) .le. 50.0_wp, nfail)
        call check_close("polar night: n = lw_night_cloud_fraction", &
                         cloud_proxy_emissivity(c,fc), eps_ref, 1.0e-6_wp, nfail)

        fc = fp
        fc%dt_days = 0.5_wp
        call check_close("sub-daily step: n = lw_night_cloud_fraction", &
                         cloud_proxy_emissivity(c,fc), eps_ref, 1.0e-6_wp, nfail)

        nan_wp = ieee_value(nan_wp,ieee_quiet_nan)
        fc = fp
        fc%latitude_deg = nan_wp
        call check_close("no latitude: n = lw_night_cloud_fraction", &
                         cloud_proxy_emissivity(c,fc), eps_ref, 1.0e-6_wp, nfail)

        ! A missing (NaN) or negative surface height is sea level.
        fc = fp
        fc%surface_height = 0.0_wp
        eps_ref = cloud_proxy_emissivity(c,fc)
        call check("the same SW reads as cloudier aloft (clearer clear sky)", &
                   cloud_proxy_emissivity(c,fp) .gt. eps_ref, nfail)
        fc%surface_height = nan_wp
        call check("NaN surface height = sea level", &
                   cloud_proxy_emissivity(c,fc) .eq. eps_ref, nfail)
        fc%surface_height = -200.0_wp
        call check("negative surface height = sea level", &
                   cloud_proxy_emissivity(c,fc) .eq. eps_ref, nfail)

        ! The host's TOA replaces chion's own (D33), and needs no latitude.
        fc = fp
        fc%has_toa_shortwave = .TRUE.
        fc%toa_shortwave     = 400.0_wp
        fc%latitude_deg      = nan_wp
        eps_ref = 0.624_wp + 0.0032_wp*(268.0_wp - c%T0) &
                  + 0.613_wp*(1.0_wp - 250.0_wp/(400.0_wp*tau))
        call check_close("host TOA replaces the internal one", &
                         cloud_proxy_emissivity(c,fc), eps_ref, 1.0e-6_wp, nfail)
        fc%toa_shortwave = 40.0_wp
        eps_ref = 0.624_wp + 0.0032_wp*(268.0_wp - c%T0) + 0.613_wp*0.389_wp
        call check_close("host TOA below the night limit: night cloudiness", &
                         cloud_proxy_emissivity(c,fc), eps_ref, 1.0e-6_wp, nfail)

        ! Cloudiness and emissivity clamps.
        fc = fp
        fc%shortwave_down = 0.0_wp
        fc%air_temperature = 330.0_wp
        call check_close("emissivity clamped at 1.3", cloud_proxy_emissivity(c,fc), &
                         1.3_wp, 1.0e-6_wp, nfail)
        fc%shortwave_down = 2000.0_wp
        fc%air_temperature = 120.0_wp
        call check_close("emissivity clamped at 0.4 (cloudiness clamped at 0)", &
                         cloud_proxy_emissivity(c,fc), 0.4_wp, 1.0e-6_wp, nfail)

        return

    end subroutine test_cloud_proxy

    subroutine forcing_init(forc)
        ! Neutral per-column forcing: nothing prescribed, no precipitation,
        ! sea-level pressure. Mirrors chion_forcing_alloc's defaults.

        implicit none

        type(chion_step_forcing_class), intent(OUT) :: forc

        forc%air_temperature     = 273.15_wp
        forc%dt_days             = 1.0_wp
        forc%snowfall_rate       = 0.0_wp
        forc%rainfall_rate       = 0.0_wp
        forc%shortwave_down      = 0.0_wp
        forc%wind_speed          = 0.0_wp

        forc%q_sw_net            = 0.0_wp
        forc%q_lw_down           = 0.0_wp
        forc%q_sh                = 0.0_wp
        forc%q_lh                = 0.0_wp

        forc%has_q_sw_net        = .FALSE.
        forc%has_q_lw_down       = .FALSE.
        forc%has_q_sh            = .FALSE.
        forc%has_q_lh            = .FALSE.

        forc%relative_humidity     = 0.0_wp
        forc%has_relative_humidity = .FALSE.

        forc%air_pressure          = DEF_SEA_LEVEL_AIR_PRESSURE
        forc%prescribed_albedo     = 0.0_wp
        forc%has_prescribed_albedo = .FALSE.

        forc%latitude_deg          = 0.0_wp
        forc%day_of_year           = 1.0_wp
        forc%solar_longitude_deg   = 0.0_wp

        return

    end subroutine forcing_init

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

    subroutine check_close(label,value,expected,rtol,nfail)
        ! Relative comparison with an absolute floor, for quantities that come
        ! out of exp() and are therefore sp-limited well above 8*epsilon.

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value
        real(wp),         intent(IN)    :: expected
        real(wp),         intent(IN)    :: rtol
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp) :: tol

        tol = max(abs(expected)*rtol,1.0e-20_wp)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_close

end program test_surface
