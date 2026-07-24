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
    !
    ! apply_snow_surface_vapor_mass_flux is NOT exercised here: it calls
    ! snow_layers (WP4), which does not exist yet. See the WP4 INTEGRATION
    ! POINT markers in src/physics/snow_surface_fluxes.f90.

    use chion_defs,          only : wp, wp_acc, chion_const_class, &
                                    chion_step_forcing_class, chion_const_init, &
                                    CHION_ALBEDO_PRESCRIBED, DEF_SEA_LEVEL_AIR_PRESSURE, &
                                    CHION_SEB_BESSI, CHION_SEB_SEMIX
    use snow_surface_fluxes
    use snow_vapor
    use snow_seb_semix

    implicit none

    type(chion_const_class)        :: c
    type(chion_step_forcing_class) :: forc
    type(nonshortwave_flux_class)  :: nsw
    type(bare_ice_flux_class)      :: bif
    type(bare_ice_ablation_class)  :: abl
    type(latent_vapor_flux_lin_class) :: lin
    type(latent_heat_coeff_class)     :: coef

    ! Snow depth, needed only by seb_scheme = "semix" for its roughness blend.
    ! Every check below except the SEMIX section runs the BESSI scheme, which
    ! ignores it entirely.
    real(wp), parameter :: H_NONE = 0.0_wp
    real(wp), parameter :: H_DEEP = 1.0_wp

    real(wp) :: Tn, q_exact, q_lin, dt_seconds, q_net, expected
    real(wp) :: f_sh_semix
    integer  :: nfail

    type(semix_exchange_class) :: sx

    nfail = 0

    call chion_const_init(c)

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

    q_exact = latent_vapor_flux(Tn,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE)
    lin     = latent_vapor_flux_linearized(Tn,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE)
    q_lin   = lin%constant - lin%linear*Tn

    call check_close("Q_lin(T^n) = Q_exact(T^n)", q_lin, q_exact, 1.0e-5_wp, nfail)

    ! The linear coefficient must be the exchange-weighted des/dT, which is
    ! strictly positive: a warmer surface loses more vapor.
    call check("linear coefficient is positive", lin%linear .gt. 0.0_wp, nfail)

    ! Away from T^n the two DO differ -- this is trap 2, and it is deliberate.
    call check("linearized and exact differ away from T^n (trap 2 preserved)", &
               abs((lin%constant - lin%linear*(Tn+5.0_wp)) &
                   - latent_vapor_flux(Tn+5.0_wp,c,270.0_wp,0.75_wp,DEF_SEA_LEVEL_AIR_PRESSURE)) &
               .gt. 1.0e-4_wp, nfail)

    ! === has_* flag branches ============================================
    write(*,*)
    write(*,"(a)") "--- has_* flags select prescribed over parameterized ---"

    call forcing_init(forc)
    forc%air_temperature = 268.0_wp
    forc%rainfall_rate   = 0.0_wp

    ! Baseline: nothing prescribed, no humidity.
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)

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
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)
    call check_close("has_q_lw_down -> q_lw_down - sigma*eps_s*Ts^4", &
                     nsw%longwave, 250.0_wp - c%sigma_sb*c%eps_snow*265.0_wp**4, &
                     1.0e-5_wp, nfail)

    ! has_q_sh: taken verbatim, no dependence on Ta or Ts at all.
    forc%has_q_sh = .TRUE.
    forc%q_sh     = -17.5_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)
    call check_close("has_q_sh -> prescribed value verbatim", &
                     nsw%sensible, -17.5_wp, 1.0e-6_wp, nfail)

    ! has_relative_humidity alone activates the internal vapor flux.
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)
    call check_close("has_relative_humidity -> BESSI vapor flux", &
                     nsw%latent, &
                     latent_vapor_flux(265.0_wp,c,268.0_wp,0.75_wp,forc%air_pressure), &
                     1.0e-5_wp, nfail)

    ! has_q_lh takes precedence over has_relative_humidity.
    forc%has_q_lh = .TRUE.
    forc%q_lh     = -8.25_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)
    call check_close("has_q_lh beats has_relative_humidity", &
                     nsw%latent, -8.25_wp, 1.0e-6_wp, nfail)
    call check_close("resolved_turbulent_latent_heat_flux agrees", &
                     resolved_turbulent_latent_heat_flux(c,forc,265.0_wp,H_NONE), -8.25_wp, &
                     1.0e-6_wp, nfail)

    ! Rain heat flux is always parameterized; there is no has_* flag for it.
    forc%rainfall_rate = 1.0e-4_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)
    call check_close("Q_rain = P_rain*cw*(Ta - T0)", &
                     nsw%rain, 1.0e-4_wp*c%cw*(268.0_wp - c%T0), 1.0e-5_wp, nfail)

    ! === seb_scheme = "semix" ===========================================
    ! The SEMIX aerodynamic scheme replaces the sensible and turbulent-latent
    ! terms at BOTH exact-flux sites. Its own physics is pinned in test_seb;
    ! what is checked here is the DISPATCH: that these two entry points route
    ! to it, that longwave and rain are untouched by the switch, and that a
    ! prescribed flux still wins.
    write(*,*)
    write(*,"(a)") "--- seb_scheme = semix dispatch ---"

    call forcing_init(forc)
    forc%air_temperature       = 268.0_wp
    forc%wind_speed            = 5.0_wp
    forc%rainfall_rate         = 0.0_wp
    forc%has_relative_humidity = .TRUE.
    forc%relative_humidity     = 0.75_wp

    c%seb_scheme = CHION_SEB_SEMIX

    sx  = semix_turbulent_exchange(c,H_DEEP,268.0_wp,265.0_wp,5.0_wp, &
                                   forc%air_pressure,0.75_wp,.TRUE.)
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP)

    call check_close("semix SH = f_sh*(Ta - Ts)", nsw%sensible, &
                     semix_sensible_heat_flux(sx,268.0_wp,265.0_wp), 1.0e-5_wp, nfail)
    call check_close("semix LH = -f_lh*(qsat - q_air)", nsw%latent, &
                     semix_latent_heat_flux(sx), 1.0e-5_wp, nfail)
    call check_close("semix latent site agrees with the flux-components site", &
                     resolved_turbulent_latent_heat_flux(c,forc,265.0_wp,H_DEEP), &
                     nsw%latent, 1.0e-6_wp, nfail)

    ! The switch is confined to the turbulent terms.
    call check_close("semix leaves longwave alone", nsw%longwave, &
                     c%sigma_sb*(c%eps_air*268.0_wp**4 - c%eps_snow*265.0_wp**4), &
                     1.0e-5_wp, nfail)
    call check_close("semix leaves the rain flux alone", nsw%rain, 0.0_wp, &
                     1.0e-6_wp, nfail)

    ! Snow depth actually reaches the roughness blend: deeper snow is rougher
    ! (z0m_snow > z0m_ice), so it exchanges more strongly.
    f_sh_semix = nsw%sensible
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_NONE)
    call check("h_snow reaches the roughness blend", &
               nsw%sensible .lt. f_sh_semix, nfail)

    ! Prescribed fluxes still beat the scheme, exactly as under BESSI.
    forc%has_q_sh = .TRUE.
    forc%q_sh     = -17.5_wp
    forc%has_q_lh = .TRUE.
    forc%q_lh     = -8.25_wp
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,265.0_wp,H_DEEP)
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

    c%seb_scheme = CHION_SEB_BESSI

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
    nsw = resolved_nonshortwave_surface_flux_components(c,forc,c%T0,H_NONE)
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
    abl = bare_ice_ablation_mass(c,forc,dt_seconds)

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
    abl = bare_ice_ablation_mass(c,forc,dt_seconds)
    call check("dry air -> negative latent flux", abl%latent_heat_flux .lt. 0.0_wp, nfail)
    call check("dry air -> vapor_mass < 0 (mass loss)", abl%vapor_mass .lt. 0.0_wp, nfail)
    call check_close("sublimation_mass = -vapor_mass", &
                     abl%sublimation_mass, -abl%vapor_mass, 1.0e-6_wp, nfail)

    ! Sign convention, saturated warm air: ea > es -> deposition -> mass gain,
    ! and sublimation_mass is clipped at zero.
    forc%relative_humidity = 1.0_wp
    abl = bare_ice_ablation_mass(c,forc,dt_seconds)
    call check("saturated warm air -> positive latent flux", &
               abl%latent_heat_flux .gt. 0.0_wp, nfail)
    call check("deposition -> vapor_mass > 0 (mass gain)", abl%vapor_mass .gt. 0.0_wp, nfail)
    call check_close("sublimation_mass = 0 under deposition", &
                     abl%sublimation_mass, 0.0_wp, 1.0e-6_wp, nfail)

    ! Melt is floored at zero when the net flux is negative.
    call forcing_init(forc)
    forc%air_temperature = 240.0_wp
    forc%shortwave_down  = 0.0_wp
    abl = bare_ice_ablation_mass(c,forc,dt_seconds)
    call check_close("melt_mass = 0 when Q_net < 0", abl%melt_mass, 0.0_wp, 1.0e-6_wp, nfail)

    ! Prescribed albedo is used only when the scheme is PRESCRIBED *and* the
    ! forcing carries one; otherwise bare ice uses alpha_ice.
    call forcing_init(forc)
    forc%air_temperature        = 275.0_wp
    forc%shortwave_down         = 400.0_wp
    forc%prescribed_albedo      = 0.90_wp
    forc%has_prescribed_albedo  = .TRUE.

    abl = bare_ice_ablation_mass(c,forc,dt_seconds)                ! scheme = dynamic
    bif = resolved_bare_ice_surface_flux_components(c,forc,c%alpha_ice)
    q_net = bif%shortwave_absorbed + bif%longwave + bif%sensible + bif%latent + bif%rain
    call check_close("has_prescribed_albedo ignored unless scheme is PRESCRIBED", &
                     abl%melt_mass*c%Lm, max(q_net,0.0_wp)*dt_seconds, 1.0e-5_wp, nfail)

    c%albedo_scheme = CHION_ALBEDO_PRESCRIBED
    abl = bare_ice_ablation_mass(c,forc,dt_seconds)
    bif = resolved_bare_ice_surface_flux_components(c,forc,0.90_wp)
    q_net = bif%shortwave_absorbed + bif%longwave + bif%sensible + bif%latent + bif%rain
    call check_close("scheme PRESCRIBED + flag -> prescribed albedo used", &
                     abl%melt_mass*c%Lm, max(q_net,0.0_wp)*dt_seconds, 1.0e-5_wp, nfail)
    c%albedo_scheme = 2                                            ! back to dynamic

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
    write(*,"(a)") " NOTE: apply_snow_surface_vapor_mass_flux is NOT tested."
    write(*,"(a)") "       It depends on snow_layers (WP4), which does not"
    write(*,"(a)") "       exist yet. See the WP4 INTEGRATION POINT markers."
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
