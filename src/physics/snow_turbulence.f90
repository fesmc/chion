module snow_turbulence
    ! Chion.jl's bulk turbulent exchange, turbulent_flux_scheme = "semix".
    !
    ! Port of Chion.jl/src/processes/energy_flux.jl (9ec6cc7 = 03bb445):
    !   _semix_aerodynamic_resistance     -> turb_semix_resistance
    !   _semix_air_density                -> turb_semix_air_density
    !   _semix_turbulent_flux_linearized  -> turb_semix_flux_linearized
    !
    ! Chion.jl calls this scheme :semix, but it is NOT CLIMBER-X SEMIX's
    ! exchange (that is turbulent_flux_scheme = "climberx", snow_seb_semix;
    ! docs/porting_notes.md D37). It is a bulk scheme of its own, calibrated
    ! against MAR sensible heat over Greenland (03bb445):
    !
    !     C_hn  = k^2/(ln(z/z0m) ln(z/z0h)),   z0h = z0m/zm_to_zh
    !     Ri_b  = g z (T_a - T_s)/(T_a u^2),   u = max(wind, 0.1)
    !     f(Ri) = sqrt(max(1 - 16 Ri, 1))      Ri < 0 (unstable)
    !           = 1/(1 + b Ri)                  Ri >= 0 (stable)
    !     r_a   = 1/(C_hn f(Ri) u)
    !     Q_sh  = a_sh rho_a cp (T_a - T_s)/r_a
    !     Q_lh  = a_lh L rho_a (q_a - q_s(T_s))/r_a
    !
    ! with z = semix_surface_height (10 m), z0m = semix_z0m_snow (0.001 m) or
    ! semix_z0m_ice (0.01 m) on a bare-ice surface, zm_to_zh = 10, the
    ! sensible exchange factor a_sh = 2.5 and the stable coefficient b = 40
    ! (both calibrated), a_lh = 1. The air humidity is relative to saturation
    ! over ICE at the air temperature, q_a = 0.622 rh e_i(T_a)/p, and the
    ! surface humidity q_s = 0.622 e_i(T_s)/p, both from BESSI's ice
    ! saturation vapour pressure (snow_vapor). Both fluxes are linearized
    ! about the surface temperature they are built at, Q(T) = constant -
    ! linear*T, the form the energy solve and the resolved fluxes take.
    !
    ! Constants: g is standard gravity 9.80665, a literal in Julia and
    ! chion's DEF_GRAVITY (D25), not the host's g. The dry-air gas constant is
    ! c%R_dry (287.058); Julia writes the literal 287.05, which legacy_chion
    ! restores (TURB_SEMIX_R_AIR_LITERAL, docs/porting_notes.md D38).
    !
    ! Latent heat: the caller passes it, turb_semix_latent_heat chooses it.
    ! On snow it is the phase's at T_s (Lv+Lm below T0, Lv at it), as in
    ! Julia. On a bare-ice surface (no snow; the substrate's top layer, or ice
    ! held at T0) chion uses Lv+Lm whatever T_s: the ice is solid, and the
    ! vapour MASS is converted with Lv+Lm in both Julia and chion. Julia
    ! carries Lv at T0 there, so its bare-ice energy and mass budgets differ by
    ! Lm per kg sublimated; legacy_chion restores it
    ! (TURB_SEMIX_ICE_SUBLIMATION, docs/porting_notes.md D35).
    !
    ! Evaluated in Julia's order and grouping.

    use chion_defs, only : wp, DEF_GRAVITY, chion_const_class, &
                           TURB_SEMIX_R_AIR_LITERAL, TURB_SEMIX_ICE_SUBLIMATION
    use snow_vapor, only : safe_positive, relative_humidity_fraction, &
                           ice_saturation_vapor_pressure, &
                           ice_saturation_vapor_pressure_derivative, &
                           surface_vapor_latent_heat

    implicit none

    private

    ! Chion.jl's literal dry-air gas constant in _semix_air_density (legacy).
    real(wp), parameter :: TURB_SEMIX_R_AIR_JULIA = 287.05_wp   ! [J kg-1 K-1]

    ! Wind floor of _semix_aerodynamic_resistance.
    real(wp), parameter :: TURB_SEMIX_WIND_MIN = 0.1_wp          ! [m s-1]

    ! Dry-air / water-vapour molar mass ratio, as Julia spells it.
    real(wp), parameter :: TURB_SEMIX_EPS_VAPOR = 0.622_wp

    type turb_semix_lin_class
        ! _semix_turbulent_flux_linearized returns (sensible_constant,
        ! sensible_coefficient, latent_constant, latent_linear) IN THAT ORDER,
        ! each pair such that Q(T) = constant - linear*T.
        real(wp) :: sensible_constant  ! [W m-2]
        real(wp) :: sensible_linear    ! [W m-2 K-1]
        real(wp) :: latent_constant    ! [W m-2]
        real(wp) :: latent_linear      ! [W m-2 K-1]
    end type turb_semix_lin_class

    public :: turb_semix_lin_class
    public :: turb_semix_resistance
    public :: turb_semix_air_density
    public :: turb_semix_latent_heat
    public :: turb_semix_roughness
    public :: turb_semix_flux_linearized

contains

    pure function turb_semix_resistance(c,surface_temperature,air_temperature, &
                                        wind_speed,z0m) result(r_a)
        ! energy_flux.jl _semix_aerodynamic_resistance. Positive Ri is warm
        ! air over a colder surface (stable): exchange damped by 1/(1 + b Ri).
        ! Negative Ri is unstable: enhanced by sqrt(1 - 16 Ri).

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: surface_temperature   ! [K]
        real(wp),                intent(IN) :: air_temperature       ! [K]
        real(wp),                intent(IN) :: wind_speed            ! [m s-1]
        real(wp),                intent(IN) :: z0m                   ! [m]
        real(wp) :: r_a                                              ! [s m-1]

        ! Local variables
        real(wp) :: z0h, neutral_ch, wind, bulk_richardson, stability_factor

        z0h        = z0m/c%semix_zm_to_zh
        neutral_ch = (c%semix_karman*c%semix_karman) &
                     /safe_positive(log(c%semix_surface_height/z0m) &
                                    *log(c%semix_surface_height/z0h))

        wind = max(wind_speed,TURB_SEMIX_WIND_MIN)

        bulk_richardson = DEF_GRAVITY*c%semix_surface_height &
                          *(air_temperature - surface_temperature) &
                          /safe_positive(air_temperature*(wind*wind))

        if (bulk_richardson .lt. 0.0_wp) then
            stability_factor = sqrt(max(1.0_wp - 16.0_wp*bulk_richardson,1.0_wp))
        else
            stability_factor = 1.0_wp/(1.0_wp + c%semix_stable_coefficient*bulk_richardson)
        end if

        r_a = 1.0_wp/safe_positive(neutral_ch*stability_factor*wind)

        return

    end function turb_semix_resistance

    pure function turb_semix_air_density(c,air_temperature,air_pressure) result(rho_a)
        ! energy_flux.jl _semix_air_density: p/(R T_a), R = c%R_dry
        ! (Julia's literal 287.05 under legacy_chion, D38).

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature   ! [K]
        real(wp),                intent(IN) :: air_pressure      ! [Pa]
        real(wp) :: rho_a                                        ! [kg m-3]

        ! Local variables
        real(wp) :: R_air

        if (TURB_SEMIX_R_AIR_LITERAL) then
            R_air = TURB_SEMIX_R_AIR_JULIA
        else
            R_air = c%R_dry
        end if

        rho_a = air_pressure/safe_positive(R_air*air_temperature)

        return

    end function turb_semix_air_density

    pure function turb_semix_latent_heat(c,surface_temperature,surface_is_ice) result(L)
        ! The latent heat of the surface's vapour exchange: the phase's at
        ! surface_temperature on snow (energy_flux.jl
        ! _surface_vapor_latent_heat); Lv+Lm on bare ice, which is solid at any
        ! temperature (chion, D35; legacy_chion: the snow rule, as Julia).

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: surface_temperature   ! [K]
        logical,                 intent(IN) :: surface_is_ice
        real(wp) :: L                                                ! [J kg-1]

        if (surface_is_ice .and. TURB_SEMIX_ICE_SUBLIMATION) then
            L = c%Lv + c%Lm
        else
            L = surface_vapor_latent_heat(surface_temperature,c)
        end if

        return

    end function turb_semix_latent_heat

    pure function turb_semix_roughness(c,surface_is_ice) result(z0m)
        ! Momentum roughness of the surface: semix_z0m_ice on bare ice (and
        ! the substrate's top layer), semix_z0m_snow otherwise. No snow-depth
        ! blend, unlike CLIMBER-X.

        implicit none

        type(chion_const_class), intent(IN) :: c
        logical,                 intent(IN) :: surface_is_ice
        real(wp) :: z0m                                       ! [m]

        if (surface_is_ice) then
            z0m = c%semix_z0m_ice
        else
            z0m = c%semix_z0m_snow
        end if

        return

    end function turb_semix_roughness

    pure function turb_semix_flux_linearized(c,surface_temperature,air_temperature, &
                                             relative_humidity,air_pressure,wind_speed, &
                                             z0m,latent_heat) result(x)
        ! energy_flux.jl _semix_turbulent_flux_linearized, with the latent heat
        ! passed in (turb_semix_latent_heat) instead of taken at
        ! surface_temperature. Both fluxes linearized about
        ! surface_temperature:
        !     Q_sh(T) = sensible_constant - sensible_linear*T
        !     Q_lh(T) = latent_constant   - latent_linear*T
        ! Exact at T = surface_temperature. The caller zeroes the latent flux
        ! without humidity forcing, as Julia does.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: surface_temperature   ! [K]
        real(wp),                intent(IN) :: air_temperature       ! [K]
        real(wp),                intent(IN) :: relative_humidity     ! [1] or [%]
        real(wp),                intent(IN) :: air_pressure          ! [Pa]
        real(wp),                intent(IN) :: wind_speed            ! [m s-1]
        real(wp),                intent(IN) :: z0m                   ! [m]
        real(wp),                intent(IN) :: latent_heat           ! [J kg-1]
        type(turb_semix_lin_class) :: x

        ! Local variables
        real(wp) :: air_density, resistance, sensible_coefficient
        real(wp) :: q_air, q_surface_pressure, q_surface, dq_surface
        real(wp) :: latent_exchange

        air_density = turb_semix_air_density(c,air_temperature,air_pressure)
        resistance  = turb_semix_resistance(c,surface_temperature,air_temperature, &
                                            wind_speed,z0m)

        sensible_coefficient = c%semix_sensible_exchange_factor*air_density*c%cp_air &
                               /resistance
        x%sensible_constant  = sensible_coefficient*air_temperature
        x%sensible_linear    = sensible_coefficient

        q_air = TURB_SEMIX_EPS_VAPOR*relative_humidity_fraction(relative_humidity) &
                *ice_saturation_vapor_pressure(air_temperature,c%T0) &
                /safe_positive(air_pressure)
        q_surface_pressure = ice_saturation_vapor_pressure(surface_temperature,c%T0)
        q_surface  = TURB_SEMIX_EPS_VAPOR*q_surface_pressure/safe_positive(air_pressure)
        dq_surface = TURB_SEMIX_EPS_VAPOR &
                     *ice_saturation_vapor_pressure_derivative(surface_temperature,c%T0, &
                                                               q_surface_pressure) &
                     /safe_positive(air_pressure)

        latent_exchange   = c%semix_latent_exchange_factor*latent_heat*air_density/resistance
        x%latent_constant = latent_exchange*(q_air - q_surface + dq_surface*surface_temperature)
        x%latent_linear   = latent_exchange*dq_surface

        return

    end function turb_semix_flux_linearized

end module snow_turbulence
