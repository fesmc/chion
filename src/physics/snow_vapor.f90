module snow_vapor
    ! Vapour-pressure parameterizations and the BESSI turbulent latent flux
    ! built from them.
    !
    ! Port of Chion.jl/src/processes/energy_flux.jl:40-95. These live in their
    ! own module because they are the lowest layer of the surface physics:
    ! snow_surface_fluxes (exact fluxes at a known temperature), snow_energy
    ! (the linearized surface row) and snow_seb_semix (the SEMIX aerodynamic
    ! scheme) all sit on top of them, and snow_seb_semix must compile BELOW
    ! snow_surface_fluxes so the latter can dispatch on seb_scheme.
    !
    ! All coefficients here are magic numbers carried over verbatim from
    ! Chion.jl; none of them live in the constants struct. See docs/PLAN.md
    ! section 5, item 3.

    use chion_defs, only : wp, wp_acc, TOL_TINY, chion_const_class

    implicit none

    private

    type latent_vapor_flux_lin_class
        ! _bessi_latent_vapor_flux_linearized (energy_flux.jl:80) returns
        ! (constant, linear) IN THAT ORDER, such that
        !     Q_latent(T) = constant - linear*T
        real(wp) :: constant           ! [W m-2]
        real(wp) :: linear             ! [W m-2 K-1]
    end type latent_vapor_flux_lin_class

    public :: latent_vapor_flux_lin_class

    public :: safe_positive
    public :: relative_humidity_fraction
    public :: water_saturation_vapor_pressure
    public :: air_vapor_pressure
    public :: ice_saturation_vapor_pressure
    public :: ice_saturation_vapor_pressure_derivative
    public :: vapor_exchange_coefficient
    public :: surface_vapor_latent_heat
    public :: vapor_mass_flux
    public :: latent_vapor_flux
    public :: latent_vapor_flux_linearized

contains

    pure function safe_positive(x) result(y)
        ! Chion.jl/src/processes/energy_flux.jl:5.
        ! Floor at EPS_TINY, used to protect divisions. Note TOL_TINY is
        ! declared dp (docs/porting_notes.md D1b); the comparison promotes and
        ! the floor is converted back to wp, which is exact for 1e-12.

        implicit none

        real(wp), intent(IN) :: x
        real(wp) :: y

        if (real(x,wp_acc) .gt. TOL_TINY) then
            y = x
        else
            y = real(TOL_TINY,wp)
        end if

        return

    end function safe_positive

    pure function relative_humidity_fraction(relative_humidity) result(rh)
        ! Chion.jl/src/processes/energy_flux.jl:46-49.
        ! Values greater than 1 are interpreted as percent, then clamped to
        ! [0,1]. Note a genuine fraction of exactly 1.0 is NOT rescaled.

        implicit none

        real(wp), intent(IN) :: relative_humidity   ! [1] or [%]
        real(wp) :: rh

        if (relative_humidity .gt. 1.0_wp) then
            rh = relative_humidity/100.0_wp
        else
            rh = relative_humidity
        end if

        rh = min(max(rh,0.0_wp),1.0_wp)

        return

    end function relative_humidity_fraction

    pure function water_saturation_vapor_pressure(temperature,T0) result(es)
        ! Saturation vapor pressure over WATER, used for the air.
        ! Chion.jl/src/processes/energy_flux.jl:51-55. Magic: 611.2/17.27/243.12.

        implicit none

        real(wp), intent(IN) :: temperature      ! [K]
        real(wp), intent(IN) :: T0               ! [K] freezing point
        real(wp) :: es                           ! [Pa]

        ! Local variables
        real(wp) :: tc

        tc = temperature - T0
        es = 611.2_wp*exp(17.27_wp*tc/(tc + 243.12_wp))

        return

    end function water_saturation_vapor_pressure

    pure function air_vapor_pressure(air_temperature,relative_humidity,T0) result(ea)
        ! Chion.jl/src/processes/energy_flux.jl:57-60.

        implicit none

        real(wp), intent(IN) :: air_temperature    ! [K]
        real(wp), intent(IN) :: relative_humidity  ! [1] or [%]
        real(wp), intent(IN) :: T0                 ! [K]
        real(wp) :: ea                             ! [Pa]

        ea = relative_humidity_fraction(relative_humidity) &
             *water_saturation_vapor_pressure(air_temperature,T0)

        return

    end function air_vapor_pressure

    pure function ice_saturation_vapor_pressure(surface_temperature,T0) result(es)
        ! Saturation vapor pressure over ICE, used for the surface.
        ! Chion.jl/src/processes/energy_flux.jl:62-66. Magic: 611.2/22.46/272.62.

        implicit none

        real(wp), intent(IN) :: surface_temperature   ! [K]
        real(wp), intent(IN) :: T0                    ! [K]
        real(wp) :: es                                ! [Pa]

        ! Local variables
        real(wp) :: tc

        tc = surface_temperature - T0
        es = 611.2_wp*exp(22.46_wp*tc/(tc + 272.62_wp))

        return

    end function ice_saturation_vapor_pressure

    pure function ice_saturation_vapor_pressure_derivative(surface_temperature,T0,es) &
                                                                    result(des_dT)
        ! d(es_ice)/dT. Chion.jl/src/processes/energy_flux.jl:68-71.
        ! es is passed in rather than recomputed, exactly as in Julia.

        implicit none

        real(wp), intent(IN) :: surface_temperature   ! [K]
        real(wp), intent(IN) :: T0                    ! [K]
        real(wp), intent(IN) :: es                    ! [Pa]
        real(wp) :: des_dT                            ! [Pa K-1]

        ! Local variables
        real(wp) :: denom

        denom  = surface_temperature - T0 + 272.62_wp
        des_dT = es*22.46_wp*272.62_wp/safe_positive(denom*denom)

        return

    end function ice_saturation_vapor_pressure_derivative

    pure function vapor_exchange_coefficient(c) result(D_v)
        ! Chion.jl/src/processes/energy_flux.jl:62-63 (dev_nils d0146e1):
        !     D_v = latent_heat_flux_ratio * D_sh/cp_air * 0.622
        ! the bulk vapour-MASS transfer coefficient; divided by the air
        ! pressure it turns a vapour-pressure difference into a mass flux. The
        ! 0.622 is the dry-air/water-vapor molar mass ratio, hard-coded.
        !
        ! It no longer carries a latent heat: the energy flux multiplies the
        ! mass flux by the phase's latent heat (surface_vapor_latent_heat),
        ! while the vapour mass is the mass flux itself.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp) :: D_v                 ! [kg m-2 s-1 (Pa/Pa)]

        D_v = c%latent_heat_flux_ratio*c%D_sh/c%cp_air*0.622_wp

        return

    end function vapor_exchange_coefficient

    pure function surface_vapor_latent_heat(surface_temperature,c) result(L)
        ! Chion.jl/src/processes/energy_flux.jl:67-68 (d0146e1): vapour
        ! exchange with a subfreezing (solid) surface is sublimation/deposition
        ! and carries Lv + Lm; at the melting point it is evaporation/
        ! condensation of surface water and carries Lv.

        implicit none

        real(wp),                intent(IN) :: surface_temperature   ! [K]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: L                                                ! [J kg-1]

        if (surface_temperature .lt. c%T0) then
            L = c%Lv + c%Lm
        else
            L = c%Lv
        end if

        return

    end function surface_vapor_latent_heat

    pure function vapor_mass_flux(surface_temperature,c,air_temperature, &
                                  relative_humidity,air_pressure) result(E)
        ! BESSI turbulent vapour-mass flux from the vapour-pressure gradient,
        !     E = D_v/p * (e_a - e_s(T_s)),
        ! Chion.jl/src/processes/energy_flux.jl:97-102 (d0146e1).
        !
        ! Positive = flux into the surface (deposition).

        implicit none

        real(wp),                intent(IN) :: surface_temperature   ! [K]
        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature       ! [K]
        real(wp),                intent(IN) :: relative_humidity     ! [1] or [%]
        real(wp),                intent(IN) :: air_pressure          ! [Pa]
        real(wp) :: E                                                ! [kg m-2 s-1]

        ! Local variables
        real(wp) :: exchange, ea, es

        exchange = vapor_exchange_coefficient(c)/safe_positive(air_pressure)
        ea       = air_vapor_pressure(air_temperature,relative_humidity,c%T0)
        es       = ice_saturation_vapor_pressure(surface_temperature,c%T0)

        E = exchange*(ea - es)

        return

    end function vapor_mass_flux

    pure function latent_vapor_flux(surface_temperature,c,air_temperature, &
                                    relative_humidity,air_pressure,latent_heat) result(q_lh)
        ! EXACT turbulent latent heat flux at a known surface temperature,
        ! Chion.jl/src/processes/energy_flux.jl:104-107: the vapour-mass flux
        ! times the latent heat of the exchange. Julia defaults latent_heat to
        ! surface_vapor_latent_heat(surface_temperature); here every caller
        ! passes it, since bare ice (solid at T0) uses Lv + Lm.
        !
        ! Positive = flux into the surface (deposition).

        implicit none

        real(wp),                intent(IN) :: surface_temperature   ! [K]
        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature       ! [K]
        real(wp),                intent(IN) :: relative_humidity     ! [1] or [%]
        real(wp),                intent(IN) :: air_pressure          ! [Pa]
        real(wp),                intent(IN) :: latent_heat           ! [J kg-1]
        real(wp) :: q_lh                                             ! [W m-2]

        q_lh = latent_heat*vapor_mass_flux(surface_temperature,c,air_temperature, &
                                           relative_humidity,air_pressure)

        return

    end function latent_vapor_flux

    pure function latent_vapor_flux_linearized(surface_temperature,c,air_temperature, &
                                               relative_humidity,air_pressure) result(coef)
        ! Latent heat flux linearized about surface_temperature = T^n:
        !     Q(T) = coef%constant - coef%linear*T
        ! Chion.jl/src/processes/energy_flux.jl:109-119, which returns
        ! (constant, linear) in that order. The latent heat is that of the
        ! phase at T^n (surface_vapor_latent_heat), held fixed over the step.
        !
        ! At T = T^n this reproduces latent_vapor_flux exactly; away from T^n
        ! it does not, and that is the deliberate inconsistency of trap 2
        ! (see the banner in snow_surface_fluxes).

        implicit none

        real(wp),                intent(IN) :: surface_temperature   ! [K] = T^n
        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature       ! [K]
        real(wp),                intent(IN) :: relative_humidity     ! [1] or [%]
        real(wp),                intent(IN) :: air_pressure          ! [Pa]
        type(latent_vapor_flux_lin_class) :: coef

        ! Local variables
        real(wp) :: exchange, ea, es, des_dT, latent_heat

        exchange    = vapor_exchange_coefficient(c)/safe_positive(air_pressure)
        ea          = air_vapor_pressure(air_temperature,relative_humidity,c%T0)
        es          = ice_saturation_vapor_pressure(surface_temperature,c%T0)
        des_dT      = ice_saturation_vapor_pressure_derivative(surface_temperature,c%T0,es)
        latent_heat = surface_vapor_latent_heat(surface_temperature,c)

        coef%linear   = latent_heat*exchange*des_dT
        coef%constant = latent_heat*exchange*(ea - es + des_dT*surface_temperature)

        return

    end function latent_vapor_flux_linearized

end module snow_vapor
