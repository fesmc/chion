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
    public :: latent_exchange_coefficient
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

    pure function latent_exchange_coefficient(c) result(D_lf)
        ! Chion.jl/src/processes/energy_flux.jl:43-44:
        !     D_lf = latent_heat_flux_ratio * D_sh/cp_air * 0.622 * (Lv + Lm)
        ! The 0.622 is the dry-air/water-vapor molar mass ratio, hard-coded.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp) :: D_lf

        D_lf = c%latent_heat_flux_ratio*c%D_sh/c%cp_air*0.622_wp*(c%Lv + c%Lm)

        return

    end function latent_exchange_coefficient

    pure function latent_vapor_flux(surface_temperature,c,air_temperature, &
                                    relative_humidity,air_pressure) result(q_lh)
        ! EXACT turbulent latent heat flux at a known surface temperature.
        ! Chion.jl/src/processes/energy_flux.jl:73-78.
        !
        ! Positive = flux into the surface (deposition).

        implicit none

        real(wp),                intent(IN) :: surface_temperature   ! [K]
        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature       ! [K]
        real(wp),                intent(IN) :: relative_humidity     ! [1] or [%]
        real(wp),                intent(IN) :: air_pressure          ! [Pa]
        real(wp) :: q_lh                                             ! [W m-2]

        ! Local variables
        real(wp) :: exchange, ea, es

        exchange = latent_exchange_coefficient(c)/safe_positive(air_pressure)
        ea       = air_vapor_pressure(air_temperature,relative_humidity,c%T0)
        es       = ice_saturation_vapor_pressure(surface_temperature,c%T0)

        q_lh = exchange*(ea - es)

        return

    end function latent_vapor_flux

    pure function latent_vapor_flux_linearized(surface_temperature,c,air_temperature, &
                                               relative_humidity,air_pressure) result(coef)
        ! Latent heat flux linearized about surface_temperature = T^n:
        !     Q(T) = coef%constant - coef%linear*T
        ! Chion.jl/src/processes/energy_flux.jl:80-88, which returns
        ! (constant, linear) in that order.
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
        real(wp) :: exchange, ea, es, des_dT

        exchange = latent_exchange_coefficient(c)/safe_positive(air_pressure)
        ea       = air_vapor_pressure(air_temperature,relative_humidity,c%T0)
        es       = ice_saturation_vapor_pressure(surface_temperature,c%T0)
        des_dT   = ice_saturation_vapor_pressure_derivative(surface_temperature,c%T0,es)

        coef%linear   = exchange*des_dT
        coef%constant = exchange*(ea - es + des_dT*surface_temperature)

        return

    end function latent_vapor_flux_linearized

end module snow_vapor
