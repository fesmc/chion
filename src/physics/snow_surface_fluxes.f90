module snow_surface_fluxes
    ! Surface energy and vapor-mass fluxes evaluated at a KNOWN surface
    ! temperature, plus the vapor-pressure parameterizations shared with the
    ! energy solver.
    !
    ! Port of Chion.jl/src/processes/surface_fluxes.jl, together with the
    ! vapor-pressure helpers at Chion.jl/src/processes/energy_flux.jl:40-95 and
    ! _diagnose_latent_heat_flux_coefficients (energy_flux.jl:277-296).
    !
    ! ----------------------------------------------------------------------
    ! TRAP 2 -- TWO INCONSISTENT SURFACE-FLUX EVALUATIONS, BY DESIGN
    ! ----------------------------------------------------------------------
    ! Everything in this module evaluates the surface fluxes EXACTLY, at a
    ! surface temperature that is already known:
    !     * at c%T0            for bare ice (bare ice is ASSUMED to be at the
    !                          melting point; it is never solved for)
    !     * at T^{n+1}         for the post-solve vapor mass flux
    ! snow_energy instead uses fluxes LINEARIZED about T^n. The energy budget
    ! and the mass budget therefore see slightly different latent heat fluxes.
    ! This is deliberate in Chion.jl and is preserved here. It is listed in
    ! docs/PLAN.md section 4.1 as "not allowed without asking" -- do NOT
    ! reconcile the two. See docs/PLAN.md section 5, item 2.
    !
    ! A second, smaller inconsistency is preserved for the same reason: the
    ! shortwave absorption HERE applies max(shortwave_down,0)
    ! (surface_fluxes.jl:74), while the one inside the energy solve does NOT
    ! (energy_flux.jl:99). Both are reproduced as written.
    ! ----------------------------------------------------------------------
    !
    ! CALLING CONVENTION: routines that touch the layered state take contiguous
    ! column slices, e.g. mass(:,icol), plus the active layer count n. See
    ! docs/porting_notes.md D8.
    !
    ! PRECISION: state and fluxes are wp; sums of energy fluxes that feed a
    ! mass conversion are accumulated in real(wp_acc) locals. See
    ! docs/PLAN.md section 3.1.

    use chion_defs, only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER, io_unit_err, &
                           CHION_ALBEDO_PRESCRIBED, &
                           chion_const_class, chion_step_forcing_class
    use snow_column_utils, only : surface_has_snow

    ! snow_layers supplies the layer removal/merge used by
    ! apply_snow_surface_vapor_mass_flux when sublimation empties the surface
    ! layer. Julia passes typemax(Int) as a max-index argument to
    ! _merge_surface_layer!; the Fortran signature has no such argument, so it
    ! is simply dropped.
    use snow_layers, only : remove_depleted_surface_and_route_water, &
                            merge_surface_layer

    implicit none

    private

    ! === Return types ========================================================
    ! Julia returns named tuples; Fortran gets small derived types. FIELD ORDER
    ! MATTERS and mirrors the Julia tuple order in every case.

    type nonshortwave_flux_class
        ! _resolved_nonshortwave_surface_flux_components (surface_fluxes.jl:20)
        real(wp) :: longwave         ! [W m-2] net longwave
        real(wp) :: sensible         ! [W m-2] sensible heat flux
        real(wp) :: latent           ! [W m-2] turbulent latent heat flux
        real(wp) :: rain             ! [W m-2] sensible heat carried by rain
    end type nonshortwave_flux_class

    type bare_ice_flux_class
        ! _resolved_bare_ice_surface_flux_components (surface_fluxes.jl:53)
        real(wp) :: shortwave_absorbed  ! [W m-2]
        real(wp) :: longwave            ! [W m-2]
        real(wp) :: sensible            ! [W m-2]
        real(wp) :: latent              ! [W m-2]
        real(wp) :: rain                ! [W m-2]
    end type bare_ice_flux_class

    type bare_ice_ablation_class
        ! _bare_ice_surface_mass_fluxes_resolved (surface_fluxes.jl:104).
        !
        ! SIGN CONVENTION (BESSI): POSITIVE vapor_mass = deposition, i.e. a
        ! mass GAIN. Negative vapor_mass is sublimation, and sublimation_mass
        ! is its positive-definite mirror, max(-vapor_mass,0).
        real(wp) :: melt_mass          ! [kg m-2] >= 0
        real(wp) :: vapor_mass         ! [kg m-2] signed; + = deposition
        real(wp) :: sublimation_mass   ! [kg m-2] >= 0
        real(wp) :: latent_heat_flux   ! [W m-2] signed
        real(wp) :: net_mass_change    ! [kg m-2] vapor_mass - melt_mass
    end type bare_ice_ablation_class

    type surface_vapor_flux_class
        ! _surface_vapor_fluxes (surface_fluxes.jl:5)
        real(wp) :: vapor_mass         ! [kg m-2] signed; + = deposition
        real(wp) :: sublimation_mass   ! [kg m-2] >= 0
        real(wp) :: latent_heat_flux   ! [W m-2] signed
    end type surface_vapor_flux_class

    type latent_vapor_flux_lin_class
        ! _bessi_latent_vapor_flux_linearized (energy_flux.jl:80) returns
        ! (constant, linear) IN THAT ORDER, such that
        !     Q_latent(T) = constant - linear*T
        real(wp) :: constant           ! [W m-2]
        real(wp) :: linear             ! [W m-2 K-1]
    end type latent_vapor_flux_lin_class

    type latent_heat_coeff_class
        ! _diagnose_latent_heat_flux_coefficients (energy_flux.jl:282) returns
        ! (linear, constant) IN THAT ORDER -- note this is the OPPOSITE order
        ! to latent_vapor_flux_lin_class above. Both orders are kept as the
        ! Julia source has them, so the two are separate types on purpose.
        real(wp) :: linear             ! [W m-2 K-1]
        real(wp) :: constant           ! [W m-2]
    end type latent_heat_coeff_class

    public :: nonshortwave_flux_class
    public :: bare_ice_flux_class
    public :: bare_ice_ablation_class
    public :: surface_vapor_flux_class
    public :: latent_vapor_flux_lin_class
    public :: latent_heat_coeff_class

    ! Vapor-pressure parameterizations (also used by snow_energy)
    public :: safe_positive
    public :: relative_humidity_fraction
    public :: water_saturation_vapor_pressure
    public :: air_vapor_pressure
    public :: ice_saturation_vapor_pressure
    public :: ice_saturation_vapor_pressure_derivative
    public :: latent_exchange_coefficient
    public :: latent_vapor_flux
    public :: latent_vapor_flux_linearized

    ! Resolved (exact-at-known-T) surface fluxes
    public :: resolved_nonshortwave_surface_flux_components
    public :: resolved_bare_ice_surface_flux_components
    public :: resolved_turbulent_latent_heat_flux
    public :: bare_ice_ablation_mass

    ! Coefficients and post-solve vapor mass bookkeeping
    public :: diagnose_latent_heat_flux_coefficients
    public :: apply_snow_surface_vapor_mass_flux

contains

    ! =====================================================================
    ! Vapor-pressure parameterizations
    !
    ! All coefficients are magic numbers carried over verbatim from
    ! Chion.jl; none of them live in the constants struct. See
    ! docs/PLAN.md section 5, item 3.
    ! =====================================================================

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
        ! it does not, and that is the deliberate inconsistency of trap 2.

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

    ! =====================================================================
    ! Resolved (exact-at-known-T) surface fluxes
    ! =====================================================================

    pure function resolved_nonshortwave_surface_flux_components(c,forc,surface_temperature) &
                                                                        result(flx)
        ! Chion.jl/src/processes/surface_fluxes.jl:20-45.
        !
        ! Every has_* flag selects a PRESCRIBED value over the internal
        ! parameterization. Note the latent flux has three cases, in order:
        !   has_q_lh                 -> prescribed q_lh
        !   has_relative_humidity    -> BESSI vapor flux
        !   otherwise                -> zero
        !
        ! Julia carries a dt_seconds argument here purely to spell zero() in
        ! the right type; it is never used numerically. Dropped (cleanup:
        ! Julia type-dispatch artifact, no behaviour change).

        implicit none

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),                       intent(IN) :: surface_temperature  ! [K]
        type(nonshortwave_flux_class) :: flx

        if (forc%has_q_lw_down) then
            flx%longwave = forc%q_lw_down &
                           - c%sigma_sb*c%eps_snow*surface_temperature**4
        else
            flx%longwave = c%sigma_sb*(c%eps_air*forc%air_temperature**4 &
                                       - c%eps_snow*surface_temperature**4)
        end if

        if (forc%has_q_sh) then
            flx%sensible = forc%q_sh
        else
            flx%sensible = c%D_sh*(forc%air_temperature - surface_temperature)
        end if

        if (forc%has_q_lh) then
            flx%latent = forc%q_lh
        else if (forc%has_relative_humidity) then
            flx%latent = latent_vapor_flux(surface_temperature,c,forc%air_temperature, &
                                           forc%relative_humidity,forc%air_pressure)
        else
            flx%latent = 0.0_wp
        end if

        flx%rain = forc%rainfall_rate*c%cw*(forc%air_temperature - c%T0)

        return

    end function resolved_nonshortwave_surface_flux_components

    pure function resolved_bare_ice_surface_flux_components(c,forc,surface_albedo) result(flx)
        ! Chion.jl/src/processes/surface_fluxes.jl:53-95.
        !
        ! IMPORTANT: the non-shortwave components are evaluated at
        ! surface_temperature = c%T0. Bare ice is ASSUMED to sit at the melting
        ! point; its surface temperature is never solved for. Any excess energy
        ! becomes melt rather than warming.
        !
        ! NOTE the max(shortwave_down,0) below. The twin expression inside the
        ! energy solve (energy_flux.jl:99) has NO such clamp. Both are
        ! preserved -- see the trap 2 banner at the top of this module.

        implicit none

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),                       intent(IN) :: surface_albedo   ! [1]
        type(bare_ice_flux_class) :: flx

        ! Local variables
        type(nonshortwave_flux_class) :: nsw

        if (forc%has_q_sw_net) then
            flx%shortwave_absorbed = forc%q_sw_net
        else
            flx%shortwave_absorbed = max(forc%shortwave_down,0.0_wp) &
                                     *(1.0_wp - min(max(surface_albedo,0.0_wp),1.0_wp))
        end if

        nsw = resolved_nonshortwave_surface_flux_components(c,forc,c%T0)

        flx%longwave = nsw%longwave
        flx%sensible = nsw%sensible
        flx%latent   = nsw%latent
        flx%rain     = nsw%rain

        return

    end function resolved_bare_ice_surface_flux_components

    pure function resolved_turbulent_latent_heat_flux(c,forc,surface_temperature) result(q_lh)
        ! Chion.jl/src/processes/surface_fluxes.jl:187-200.
        ! The latent-flux-only subset of the three-case selection above.

        implicit none

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),                       intent(IN) :: surface_temperature   ! [K]
        real(wp) :: q_lh                                                    ! [W m-2]

        if (forc%has_q_lh) then
            q_lh = forc%q_lh
        else if (forc%has_relative_humidity) then
            q_lh = latent_vapor_flux(surface_temperature,c,forc%air_temperature, &
                                     forc%relative_humidity,forc%air_pressure)
        else
            q_lh = 0.0_wp
        end if

        return

    end function resolved_turbulent_latent_heat_flux

    pure function bare_ice_ablation_mass(c,forc,dt_seconds) result(abl)
        ! Chion.jl/src/processes/surface_fluxes.jl:104-185
        ! (_bare_ice_surface_mass_fluxes_resolved + _bare_ice_ablation_mass).
        !
        ! Bare ice stays on the daily-mean surface-energy path: the diurnal
        ! shortwave correction is snow-column-only.
        !
        ! IMPORTANT: melt and the vapor exchange are computed from the SAME
        ! Q_net INDEPENDENTLY. The melt energy is NOT reduced by the energy that
        ! went into sublimation/deposition, and the vapor mass always uses
        ! (Lv + Lm) on bare ice regardless of temperature. Preserve both.
        !
        ! Albedo selection mirrors surface_fluxes.jl:181-183: the prescribed
        ! albedo is used only when the scheme is PRESCRIBED *and* the forcing
        ! actually carries one; otherwise bare ice uses c%alpha_ice.

        implicit none

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),                       intent(IN) :: dt_seconds   ! [s]
        type(bare_ice_ablation_class) :: abl

        ! Local variables
        real(wp)                  :: surface_albedo
        type(bare_ice_flux_class) :: flx
        real(wp_acc)              :: q_net

        surface_albedo = c%alpha_ice
        if (c%albedo_scheme .eq. CHION_ALBEDO_PRESCRIBED) then
            if (forc%has_prescribed_albedo) then
                surface_albedo = min(max(forc%prescribed_albedo,0.0_wp),1.0_wp)
            end if
        end if

        flx = resolved_bare_ice_surface_flux_components(c,forc,surface_albedo)

        ! Accumulated in wp_acc: a sum of five fluxes feeding a mass conversion.
        q_net = real(flx%shortwave_absorbed,wp_acc) + real(flx%longwave,wp_acc) &
              + real(flx%sensible,wp_acc)           + real(flx%latent,wp_acc)   &
              + real(flx%rain,wp_acc)

        abl%melt_mass = real(max(q_net,0.0_wp_acc)*real(dt_seconds,wp_acc) &
                             /real(c%Lm,wp_acc),wp)

        abl%latent_heat_flux = flx%latent
        abl%vapor_mass       = real(real(flx%latent,wp_acc)*real(dt_seconds,wp_acc) &
                                    /real(c%Lv + c%Lm,wp_acc),wp)
        abl%sublimation_mass = max(-abl%vapor_mass,0.0_wp)

        abl%net_mass_change = abl%vapor_mass - abl%melt_mass

        return

    end function bare_ice_ablation_mass

    ! =====================================================================
    ! Latent-heat coefficients from precipitation
    ! =====================================================================

    pure function diagnose_latent_heat_flux_coefficients(has_surface_snow,c,air_temperature, &
                                                         snowfall_rate,rainfall_rate) result(coef)
        ! Chion.jl/src/processes/energy_flux.jl:282-296, returning
        ! (linear, constant) in that order.
        !
        ! These are the effective surface-flux terms contributed by falling
        ! precipitation, added to the linearized surface row of the energy
        ! solve. SNOWFALL TAKES PRECEDENCE over rainfall: if any snow is
        ! falling, the rain term is not applied at all, even when it is also
        ! raining. Rain contributes only when the surface already has snow.

        implicit none

        logical,                 intent(IN) :: has_surface_snow
        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature    ! [K]
        real(wp),                intent(IN) :: snowfall_rate      ! [kg m-2 s-1]
        real(wp),                intent(IN) :: rainfall_rate      ! [kg m-2 s-1]
        type(latent_heat_coeff_class) :: coef

        if (snowfall_rate .gt. 0.0_wp) then
            coef%linear   = snowfall_rate*c%ci
            coef%constant = snowfall_rate*c%ci*air_temperature
        else if (has_surface_snow .and. rainfall_rate .gt. 0.0_wp) then
            coef%linear   = 0.0_wp
            coef%constant = rainfall_rate*c%cw*(air_temperature - c%T0)
        else
            coef%linear   = 0.0_wp
            coef%constant = 0.0_wp
        end if

        return

    end function diagnose_latent_heat_flux_coefficients

    ! =====================================================================
    ! Post-solve surface vapor mass flux
    ! =====================================================================

    subroutine apply_snow_surface_vapor_mass_flux(mass,mass_w,density,temperature,n, &
                                                  runoff,t_srf,albedo,c,forc,dt_seconds, &
                                                  mass_split,mass_min,vflux)
        ! Chion.jl/src/processes/surface_fluxes.jl:210-269.
        !
        ! Called AFTER the implicit energy solve, so temperature(1) is the NEW
        ! surface temperature T^{n+1}. The latent heat flux is therefore
        ! re-evaluated EXACTLY at T^{n+1}, while the energy solve used the
        ! linearization about T^n -- trap 2 again, deliberate.
        !
        ! Two branches, on the NEW surface temperature:
        !   Ts <  T0  solid exchange:  vapor = Q*dt/(Lv+Lm), applied to mass(1)
        !   Ts >= T0  liquid exchange: vapor = Q*dt/Lv,      applied to mass_w(1)
        ! Both take max(...,0) on the updated layer value, so a sublimation
        ! demand larger than the available surface mass is silently truncated
        ! and the vapor_mass returned is NOT reduced to match. Preserved as-is.
        !
        ! Only the solid branch can empty the surface layer, so only it runs
        ! the depleted-surface removal and surface-merge loops.

        implicit none

        real(wp),                       intent(INOUT) :: mass(:)        ! (Ntot) [kg m-2]
        real(wp),                       intent(INOUT) :: mass_w(:)      ! (Ntot) [kg m-2]
        real(wp),                       intent(INOUT) :: density(:)     ! (Ntot) [kg m-3]
        real(wp),                       intent(INOUT) :: temperature(:) ! (Ntot) [K]
        integer,                        intent(INOUT) :: n              ! active layer count
        real(wp_acc),                   intent(INOUT) :: runoff         ! [kg m-2] cumulative
        real(wp),                       intent(INOUT) :: t_srf          ! [K]
        real(wp),                       intent(INOUT) :: albedo         ! [1] dynamic albedo
        type(chion_const_class),        intent(IN)    :: c
        type(chion_step_forcing_class), intent(IN)    :: forc
        real(wp),                       intent(IN)    :: dt_seconds     ! [s]
        real(wp),                       intent(IN)    :: mass_split     ! [kg m-2]
        real(wp),                       intent(IN)    :: mass_min       ! [kg m-2]
        type(surface_vapor_flux_class), intent(OUT)   :: vflux

        ! Local variables
        real(wp)     :: surface_temperature, q_lh
        real(wp_acc) :: vapor

        vflux%vapor_mass       = 0.0_wp
        vflux%sublimation_mass = 0.0_wp
        vflux%latent_heat_flux = 0.0_wp

        if (.not. surface_has_snow(mass,n)) return

        surface_temperature = temperature(1)

        q_lh = resolved_turbulent_latent_heat_flux(c,forc,surface_temperature)

        vflux%latent_heat_flux = q_lh

        ! Julia returns early on an exactly-zero flux, so that a column with no
        ! humidity forcing never touches the layer structure. Reproduced.
        if (q_lh .eq. 0.0_wp) return

        if (surface_temperature .lt. c%T0) then
            vapor = real(q_lh,wp_acc)*real(dt_seconds,wp_acc)/real(c%Lv + c%Lm,wp_acc)
        else
            vapor = real(q_lh,wp_acc)*real(dt_seconds,wp_acc)/real(c%Lv,wp_acc)
        end if

        vflux%vapor_mass       = real(vapor,wp)
        vflux%sublimation_mass = max(-vflux%vapor_mass,0.0_wp)

        if (surface_temperature .lt. c%T0) then

            ! Solid exchange: sublimation/deposition of the surface snow layer.
            mass(1) = max(mass(1) + vflux%vapor_mass,0.0_wp)

            ! Peel off any surface layer that sublimation has emptied, routing
            ! its liquid water to runoff.
            do while (n .gt. 0)
                if (mass(1) .gt. real(TOL_EMPTY_LAYER,wp)) exit
                call remove_depleted_surface_and_route_water(mass,mass_w,density, &
                                                             temperature,n,runoff,c)
            end do

            ! Re-merge a surface layer that has become too thin.
            do while (n .gt. 1)
                if (mass(1) .ge. mass_min) exit
                call merge_surface_layer(mass,mass_w,density,temperature,n, &
                                         mass_split,mass_min,c)
            end do

        else

            ! Liquid exchange: evaporation/condensation of surface liquid water.
            ! This branch cannot change the layer structure.
            mass_w(1) = max(mass_w(1) + vflux%vapor_mass,0.0_wp)

        end if

        if (n .gt. 0) then
            t_srf = temperature(1)
        else
            t_srf = c%T0
        end if

        if (n .eq. 0) albedo = c%alpha_ice

        return

    end subroutine apply_snow_surface_vapor_mass_flux

end module snow_surface_fluxes
