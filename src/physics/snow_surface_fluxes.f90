module snow_surface_fluxes
    ! Surface energy and vapor-mass fluxes evaluated at a KNOWN surface
    ! temperature.
    !
    ! Port of Chion.jl/src/processes/surface_fluxes.jl, together with
    ! _diagnose_latent_heat_flux_coefficients (energy_flux.jl:277-296). The
    ! vapor-pressure helpers of energy_flux.jl:40-95 live in snow_vapor.
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

    use chion_defs, only : wp, wp_acc, TOL_EMPTY_LAYER, io_unit_err, &
                           CHION_ALBEDO_PRESCRIBED, CHION_SEB_SEMIX, &
                           chion_const_class, chion_step_forcing_class
    use snow_column_utils, only : surface_has_snow

    ! Vapour-pressure parameterizations and the BESSI latent flux built from
    ! them. Extracted to snow_vapor so that snow_seb_semix can share them
    ! without a circular dependency on this module.
    use snow_vapor, only : latent_vapor_flux

    ! SEMIX aerodynamic surface scheme, selected by c%seb_scheme. It supplies
    ! the exact-at-known-T turbulent fluxes here, as it supplies the linearized
    ! ones in snow_energy: all three flux sites move together, so the bare-ice
    ! branch and the vapour-mass budget never fall back to BESSI's D_sh while
    ! the energy solve uses r_a (docs/semix_port_scope.md).
    use snow_seb_semix, only : semix_exchange_class, semix_snow_depth, &
                               semix_turbulent_exchange, &
                               semix_sensible_heat_flux, semix_latent_heat_flux, &
                               semix_surface_emissivity, semix_longwave_down, &
                               semix_longwave_flux

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

    type latent_heat_coeff_class
        ! _diagnose_latent_heat_flux_coefficients (energy_flux.jl:282) returns
        ! (linear, constant) IN THAT ORDER -- note this is the OPPOSITE order
        ! to snow_vapor's latent_vapor_flux_lin_class. Both orders are kept as
        ! the Julia source has them, so the two are separate types on purpose.
        real(wp) :: linear             ! [W m-2 K-1]
        real(wp) :: constant           ! [W m-2]
    end type latent_heat_coeff_class

    public :: nonshortwave_flux_class
    public :: bare_ice_flux_class
    public :: bare_ice_ablation_class
    public :: surface_vapor_flux_class
    public :: latent_heat_coeff_class

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
    ! Resolved (exact-at-known-T) surface fluxes
    ! =====================================================================

    pure function resolved_nonshortwave_surface_flux_components(c,forc,surface_temperature, &
                                                               h_snow,has_snow) result(flx)
        ! Chion.jl/src/processes/surface_fluxes.jl:20-45.
        !
        ! Every has_* flag selects a PRESCRIBED value over the internal
        ! parameterization. Note the latent flux has three cases, in order:
        !   has_q_lh                 -> prescribed q_lh
        !   has_relative_humidity    -> BESSI vapor flux
        !   otherwise                -> zero
        !
        ! Under seb_scheme = semix the longwave, sensible and latent terms all
        ! come from the SEMIX scheme instead. That needs two things BESSI does
        ! not: the snow depth, for the roughness blend, and whether the surface
        ! is snow or bare ice, for the emissivity (ebal's mask_snow). Both are
        ! arguments rather than inferred from each other -- a snow column can be
        ! arbitrarily thin without ceasing to be snow. The BESSI scheme ignores
        ! both.
        !
        ! Julia carries a dt_seconds argument here purely to spell zero() in
        ! the right type; it is never used numerically. Dropped (cleanup:
        ! Julia type-dispatch artifact, no behaviour change).

        implicit none

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),                       intent(IN) :: surface_temperature  ! [K]
        real(wp),                       intent(IN) :: h_snow               ! [m]
        logical,                        intent(IN) :: has_snow
        type(nonshortwave_flux_class) :: flx

        ! Local variables
        logical                    :: uses_semix_seb
        type(semix_exchange_class) :: sx

        uses_semix_seb = (c%seb_scheme .eq. CHION_SEB_SEMIX)

        if (uses_semix_seb) then
            sx = semix_turbulent_exchange(c,h_snow,forc%air_temperature, &
                                          surface_temperature,forc%wind_speed, &
                                          forc%air_pressure,forc%relative_humidity, &
                                          forc%has_relative_humidity)
        end if

        ! Longwave. semix absorbs the downwelling flux with the surface
        ! emissivity (snow or ice) where BESSI takes it at face value. The
        ! BESSI branches keep their original grouping so their answer is
        ! unchanged to the last bit.
        if (uses_semix_seb) then
            flx%longwave = semix_longwave_flux(c, &
                               semix_surface_emissivity(c,has_snow), &
                               semix_longwave_down(c,forc%q_lw_down,forc%has_q_lw_down, &
                                                   forc%air_temperature), &
                               surface_temperature)
        else if (forc%has_q_lw_down) then
            flx%longwave = forc%q_lw_down &
                           - c%sigma_sb*c%eps_snow*surface_temperature**4
        else
            flx%longwave = c%sigma_sb*(c%eps_air*forc%air_temperature**4 &
                                       - c%eps_snow*surface_temperature**4)
        end if

        if (forc%has_q_sh) then
            flx%sensible = forc%q_sh
        else if (uses_semix_seb) then
            flx%sensible = semix_sensible_heat_flux(sx,forc%air_temperature, &
                                                    surface_temperature)
        else
            flx%sensible = c%D_sh*(forc%air_temperature - surface_temperature)
        end if

        if (forc%has_q_lh) then
            flx%latent = forc%q_lh
        else if (uses_semix_seb) then
            ! f_lh is already zero without humidity forcing, so this covers the
            ! third case of the BESSI selection too.
            flx%latent = semix_latent_heat_flux(sx)
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

        ! h_snow = 0 and has_snow = .FALSE.: this branch runs only when the
        ! column has no surface snow, so the SEMIX roughness blend collapses to
        ! the bare-ice value and the emissivity is eps_ice.
        nsw = resolved_nonshortwave_surface_flux_components(c,forc,c%T0,0.0_wp,.FALSE.)

        flx%longwave = nsw%longwave
        flx%sensible = nsw%sensible
        flx%latent   = nsw%latent
        flx%rain     = nsw%rain

        return

    end function resolved_bare_ice_surface_flux_components

    pure function resolved_turbulent_latent_heat_flux(c,forc,surface_temperature,h_snow) &
                                                                            result(q_lh)
        ! Chion.jl/src/processes/surface_fluxes.jl:187-200.
        ! The latent-flux-only subset of the three-case selection above, with
        ! the same seb_scheme branch and the same h_snow argument.

        implicit none

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),                       intent(IN) :: surface_temperature   ! [K]
        real(wp),                       intent(IN) :: h_snow                ! [m]
        real(wp) :: q_lh                                                    ! [W m-2]

        ! Local variables
        type(semix_exchange_class) :: sx

        if (forc%has_q_lh) then
            q_lh = forc%q_lh
        else if (c%seb_scheme .eq. CHION_SEB_SEMIX) then
            sx = semix_turbulent_exchange(c,h_snow,forc%air_temperature, &
                                          surface_temperature,forc%wind_speed, &
                                          forc%air_pressure,forc%relative_humidity, &
                                          forc%has_relative_humidity)
            q_lh = semix_latent_heat_flux(sx)
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
        !
        ! seb_scheme = semix splits those two roles apart. The RESERVOIR choice
        ! (solid mass(1) against liquid mass_w(1)) still turns on T0, but the
        ! LATENT HEAT used to convert the flux into mass no longer does: SEMIX
        ! builds f_lh with the latent heat of sublimation at every temperature
        ! (smb_ebal.f90:107), so converting with Lv above the melting point
        ! would overstate the mass by (Lv+Lm)/Lv. Bare ice already uses
        ! (Lv+Lm) unconditionally (bare_ice_ablation_mass), so this also makes
        ! the snow and bare-ice budgets agree under the SEMIX scheme.

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
        real(wp)     :: surface_temperature, q_lh, h_snow, L_exchange
        real(wp_acc) :: vapor

        vflux%vapor_mass       = 0.0_wp
        vflux%sublimation_mass = 0.0_wp
        vflux%latent_heat_flux = 0.0_wp

        if (.not. surface_has_snow(mass,n)) return

        surface_temperature = temperature(1)
        h_snow              = semix_snow_depth(mass,density,n)

        q_lh = resolved_turbulent_latent_heat_flux(c,forc,surface_temperature,h_snow)

        vflux%latent_heat_flux = q_lh

        ! Julia returns early on an exactly-zero flux, so that a column with no
        ! humidity forcing never touches the layer structure. Reproduced.
        if (q_lh .eq. 0.0_wp) return

        if (c%seb_scheme .eq. CHION_SEB_SEMIX) then
            L_exchange = c%Lv + c%Lm
        else if (surface_temperature .lt. c%T0) then
            L_exchange = c%Lv + c%Lm
        else
            L_exchange = c%Lv
        end if

        vapor = real(q_lh,wp_acc)*real(dt_seconds,wp_acc)/real(L_exchange,wp_acc)

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
