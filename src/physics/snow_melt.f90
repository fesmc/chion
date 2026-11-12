module snow_melt
    ! Application of surface melt to a single snowpack column.
    !
    ! Port of Chion.jl/src/processes/melt.jl (_melt_surface_layers!, 38cb742).
    ! The melt ENERGY is computed by the energy solve (WP5); this module
    ! converts the melt potential it stands for into liquid water, consuming
    ! snow layers from the top and restructuring the column as layers are
    ! depleted.
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n, which this routine may DECREASE. See docs/porting_notes.md D8.
    !
    ! CLOSED ENERGY BUDGET: melt_mass is the melt potential, the melt energy
    ! divided by Lm [kg m-2 of ice at T0]. Ice taken from a surface layer below
    ! T0 must first be warmed to T0, so each kilogram of it costs
    ! 1 + ci*(T0 - T(1))/Lm of the potential, and less than melt_mass melts.
    ! Chion.jl charges this only for melt from the surface energy balance, not
    ! for a prescribed melt (a mass forcing); chion's melt always comes from
    ! the energy balance, so it always charges.
    !
    ! The loop consumes layers from the top. Each iteration either removes an
    ! already-empty surface layer (routing its water down, or to runoff if it
    ! was the last layer), or melts what the potential affords of the surface
    ! layer: all of it unless the affordable mass is more than TOL_EMPTY_LAYER
    ! below its mass (the potential is then used up). After each bite the
    ! surface layer is either removed (if depleted below TOL_EMPTY_LAYER) or
    ! merged with the layer below (if it has fallen under mass_min and is not
    ! the only layer). The potential left once the snow is exhausted is
    ! returned: the caller melts the ice beneath with it
    ! (Chion.jl/src/step.jl column_step_core!).
    !
    ! NOTE the two thresholds in the loop are different and deliberately so:
    ! TOL_TINY gates "this layer is already empty, drop it", TOL_EMPTY_LAYER
    ! gates "this layer has just been depleted, drop it". See docs/PLAN.md
    ! section 5, item 1.
    !
    ! On exit the interface temperature is set to T0 (the surface is melting),
    ! and a column that has lost all its snow has its albedo reset to bare ice.
    ! Neither happens on the early return -- deliberate: the energy solve has
    ! already set t_srf = T0 whenever it requested melt.
    !
    ! The two restructuring operations live in snow_layers (WP4):
    ! remove_depleted_surface_and_route_water and merge_layer (k = 1). Julia
    ! passes Ntot = typemax(Int) to the latter; the Fortran routine drops that
    ! argument, because the Julia body never reads it.

    use chion_defs,        only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER, &
                                  chion_const_class
    use snow_layers,       only : remove_depleted_surface_and_route_water, &
                                  merge_layer

    implicit none

    private

    public :: apply_melt

contains

    subroutine apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo_dyn, &
                          mass_split,mass_min,melt_mass,c,melted,remaining)
        ! Chion.jl/src/processes/melt.jl:21-80 (_melt_surface_layers!,
        ! charge_sensible_heat = true).

        implicit none

        real(wp),     intent(INOUT) :: mass(:)         ! (Ntot) [kg m-2] solid
        real(wp),     intent(INOUT) :: mass_w(:)       ! (Ntot) [kg m-2] liquid
        real(wp),     intent(INOUT) :: density(:)      ! (Ntot) [kg m-3]
        real(wp),     intent(INOUT) :: temperature(:)  ! (Ntot) [K]
        integer,      intent(INOUT) :: n               ! active layer count
        real(wp_acc), intent(INOUT) :: runoff          ! [kg m-2] column total
        real(wp),     intent(INOUT) :: t_srf           ! [K] surface diagnostic
        real(wp),     intent(INOUT) :: albedo_dyn      ! [1] dynamic albedo
        real(wp),     intent(IN)    :: mass_split      ! [kg m-2]
        real(wp),     intent(IN)    :: mass_min        ! [kg m-2]
        real(wp),     intent(IN)    :: melt_mass       ! [kg m-2] melt potential (ice at T0)
        type(chion_const_class), intent(IN) :: c
        real(wp_acc), intent(OUT)   :: melted          ! [kg m-2] snow melted
        real(wp_acc), intent(OUT)   :: remaining       ! [kg m-2] potential left (ice at T0)

        ! Local variables
        real(wp_acc) :: dm, cost, affordable
        real(wp)     :: m_srf

        ! _safe_nonnegative(melt_mass)
        remaining = max(real(melt_mass,wp_acc),0.0_wp_acc)

        melted = 0.0_wp_acc

        ! Early return leaves t_srf and albedo_dyn untouched (as in Julia).
        if (remaining .le. 0.0_wp_acc .or. n .le. 0) return

        do while (remaining .gt. 0.0_wp_acc .and. n .gt. 0)

            m_srf = mass(1)

            if (real(m_srf,wp_acc) .le. TOL_TINY) then
                ! Layer is already empty: drop it, routing its residue into
                ! layer 2, or to runoff if it was the last layer.
                call remove_depleted_surface_and_route_water(mass,mass_w,density, &
                                                             temperature,n,runoff,c)
                cycle
            end if

            ! Melt potential consumed per kilogram of ice from this layer.
            cost = 1.0_wp_acc + real(c%ci,wp_acc) &
                   *max(real(c%T0,wp_acc) - real(temperature(1),wp_acc),0.0_wp_acc) &
                   /real(c%Lm,wp_acc)

            affordable = remaining/cost

            if (affordable .lt. real(m_srf,wp_acc) - TOL_EMPTY_LAYER) then
                dm        = affordable
                remaining = 0.0_wp_acc
            else
                dm        = real(m_srf,wp_acc)
                remaining = max(remaining - dm*cost,0.0_wp_acc)
            end if

            mass(1)   = real(real(m_srf,wp_acc)     - dm,wp)
            mass_w(1) = real(real(mass_w(1),wp_acc) + dm,wp)

            melted = melted + dm

            if (real(mass(1),wp_acc) .le. TOL_EMPTY_LAYER) then
                call remove_depleted_surface_and_route_water(mass,mass_w,density, &
                                                             temperature,n,runoff,c)
            else if (n .gt. 1 .and. mass(1) .lt. mass_min) then
                call merge_layer(mass,mass_w,density,temperature,n,1, &
                                 mass_split,mass_min,c)
            end if

        end do

        t_srf = c%T0

        if (n .eq. 0) albedo_dyn = c%alpha_ice

        return

    end subroutine apply_melt

end module snow_melt
