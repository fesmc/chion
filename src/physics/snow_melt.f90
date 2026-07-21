module snow_melt
    ! Application of surface melt to a single snowpack column.
    !
    ! Port of Chion.jl/src/processes/melt.jl (_apply_melt!). The melt ENERGY is
    ! computed by the energy solve (WP5); this module only converts a given
    ! melt mass into liquid water, consuming snow layers from the top and
    ! restructuring the column as layers are depleted.
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n, which this routine may DECREASE. See docs/porting_notes.md D8.
    !
    ! TWO PATHS:
    !   Fast path -- when the requested melt is comfortably smaller than the
    !     surface layer (melt < mass(1) - TOL_EMPTY_LAYER), the whole amount is
    !     moved from mass(1) into mass_w(1) and nothing else happens. Note this
    !     path returns the FULL requested melt.
    !   General path -- otherwise a loop consumes layers from the top. Each
    !     iteration either removes an already-empty surface layer (routing its
    !     water down, or to runoff if it was the last layer), or melts up to
    !     the whole surface layer. After each bite the surface layer is either
    !     removed (if depleted below TOL_EMPTY_LAYER) or merged with the layer
    !     below (if it has fallen under mass_min and is not the only layer).
    !     This path returns only the mass actually melted, which is less than
    !     requested once the column runs out of snow -- the caller uses that
    !     shortfall to drive bare-ice melt (Chion.jl/src/step.jl:346-350).
    !
    ! NOTE the two thresholds in the loop are different and deliberately so:
    ! TOL_TINY gates "this layer is already empty, drop it", TOL_EMPTY_LAYER
    ! gates "this layer has just been depleted, drop it". See docs/PLAN.md
    ! section 5, item 1.
    !
    ! On exit the surface temperature diagnostic is refreshed, and a column
    ! that has lost all its snow has its albedo reset to bare ice. Neither
    ! happens on the early return or on the fast path -- also deliberate.
    !
    ! The two restructuring operations live in snow_layers (WP4):
    ! remove_depleted_surface_and_route_water and merge_surface_layer. Julia
    ! passes Ntot = typemax(Int) to the latter; the Fortran routine drops that
    ! argument, because the Julia body never reads it.

    use chion_defs,        only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER, &
                                  chion_const_class
    use snow_layers,       only : remove_depleted_surface_and_route_water, &
                                  merge_surface_layer

    implicit none

    private

    public :: apply_melt

contains

    subroutine apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo_dyn, &
                          mass_split,mass_min,melt_mass,c,melted)
        ! Chion.jl/src/processes/melt.jl:12-65.

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
        real(wp),     intent(IN)    :: melt_mass       ! [kg m-2] requested melt
        type(chion_const_class), intent(IN) :: c
        real(wp_acc), intent(OUT)   :: melted          ! [kg m-2] actually melted

        ! Local variables
        real(wp_acc) :: remaining, dm
        real(wp)     :: m_srf

        ! _safe_nonnegative(melt_mass)
        remaining = max(real(melt_mass,wp_acc),0.0_wp_acc)

        melted = 0.0_wp_acc

        ! Early return leaves t_srf and albedo_dyn untouched (as in Julia).
        if (remaining .le. 0.0_wp_acc .or. n .le. 0) return

        ! === Fast path ====================================================
        m_srf = mass(1)

        if (remaining .lt. real(m_srf,wp_acc) - TOL_EMPTY_LAYER) then
            mass(1)   = real(real(m_srf,wp_acc)   - remaining,wp)
            mass_w(1) = real(real(mass_w(1),wp_acc) + remaining,wp)
            melted    = remaining
            return
        end if

        ! === General path =================================================
        do while (remaining .gt. 0.0_wp_acc .and. n .gt. 0)

            m_srf = mass(1)

            if (real(m_srf,wp_acc) .le. TOL_TINY) then
                ! Layer is already empty: drop it, routing its residual water
                ! into layer 2, or to runoff if it was the last layer.
                call remove_depleted_surface_and_route_water(mass,mass_w,density, &
                                                             temperature,n,runoff,c)
                cycle
            end if

            dm = min(real(m_srf,wp_acc),remaining)

            mass(1)   = real(real(m_srf,wp_acc)     - dm,wp)
            mass_w(1) = real(real(mass_w(1),wp_acc) + dm,wp)

            remaining = remaining - dm
            melted    = melted    + dm

            if (real(mass(1),wp_acc) .le. TOL_EMPTY_LAYER) then
                call remove_depleted_surface_and_route_water(mass,mass_w,density, &
                                                             temperature,n,runoff,c)
            else if (n .gt. 1 .and. mass(1) .lt. mass_min) then
                call merge_surface_layer(mass,mass_w,density,temperature,n, &
                                         mass_split,mass_min,c)
            end if

        end do

        if (n .gt. 0) then
            t_srf = temperature(1)
        else
            t_srf = c%T0
            albedo_dyn = c%alpha_ice
        end if

        return

    end subroutine apply_melt

end module snow_melt
