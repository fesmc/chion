module snow_column_utils
    ! Shared predicates and reductions over a single snowpack column.
    !
    ! Port of Chion.jl/src/column_state_utils.jl, plus _surface_has_snow
    ! (which lives in Chion.jl/src/processes/energy_flux.jl:303) and
    ! _is_lowest_active_snow_layer (processes/percolation.jl:11).
    !
    ! Julia's accessor helpers (_get_layer, _set_layer!, _get_scalar, ...)
    ! have no Fortran analogue: they exist only to let the same kernel run over
    ! a Matrix on GPU and a Vector on CPU. In Fortran they collapse to plain
    ! array indexing. What remains, and what this module provides, are the
    ! PREDICATES and their exact thresholds.
    !
    ! CALLING CONVENTION: routines take a contiguous column slice, e.g.
    ! mass(:,icol), together with the active layer count n. Layer arrays are
    ! (Ntot,ncol) and Fortran is column-major, so a slice is contiguous and
    ! costs nothing. See docs/porting_notes.md D8.
    !
    ! THRESHOLDS: Chion.jl uses three different "is this layer empty" tests,
    ! deliberately and not interchangeably:
    !     > 0                 percolation entry, is_lowest_active_snow_layer
    !     > TOL_TINY          column_has_liquid_water, LWC, bulk_snow_density
    !     > TOL_EMPTY_LAYER   surface_has_snow, constant albedo, albedo update
    ! Do not unify them. See docs/PLAN.md section 5, item 1.

    use chion_defs, only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER

    implicit none

    private

    public :: surface_has_snow
    public :: column_has_liquid_water
    public :: is_lowest_active_snow_layer
    public :: bulk_snow_density
    public :: total_snow_water_mass
    public :: layer_pore_volume
    public :: layer_lwc

contains

    pure function surface_has_snow(mass,n) result(has_snow)
        ! Chion.jl/src/processes/energy_flux.jl:303-305.
        ! Gates the bare-ice branch of the column step, so this is one of the
        ! most consequential predicates in the model.
        !
        ! NOTE the threshold is TOL_EMPTY_LAYER, while the energy solve's own
        ! early exit uses mass(1) <= 0 (energy_flux.jl:343). The two differ.

        implicit none

        real(wp), intent(IN) :: mass(:)     ! (Ntot) [kg m-2] solid mass, layer 1 = surface
        integer,  intent(IN) :: n           ! number of active layers
        logical :: has_snow

        has_snow = .FALSE.
        if (n .gt. 0) then
            if (mass(1) .gt. TOL_EMPTY_LAYER) has_snow = .TRUE.
        end if

        return

    end function surface_has_snow

    pure function column_has_liquid_water(mass_w,n) result(has_water)
        ! Chion.jl/src/column_state_utils.jl:111-123.
        !
        ! NOTE the Julia docstring says "above the empty-layer tolerance" but
        ! the code uses EPS_TINY, not EPS_EMPTY_LAYER. The code is authoritative;
        ! reported upstream. See docs/porting_notes.md.

        implicit none

        real(wp), intent(IN) :: mass_w(:)   ! (Ntot) [kg m-2] liquid water mass
        integer,  intent(IN) :: n
        logical :: has_water

        ! Local variables
        integer :: k

        has_water = .FALSE.

        do k = 1, n
            if (mass_w(k) .gt. TOL_TINY) then
                has_water = .TRUE.
                exit
            end if
        end do

        return

    end function column_has_liquid_water

    pure function is_lowest_active_snow_layer(mass,n,k) result(is_lowest)
        ! Chion.jl/src/processes/percolation.jl:11-19:
        !     k == n_active || mass(k+1) <= 0
        ! Determines whether percolating water leaves the column as runoff or
        ! moves into the layer below.
        !
        ! IMPORTANT: Julia's `||` short-circuits, so mass(k+1) is never accessed
        ! when k == n. Fortran does NOT guarantee short-circuit evaluation of
        ! .or., and k == n == Ntot would read out of bounds. Hence the nested
        ! form below -- do not "simplify" it back to a single .or. expression.

        implicit none

        real(wp), intent(IN) :: mass(:)     ! (Ntot) [kg m-2]
        integer,  intent(IN) :: n           ! number of active layers
        integer,  intent(IN) :: k           ! layer index to test
        logical :: is_lowest

        if (k .ge. n) then
            is_lowest = .TRUE.
        else if (mass(k+1) .le. 0.0_wp) then
            is_lowest = .TRUE.
        else
            is_lowest = .FALSE.
        end if

        return

    end function is_lowest_active_snow_layer

    pure function bulk_snow_density(mass,density,n) result(rho_bulk)
        ! Chion.jl/src/column_state_utils.jl:52-78.
        ! Mass-over-thickness, summing BOTH numerator and denominator only over
        ! layers that pass the guards.
        !
        ! NOTE this differs from the diagnostics version (_column_summary,
        ! Chion.jl/src/diagnostics.jl:123), whose numerator is the
        ! UNCONDITIONAL mass sum. The two disagree whenever a layer has mass but
        ! near-zero density. Both are ported, separately and deliberately.
        !
        ! Accumulated in wp_acc: this is a sum over layers feeding a division.

        implicit none

        real(wp), intent(IN) :: mass(:)     ! (Ntot) [kg m-2]
        real(wp), intent(IN) :: density(:)  ! (Ntot) [kg m-3]
        integer,  intent(IN) :: n
        real(wp) :: rho_bulk

        ! Local variables
        integer      :: k
        real(wp_acc) :: total_mass, total_thickness

        total_mass      = 0.0_wp_acc
        total_thickness = 0.0_wp_acc

        do k = 1, n
            if (mass(k) .gt. 0.0_wp .and. real(density(k),wp_acc) .gt. TOL_TINY) then
                total_mass      = total_mass      + real(mass(k),wp_acc)
                total_thickness = total_thickness + real(mass(k),wp_acc)/real(density(k),wp_acc)
            end if
        end do

        if (total_mass .le. 0.0_wp_acc .or. total_thickness .le. TOL_TINY) then
            rho_bulk = 0.0_wp
        else
            rho_bulk = real(total_mass/total_thickness,wp)
        end if

        return

    end function bulk_snow_density

    pure function total_snow_water_mass(mass,mass_w,n) result(total)
        ! Chion.jl/src/column_state_utils.jl:86-103.
        !
        ! NOTE this CLIPS negative layer values at zero, whereas the diagnostics
        ! reduction (_column_summary) does not. The two therefore disagree if
        ! any layer mass goes negative. Both behaviours are ported as-is.
        !
        ! Returned in wp_acc: this is a conservation quantity, summed and
        ! compared against cumulative accumulators.

        implicit none

        real(wp), intent(IN) :: mass(:)     ! (Ntot) [kg m-2] solid
        real(wp), intent(IN) :: mass_w(:)   ! (Ntot) [kg m-2] liquid
        integer,  intent(IN) :: n
        real(wp_acc) :: total

        ! Local variables
        integer :: k

        total = 0.0_wp_acc

        do k = 1, n
            total = total + max(real(mass(k),wp_acc),  0.0_wp_acc) &
                          + max(real(mass_w(k),wp_acc),0.0_wp_acc)
        end do

        return

    end function total_snow_water_mass

    pure function layer_pore_volume(m_s,rho,rho_i) result(phi)
        ! Pore volume of one layer [m]:  phi = m_s/rho - m_s/rho_i
        ! Chion.jl/src/processes/percolation.jl:50 and processes/albedo.jl:47.
        !
        ! MUST be computed in wp_acc. This is a difference of two nearly equal
        ! numbers as rho -> rho_i, and it feeds a DIVISION (see layer_lwc). In
        ! sp the result is quantized at ~1e-5 m, which both destroys the LWC
        ! and makes the "phi <= TOL_TINY" guard unable to fire at all.
        ! Measured in validation/precision_eval.f90; see docs/PLAN.md section 3.1.
        !
        ! Callers must still apply the phi <= TOL_TINY guard themselves: a
        ! collapsed pore space is a physical branch (water is pushed onward),
        ! not an error, and the two call sites handle it differently.

        implicit none

        real(wp), intent(IN) :: m_s          ! [kg m-2] solid mass of the layer
        real(wp), intent(IN) :: rho          ! [kg m-3] layer density
        real(wp), intent(IN) :: rho_i        ! [kg m-3] ice density
        real(wp_acc) :: phi

        phi = real(m_s,wp_acc)/real(rho,wp_acc) - real(m_s,wp_acc)/real(rho_i,wp_acc)

        return

    end function layer_pore_volume

    pure function layer_lwc(m_s,m_w,rho,rho_i,rho_w) result(lwc)
        ! Volumetric liquid water content of one layer [1]:
        !     lwc = m_w / rho_w / phi
        ! Chion.jl/src/processes/percolation.jl:63, processes/albedo.jl:52.
        !
        ! Returns 0 when the pore space is collapsed or the layer has no solid
        ! mass, matching both Julia call sites' guards. Not clamped above 1:
        ! the albedo wetness law relies on lwc being able to exceed
        ! max_lwc_albedo, and clamps downstream.

        implicit none

        real(wp), intent(IN) :: m_s          ! [kg m-2] solid mass
        real(wp), intent(IN) :: m_w          ! [kg m-2] liquid water mass
        real(wp), intent(IN) :: rho          ! [kg m-3] layer density
        real(wp), intent(IN) :: rho_i        ! [kg m-3] ice density
        real(wp), intent(IN) :: rho_w        ! [kg m-3] water density
        real(wp_acc) :: lwc

        ! Local variables
        real(wp_acc) :: phi

        lwc = 0.0_wp_acc

        if (m_s .le. 0.0_wp) return
        if (rho .le. 0.0_wp) return

        phi = layer_pore_volume(m_s,rho,rho_i)

        if (phi .le. TOL_TINY) return

        lwc = max(real(m_w,wp_acc),0.0_wp_acc)/real(rho_w,wp_acc)/phi

        return

    end function layer_lwc

end module snow_column_utils
