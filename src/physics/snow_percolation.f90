module snow_percolation
    ! Liquid-water percolation through a single snowpack column.
    !
    ! Port of Chion.jl/src/processes/percolation.jl (_go_percolation!, the
    ! kernel that step.jl actually calls -- the `state`-level wrapper
    ! go_percolation! only adds an early-return guard and the accumulator
    ! update, both of which belong to the caller here).
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n. See docs/porting_notes.md D8.
    !
    ! SCHEME: a single top-down pass, k = 1..n. Each layer retains liquid water
    ! up to a FIXED volumetric fraction max_lwc of its pore space; anything
    ! above that is pushed into layer k+1, or leaves the column as runoff if
    ! layer k is the lowest active snow layer. Because the excess is deposited
    ! into mass_w(k+1) BEFORE the loop reaches k+1, water cascades all the way
    ! down within this one sweep -- no iteration is needed.
    !
    ! Retention is NOT the density-dependent Coleou-Lesaffre form used by some
    ! other snow models: Chion.jl uses a plain constant volumetric fraction.
    !
    ! Percolation does NOT modify mass, density or temperature. It only moves
    ! liquid water between layers and out of the column.
    !
    ! NOTE on max_lwc: this parameter and chion_const_class%max_lwc_albedo are
    ! DIFFERENT quantities that happen to share a default of 0.1. The first is
    ! the irreducible water saturation of the snow; the second is the LWC at
    ! which the albedo reaches alpha_wet. Keep them independently configurable.
    ! See docs/PLAN.md section 5, item 4.

    use chion_defs,        only : wp, wp_acc, TOL_TINY
    use snow_column_utils, only : is_lowest_active_snow_layer, &
                                  layer_pore_volume, layer_lwc

    implicit none

    private

    real(wp), parameter, public :: PERCOLATION_MAX_LWC_DEFAULT = 0.1_wp

    public :: apply_percolation

contains

    subroutine apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff,max_lwc)
        ! Chion.jl/src/processes/percolation.jl:27-75.
        !
        ! mass_w is mutated in place; runoff is the mass leaving the bottom of
        ! the column during this call (NOT a running total -- the caller adds
        ! it to its own wp_acc accumulator).

        implicit none

        real(wp),     intent(IN)    :: mass(:)      ! (Ntot) [kg m-2] solid mass
        real(wp),     intent(INOUT) :: mass_w(:)    ! (Ntot) [kg m-2] liquid water
        real(wp),     intent(IN)    :: density(:)   ! (Ntot) [kg m-3]
        integer,      intent(IN)    :: n            ! number of active layers
        real(wp),     intent(IN)    :: rho_i        ! [kg m-3] ice density
        real(wp),     intent(IN)    :: rho_w        ! [kg m-3] water density
        real(wp_acc), intent(OUT)   :: runoff       ! [kg m-2] this call only
        real(wp), optional, intent(IN) :: max_lwc   ! [1] irreducible saturation

        ! Local variables
        integer      :: k
        real(wp)     :: m_s, m_w
        real(wp_acc) :: lwc_max, phi, lwc, excess

        lwc_max = real(PERCOLATION_MAX_LWC_DEFAULT,wp_acc)
        if (present(max_lwc)) lwc_max = real(max_lwc,wp_acc)

        runoff = 0.0_wp_acc

        do k = 1, n

            m_s = mass(k)
            m_w = mass_w(k)

            ! === Case A: layer has no solid mass ==========================
            ! The layer cannot retain water, so its water leaves the column
            ! IMMEDIATELY as runoff -- it is NOT passed to layer k+1.
            !
            ! Chion.jl/docs/src/processes/percolation.md says the opposite
            ! ("routed onward"). The CODE (percolation.jl:44-48) is
            ! authoritative and is what is reproduced here. Reported upstream;
            ! see docs/porting_notes.md.
            if (m_s .le. 0.0_wp) then
                runoff    = runoff + real(m_w,wp_acc)
                mass_w(k) = 0.0_wp
                cycle
            end if

            ! === Case B: pore space has collapsed =========================
            ! phi must be evaluated in wp_acc: it is a difference of nearly
            ! equal numbers, and in sp it is quantized far above TOL_TINY so
            ! this guard could never fire. See docs/porting_notes.md D1.
            phi = layer_pore_volume(m_s,density(k),rho_i)

            if (phi .le. TOL_TINY) then
                ! All the water is pushed onward: no pore space to hold it.
                mass_w(k) = 0.0_wp
                call route_excess(mass,mass_w,n,k,real(m_w,wp_acc),runoff)
                cycle
            end if

            ! === Case C: normal retention =================================
            lwc = layer_lwc(m_s,m_w,density(k),rho_i,rho_w)

            if (lwc .gt. lwc_max) then
                excess    = (lwc - lwc_max)*phi*real(rho_w,wp_acc)
                mass_w(k) = real(real(m_w,wp_acc) - excess,wp)
                call route_excess(mass,mass_w,n,k,excess,runoff)
            end if

        end do

        return

    end subroutine apply_percolation

    subroutine route_excess(mass,mass_w,n,k,excess,runoff)
        ! Send `excess` either to layer k+1 or out of the column.
        !
        ! The "is this the lowest layer" test is
        ! snow_column_utils::is_lowest_active_snow_layer, which is written as a
        ! nested if precisely because Fortran does not short-circuit .or. and
        ! k == n == Ntot would otherwise read mass(Ntot+1). Do not inline it
        ! back into a single expression. See docs/porting_notes.md D10.

        implicit none

        real(wp),     intent(IN)    :: mass(:)
        real(wp),     intent(INOUT) :: mass_w(:)
        integer,      intent(IN)    :: n
        integer,      intent(IN)    :: k
        real(wp_acc), intent(IN)    :: excess
        real(wp_acc), intent(INOUT) :: runoff

        if (is_lowest_active_snow_layer(mass,n,k)) then
            runoff = runoff + excess
        else
            mass_w(k+1) = real(real(mass_w(k+1),wp_acc) + excess,wp)
        end if

        return

    end subroutine route_excess

end module snow_percolation
