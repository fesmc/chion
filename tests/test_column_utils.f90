program test_column_utils
    ! WP3 acceptance test: each predicate flips at exactly its documented
    ! threshold, and the three thresholds are demonstrably different.

    use chion_defs,        only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER
    use snow_column_utils

    implicit none

    integer, parameter :: Ntot = 5

    real(wp) :: mass(Ntot), mass_w(Ntot), density(Ntot)
    integer  :: nfail

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP3 acceptance test: snow_column_utils"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! === surface_has_snow: threshold is TOL_EMPTY_LAYER ==================
    write(*,"(a)") "--- surface_has_snow (threshold TOL_EMPTY_LAYER = 1e-10) ---"

    mass = 0.0_wp

    mass(1) = 100.0_wp
    call check("n=0 always false, even with mass present", &
               .not. surface_has_snow(mass,0), nfail)

    call check("n=1, mass=100 -> true", surface_has_snow(mass,1), nfail)

    mass(1) = 0.0_wp
    call check("mass=0 -> false", .not. surface_has_snow(mass,1), nfail)

    ! Straddle the threshold. In sp, TOL_EMPTY_LAYER = 1e-10 is representable
    ! as a magnitude, so this is a meaningful test of the comparison itself.
    mass(1) = real(TOL_EMPTY_LAYER,wp)*0.5_wp
    call check("mass just below TOL_EMPTY_LAYER -> false", &
               .not. surface_has_snow(mass,1), nfail)

    mass(1) = real(TOL_EMPTY_LAYER,wp)*2.0_wp
    call check("mass just above TOL_EMPTY_LAYER -> true", &
               surface_has_snow(mass,1), nfail)

    ! The distinguishing case: TOL_TINY < mass < TOL_EMPTY_LAYER.
    ! This is snow by the liquid-water threshold but NOT by the surface
    ! threshold. If someone unifies the two constants, this check fails.
    mass(1) = 1.0e-11_wp
    call check("mass between TOL_TINY and TOL_EMPTY_LAYER -> false", &
               .not. surface_has_snow(mass,1), nfail)

    ! === column_has_liquid_water: threshold is TOL_TINY ==================
    write(*,*)
    write(*,"(a)") "--- column_has_liquid_water (threshold TOL_TINY = 1e-12) ---"

    mass_w = 0.0_wp
    call check("no water -> false", .not. column_has_liquid_water(mass_w,3), nfail)

    mass_w(3) = 5.0_wp
    call check("water in layer 3, n=3 -> true", column_has_liquid_water(mass_w,3), nfail)
    call check("water in layer 3, n=2 -> false (inactive layer ignored)", &
               .not. column_has_liquid_water(mass_w,2), nfail)

    mass_w    = 0.0_wp
    mass_w(1) = 1.0e-11_wp
    call check("water above TOL_TINY, below TOL_EMPTY_LAYER -> true", &
               column_has_liquid_water(mass_w,1), nfail)

    mass_w(1) = 1.0e-13_wp
    call check("water below TOL_TINY -> false", &
               .not. column_has_liquid_water(mass_w,1), nfail)

    ! === is_lowest_active_snow_layer: threshold is exactly 0 =============
    write(*,*)
    write(*,"(a)") "--- is_lowest_active_snow_layer (threshold exactly 0) ---"

    mass    = 0.0_wp
    mass(1) = 100.0_wp
    mass(2) = 100.0_wp
    mass(3) = 100.0_wp

    call check("k<n with mass below -> not lowest", &
               .not. is_lowest_active_snow_layer(mass,3,1), nfail)
    call check("k=n -> lowest", is_lowest_active_snow_layer(mass,3,3), nfail)

    mass(3) = 0.0_wp
    call check("k<n but layer below is massless -> lowest", &
               is_lowest_active_snow_layer(mass,3,2), nfail)

    mass(3) = 1.0e-11_wp
    call check("layer below has tiny but positive mass -> not lowest", &
               .not. is_lowest_active_snow_layer(mass,3,2), nfail)

    ! Bounds safety: k = n = Ntot must not read mass(Ntot+1).
    mass = 100.0_wp
    call check("k=n=Ntot does not read out of bounds", &
               is_lowest_active_snow_layer(mass,Ntot,Ntot), nfail)

    ! === bulk_snow_density ==============================================
    write(*,*)
    write(*,"(a)") "--- bulk_snow_density ---"

    mass    = 0.0_wp
    density = 0.0_wp

    call check_val("n=0 -> 0", bulk_snow_density(mass,density,0), 0.0_wp, nfail)

    ! Two layers: 100 kg m-2 at 400, 300 kg m-2 at 600.
    ! thickness = 100/400 + 300/600 = 0.25 + 0.5 = 0.75 m
    ! bulk      = 400/0.75 = 533.3333
    mass(1)    = 100.0_wp ; density(1) = 400.0_wp
    mass(2)    = 300.0_wp ; density(2) = 600.0_wp
    call check_val("two layers", bulk_snow_density(mass,density,2), &
                   400.0_wp/0.75_wp, nfail)

    ! A layer with mass but zero density is excluded from BOTH sums here.
    ! bulk = 100/(100/400) = 400
    mass(2) = 300.0_wp ; density(2) = 0.0_wp
    call check_val("zero-density layer excluded from both sums", &
                   bulk_snow_density(mass,density,2), 400.0_wp, nfail)

    ! === total_snow_water_mass ==========================================
    write(*,*)
    write(*,"(a)") "--- total_snow_water_mass (clips negatives) ---"

    mass   = 0.0_wp ; mass_w = 0.0_wp
    mass(1) = 100.0_wp ; mass_w(1) = 10.0_wp
    mass(2) = 200.0_wp ; mass_w(2) = 20.0_wp
    call check_acc("sum of solid + liquid", &
                   total_snow_water_mass(mass,mass_w,2), 330.0_wp_acc, nfail)

    mass(2) = -50.0_wp
    call check_acc("negative layer mass is clipped to zero", &
                   total_snow_water_mass(mass,mass_w,2), 130.0_wp_acc, nfail)

    ! === layer_pore_volume: the dp requirement ==========================
    write(*,*)
    write(*,"(a)") "--- layer_pore_volume (must be dp; see PLAN section 3.1) ---"

    call check("returns wp_acc", kind(layer_pore_volume(300.0_wp,400.0_wp,917.0_wp)) &
               .eq. wp_acc, nfail)

    ! m=300, rho=400, rho_i=917: phi = 300/400 - 300/917 = 0.75 - 0.327154 = 0.422846
    call check_acc("phi for typical firn layer", &
                   layer_pore_volume(300.0_wp,400.0_wp,917.0_wp), &
                   300.0_wp_acc/400.0_wp_acc - 300.0_wp_acc/917.0_wp_acc, nfail)

    call check_acc("phi = 0 at ice density", &
                   layer_pore_volume(300.0_wp,917.0_wp,917.0_wp), 0.0_wp_acc, nfail)

    ! The reason this must be dp: computed in sp the same expression is
    ! quantized far above TOL_TINY, so the guard could never fire.
    call check("phi near ice density is resolved, not quantized to zero", &
               layer_pore_volume(300.0_wp,916.99_wp,917.0_wp) .gt. TOL_TINY, nfail)

    ! === layer_lwc ======================================================
    write(*,*)
    write(*,"(a)") "--- layer_lwc ---"

    call check_acc("no solid mass -> 0", &
                   layer_lwc(0.0_wp,10.0_wp,400.0_wp,917.0_wp,1000.0_wp), 0.0_wp_acc, nfail)

    call check_acc("collapsed pore space -> 0", &
                   layer_lwc(300.0_wp,10.0_wp,917.0_wp,917.0_wp,1000.0_wp), 0.0_wp_acc, nfail)

    ! phi = 0.422846 m; lwc = 10/1000/0.422846 = 0.023649
    call check_acc("typical wet layer", &
                   layer_lwc(300.0_wp,10.0_wp,400.0_wp,917.0_wp,1000.0_wp), &
                   10.0_wp_acc/1000.0_wp_acc &
                   /(300.0_wp_acc/400.0_wp_acc - 300.0_wp_acc/917.0_wp_acc), nfail)

    call check("lwc is NOT clamped at max_lwc_albedo (albedo law needs this)", &
               layer_lwc(300.0_wp,300.0_wp,400.0_wp,917.0_wp,1000.0_wp) .gt. 0.1_wp_acc, nfail)

    ! === Summary ========================================================
    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP3: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP3: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

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

    subroutine check_val(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value
        real(wp),         intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp) :: tol

        tol = max(abs(expected)*8.0_wp*epsilon(1.0_wp),tiny(1.0_wp))

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_val

    subroutine check_acc(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp_acc) :: tol

        ! Inputs are sp, so the achievable accuracy is sp-limited even though
        ! the arithmetic is dp.
        tol = max(abs(expected)*1.0e-6_wp_acc,1.0e-30_wp_acc)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g16.8)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_acc

end program test_column_utils
