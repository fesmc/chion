program test_diagnostics
    ! WP10 acceptance test: snow_diagnostics reproduces Chion.jl's
    ! _column_summary on a hand-computed 3-layer example, INCLUDING its
    ! mutually inconsistent inclusion criteria.
    !
    ! The central assertions are the ones that a well-meaning cleanup would
    ! break:
    !   * a zero-density layer with mass is excluded from `thickness` but
    !     INCLUDED in bulk_density's numerator, so `bulk_density` here does
    !     NOT equal snow_column_utils::bulk_snow_density on the same column;
    !   * a negative layer mass is NOT clipped by wet_mass / liquid_water,
    !     so they do NOT equal snow_column_utils::total_snow_water_mass;
    !   * stale values in layers above n are ignored, not zeroed.
    ! All three asymmetries are intentional. See docs/PLAN.md section 5.

    use chion_defs,        only : wp, wp_acc, TOL_TINY, chion_const_class, chion_const_init
    use snow_column_utils, only : bulk_snow_density, total_snow_water_mass
    use snow_diagnostics

    implicit none

    integer, parameter :: Ntot = 5
    integer, parameter :: ncol = 3

    type(chion_const_class)          :: c
    type(snow_column_snapshot_class) :: snap

    real(wp) :: mass(Ntot), mass_w(Ntot), density(Ntot), temperature(Ntot)
    real(wp) :: thk, wet, rho_b, lwm

    real(wp)     :: mass2(Ntot,ncol), mass_w2(Ntot,ncol), density2(Ntot,ncol)
    integer      :: n_lay(ncol)
    real(wp_acc) :: mass_base(ncol)
    real(wp)     :: thk_v(ncol), wet_v(ncol), rho_v(ncol), lwm_v(ncol), base_v(ncol)

    real(wp_acc) :: thk_exp, rho_exp
    integer      :: nfail

    call chion_const_init(c)

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP10 acceptance test: snow_diagnostics"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! =====================================================================
    ! Hand-computed 3-layer column with every edge case in it.
    !
    !   k  mass    mass_w   density    comment
    !   1   100.0    10.0     400.0    normal layer
    !   2   300.0     5.0       0.0    ZERO DENSITY, but has mass
    !   3   -50.0     2.0     500.0    NEGATIVE mass
    !   4   999.0   999.0     999.0    STALE, beyond n = 3
    !   5   777.0   777.0     777.0    STALE, beyond n = 3
    !
    ! Expected (Chion.jl diagnostics.jl:123-143):
    !   thickness    = 100/400                       = 0.25 m
    !                  (k=2 excluded: density <= TOL_TINY)
    !                  (k=3 excluded: mass <= 0)
    !   wet_mass     = (100+10) + (300+5) + (-50+2)  = 367.0 kg m-2
    !                  unconditional, negatives NOT clipped
    !   liquid_water = 10 + 5 + 2                    = 17.0 kg m-2
    !   solid total  = 100 + 300 - 50                = 350.0 kg m-2
    !                  UNCONDITIONAL -- includes k=2 and k=3
    !   bulk_density = 350 / 0.25                    = 1400.0 kg m-3
    !
    ! 1400 kg m-3 is denser than ice. That is the point: the numerator
    ! counts mass whose volume was never counted.
    ! =====================================================================

    write(*,"(a)") "--- column_summary: hand-computed 3-layer column ---"

    mass        = 0.0_wp
    mass_w      = 0.0_wp
    density     = 0.0_wp
    temperature = 0.0_wp

    mass(1) =  100.0_wp ; mass_w(1) = 10.0_wp ; density(1) = 400.0_wp ; temperature(1) = 270.0_wp
    mass(2) =  300.0_wp ; mass_w(2) =  5.0_wp ; density(2) =   0.0_wp ; temperature(2) = 271.0_wp
    mass(3) =  -50.0_wp ; mass_w(3) =  2.0_wp ; density(3) = 500.0_wp ; temperature(3) = 272.0_wp
    mass(4) =  999.0_wp ; mass_w(4) = 999.0_wp; density(4) = 999.0_wp ; temperature(4) = 999.0_wp
    mass(5) =  777.0_wp ; mass_w(5) = 777.0_wp; density(5) = 777.0_wp ; temperature(5) = 777.0_wp

    call column_summary(mass,mass_w,density,3,thk,wet,rho_b,lwm)

    thk_exp = 100.0_wp_acc/400.0_wp_acc
    rho_exp = 350.0_wp_acc/thk_exp

    call check_val("thickness excludes zero-density and negative-mass layers", &
                   thk, real(thk_exp,wp), nfail)
    call check_val("wet_mass is unconditional and unclipped", &
                   wet, 367.0_wp, nfail)
    call check_val("liquid_water is unconditional and unclipped", &
                   lwm, 17.0_wp, nfail)
    call check_val("bulk_density uses the UNCONDITIONAL mass sum", &
                   rho_b, real(rho_exp,wp), nfail)

    call check("bulk_density exceeds ice density here, by construction", &
               rho_b .gt. 917.0_wp, nfail)

    ! =====================================================================
    ! THE INTENTIONAL DISAGREEMENTS. A "cleanup" that unifies these
    ! reductions must fail here.
    ! =====================================================================

    write(*,*)
    write(*,"(a)") "--- intentional disagreement with snow_column_utils ---"

    ! bulk_snow_density restricts BOTH sums to qualifying layers:
    !   mass = 100 only (k=2 zero density, k=3 negative mass)
    !   thickness = 0.25
    !   -> 400 kg m-3, versus column_summary's 1400.
    call check_val("bulk_snow_density on the same column (both sums restricted)", &
                   bulk_snow_density(mass,density,3), 400.0_wp, nfail)

    call check("column_summary bulk_density /= bulk_snow_density (INTENTIONAL)", &
               abs(rho_b - bulk_snow_density(mass,density,3)) .gt. 1.0_wp, nfail)

    ! total_snow_water_mass clips negatives: (100+10) + (300+5) + (0+2) = 417
    call check_acc("total_snow_water_mass clips the negative layer", &
                   total_snow_water_mass(mass,mass_w,3), 417.0_wp_acc, nfail)

    call check("column_summary wet_mass /= total_snow_water_mass (INTENTIONAL)", &
               abs(real(wet,wp_acc) - total_snow_water_mass(mass,mass_w,3)) &
               .gt. 1.0_wp_acc, nfail)

    ! =====================================================================
    ! Stale layers beyond n
    ! =====================================================================

    write(*,*)
    write(*,"(a)") "--- layers beyond n are ignored, not zeroed ---"

    ! Layers 4 and 5 carry huge stale values. If any of them leaked into the
    ! sums, wet_mass would be ~3900 rather than 367.
    call check("stale layers 4,5 did not contaminate wet_mass", &
               abs(wet - 367.0_wp) .lt. 1.0_wp, nfail)

    ! And they are still there afterwards -- column_summary is pure.
    call check_val("layer 4 mass untouched by column_summary", mass(4), 999.0_wp, nfail)
    call check_val("layer 5 mass untouched by column_summary", mass(5), 777.0_wp, nfail)

    ! Raising n to 5 must change the answer, proving 1:n is the real bound.
    call column_summary(mass,mass_w,density,5,thk,wet,rho_b,lwm)
    call check_val("n=5 includes the stale layers (367 + 1998 + 1554)", &
                   wet, 367.0_wp + 1998.0_wp + 1554.0_wp, nfail)

    ! =====================================================================
    ! Degenerate columns
    ! =====================================================================

    write(*,*)
    write(*,"(a)") "--- degenerate columns ---"

    call column_summary(mass,mass_w,density,0,thk,wet,rho_b,lwm)
    call check_val("n=0 thickness", thk, 0.0_wp, nfail)
    call check_val("n=0 wet_mass",  wet, 0.0_wp, nfail)
    call check_val("n=0 liquid",    lwm, 0.0_wp, nfail)
    call check_val("n=0 bulk_density guarded to zero", rho_b, 0.0_wp, nfail)

    ! A single zero-density layer: thickness stays 0, so the bulk_density
    ! guard fires and no division by zero occurs (this is what makes the
    ! -ffpe-trap=zero build pass).
    mass    = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1) = 250.0_wp
    call column_summary(mass,mass_w,density,1,thk,wet,rho_b,lwm)
    call check_val("all-zero-density column: thickness 0", thk, 0.0_wp, nfail)
    call check_val("all-zero-density column: bulk_density guarded to 0", rho_b, 0.0_wp, nfail)
    call check_val("all-zero-density column: wet_mass still counts the mass", wet, 250.0_wp, nfail)

    ! Density exactly at the TOL_TINY threshold must be EXCLUDED (strict >).
    density(1) = real(TOL_TINY,wp)
    call column_summary(mass,mass_w,density,1,thk,wet,rho_b,lwm)
    call check_val("density exactly at TOL_TINY is excluded (strict >)", thk, 0.0_wp, nfail)

    ! Mass exactly zero must be EXCLUDED (strict >), even with good density.
    mass(1) = 0.0_wp ; density(1) = 400.0_wp
    call column_summary(mass,mass_w,density,1,thk,wet,rho_b,lwm)
    call check_val("mass exactly zero is excluded (strict >)", thk, 0.0_wp, nfail)

    ! =====================================================================
    ! summarize_domain_state / summarize_year_state
    ! =====================================================================

    write(*,*)
    write(*,"(a)") "--- summarize_domain_state / summarize_year_state ---"

    mass2 = 0.0_wp ; mass_w2 = 0.0_wp ; density2 = 0.0_wp

    ! Column 1: the hand-computed edge-case column.
    mass2(1,1) =  100.0_wp ; mass_w2(1,1) = 10.0_wp ; density2(1,1) = 400.0_wp
    mass2(2,1) =  300.0_wp ; mass_w2(2,1) =  5.0_wp ; density2(2,1) =   0.0_wp
    mass2(3,1) =  -50.0_wp ; mass_w2(3,1) =  2.0_wp ; density2(3,1) = 500.0_wp
    mass2(4,1) =  999.0_wp ; mass_w2(4,1) = 999.0_wp; density2(4,1) = 999.0_wp
    n_lay(1)   = 3

    ! Column 2: a clean two-layer column. thickness = 100/400 + 300/600 = 0.75
    mass2(1,2) = 100.0_wp ; density2(1,2) = 400.0_wp
    mass2(2,2) = 300.0_wp ; density2(2,2) = 600.0_wp
    n_lay(2)   = 2

    ! Column 3: empty.
    n_lay(3)   = 0

    mass_base = [1.0_wp_acc, 2.5_wp_acc, 0.0_wp_acc]

    call summarize_domain_state(mass2,mass_w2,density2,n_lay,thk_v,wet_v,rho_v,lwm_v)

    call check_val("col 1 thickness", thk_v(1), real(thk_exp,wp), nfail)
    call check_val("col 1 bulk_density", rho_v(1), real(rho_exp,wp), nfail)
    call check_val("col 2 thickness", thk_v(2), 0.75_wp, nfail)
    call check_val("col 2 bulk_density agrees with bulk_snow_density here", &
                   rho_v(2), bulk_snow_density(mass2(:,2),density2(:,2),2), nfail)
    call check_val("col 3 (empty) thickness", thk_v(3), 0.0_wp, nfail)
    call check_val("col 3 (empty) bulk_density", rho_v(3), 0.0_wp, nfail)
    call check_val("col 1 liquid_water", lwm_v(1), 17.0_wp, nfail)

    call summarize_year_state(mass2,mass_w2,density2,n_lay,mass_base, &
                              thk_v,wet_v,rho_v,base_v)

    call check_val("year subset: thickness matches domain summary", &
                   thk_v(1), real(thk_exp,wp), nfail)
    call check_val("year subset: wet_mass matches domain summary", &
                   wet_v(1), 367.0_wp, nfail)
    call check_val("year subset: base_mass passed through from wp_acc", &
                   base_v(2), 2.5_wp, nfail)

    ! =====================================================================
    ! get_state / print_state
    ! =====================================================================

    write(*,*)
    write(*,"(a)") "--- get_state (safe division; Julia _state_dict is not) ---"

    mass        = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp ; temperature = 0.0_wp
    mass(1) =  100.0_wp ; mass_w(1) = 10.0_wp ; density(1) = 400.0_wp ; temperature(1) = 270.0_wp
    mass(2) =  300.0_wp ; mass_w(2) =  5.0_wp ; density(2) =   0.0_wp ; temperature(2) = 271.0_wp
    mass(3) =  -50.0_wp ; mass_w(3) =  2.0_wp ; density(3) = 500.0_wp ; temperature(3) = 272.0_wp

    call get_state(snap,mass,mass_w,density,temperature,3,0.75_wp,12.5_wp_acc,c)

    call check("get_state n", snap%n .eq. 3, nfail)
    call check_val("get_state total_mass (unconditional solid sum)", &
                   snap%total_mass, 350.0_wp, nfail)
    call check_val("get_state total_liquid_water", snap%total_liquid_water, 17.0_wp, nfail)
    call check_val("get_state total_wet_mass", snap%total_wet_mass, 367.0_wp, nfail)
    call check_val("get_state total_thickness matches column_summary", &
                   snap%total_thickness, real(thk_exp,wp), nfail)

    ! This is the Julia hazard: _state_dict would compute 300/0 = Inf here
    ! and total_thickness would be Inf. chion guards it.
    call check_val("layer 2 thickness guarded to 0, not Inf (Julia hazard)", &
                   snap%thickness(2), 0.0_wp, nfail)
    call check("total_thickness is finite", &
               snap%total_thickness .lt. huge(1.0_wp), nfail)
    call check_val("layer 3 (negative mass) thickness guarded to 0", &
                   snap%thickness(3), 0.0_wp, nfail)
    call check_val("layer 1 thickness", snap%thickness(1), 0.25_wp, nfail)

    call check_val("get_state surface_temperature = T(1)", &
                   snap%surface_temperature, 270.0_wp, nfail)
    call check_val("get_state albedo", snap%albedo, 0.75_wp, nfail)
    call check_val("get_state smb_ice from wp_acc", snap%smb_ice, 12.5_wp, nfail)

    call get_state(snap,mass,mass_w,density,temperature,0,0.81_wp,0.0_wp_acc,c)
    call check_val("n=0: surface_temperature falls back to T0", &
                   snap%surface_temperature, c%T0, nfail)
    call check_val("n=0: total_mass", snap%total_mass, 0.0_wp, nfail)

    write(*,*)
    write(*,"(a)") "--- print_state output (visual check) ---"
    call get_state(snap,mass,mass_w,density,temperature,3,0.75_wp,12.5_wp_acc,c)
    call print_state(snap,1)

    ! === Summary ========================================================
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP10: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP10: ", nfail, " CHECK(S) FAILED"
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

end program test_diagnostics
