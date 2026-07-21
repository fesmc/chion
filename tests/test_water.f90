program test_water
    ! WP6 acceptance test: percolation, refreezing and melt.
    !
    ! Covers:
    !   percolation -- water conservation, the massless-layer-to-runoff path,
    !                  the collapsed-pore path, and retention at exactly
    !                  max_lwc;
    !   refreezing  -- energy closure refrozen*Lm == sum of consumed cold
    !                  content, both the partial and the complete branch, and
    !                  the rho_i density cap including its deliberate
    !                  volume non-conservation;
    !   melt        -- fast path, and mass conservation through layer
    !                  depletion.

    use chion_defs,        only : wp, wp_acc, TOL_TINY, &
                                  chion_const_class, chion_const_init
    use snow_column_utils, only : layer_pore_volume, layer_lwc
    use snow_percolation,  only : apply_percolation, PERCOLATION_MAX_LWC_DEFAULT
    use snow_refreezing,   only : apply_refreezing
    use snow_melt,         only : apply_melt

    implicit none

    integer, parameter :: Ntot = 5

    type(chion_const_class) :: c

    real(wp) :: mass(Ntot), mass_w(Ntot), density(Ntot), temperature(Ntot)
    real(wp) :: t_srf, albedo_dyn
    real(wp) :: rho_i, rho_w, T0, ci, Lm

    real(wp_acc) :: runoff, refrozen, melted
    real(wp_acc) :: before, after, expected, phi
    integer      :: n, nfail, k

    call chion_const_init(c)

    rho_i = c%rho_i
    rho_w = c%rho_w
    T0    = c%T0
    ci    = c%ci
    Lm    = c%Lm

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP6 acceptance test: percolation / refreezing / melt"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! ======================================================================
    ! PERCOLATION
    ! ======================================================================
    write(*,"(a)") "--- percolation: water conservation ---"

    n           = 5
    mass        = [  50.0_wp, 120.0_wp, 200.0_wp, 300.0_wp, 400.0_wp ]
    density     = [ 200.0_wp, 300.0_wp, 450.0_wp, 600.0_wp, 880.0_wp ]
    mass_w      = [  30.0_wp,  10.0_wp,   0.0_wp,   5.0_wp,   2.0_wp ]
    temperature = T0

    before = sum(real(mass_w(1:n),wp_acc))

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff)

    after = sum(real(mass_w(1:n),wp_acc))

    call check_rel("sum(mass_w) before == sum(mass_w) after + runoff", &
                   after + runoff, before, nfail)

    call check("runoff is strictly positive for a saturated column", &
               runoff .gt. 0.0_wp_acc, nfail)

    ! Every layer must now be at or below max_lwc (single sweep suffices
    ! because excess is deposited before the loop reaches the layer below).
    do k = 1, n
        call check("layer is at or below max_lwc after the sweep", &
                   layer_lwc(mass(k),mass_w(k),density(k),rho_i,rho_w) &
                   .le. real(PERCOLATION_MAX_LWC_DEFAULT,wp_acc) &
                        *(1.0_wp_acc + 1.0e-6_wp_acc), nfail)
    end do

    call check("percolation did not touch mass", &
               abs(mass(1)-50.0_wp) .le. 0.0_wp .and. &
               abs(mass(5)-400.0_wp) .le. 0.0_wp, nfail)

    write(*,*)
    write(*,"(a)") "--- percolation: massless layer routes STRAIGHT to runoff ---"

    ! Layer 2 has no solid mass but does hold water, and layer 3 below it is
    ! a perfectly good snow layer. The Julia CODE sends that water to runoff
    ! anyway; the Julia DOC page says it goes to layer 3. This check is what
    ! distinguishes the two.
    n           = 3
    mass        = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1)     = 100.0_wp ; density(1) = 400.0_wp
    mass(2)     =   0.0_wp ; density(2) = 400.0_wp ; mass_w(2) = 5.0_wp
    mass(3)     = 100.0_wp ; density(3) = 400.0_wp
    temperature = T0

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff)

    call check_rel("massless layer water becomes runoff", runoff, 5.0_wp_acc, nfail)
    call check("massless layer is emptied", mass_w(2) .eq. 0.0_wp, nfail)
    call check("water did NOT reach the layer below (code, not doc)", &
               mass_w(3) .eq. 0.0_wp, nfail)

    write(*,*)
    write(*,"(a)") "--- percolation: collapsed pore space pushes water down ---"

    n           = 2
    mass        = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1)     = 100.0_wp ; density(1) = rho_i ; mass_w(1) = 5.0_wp
    mass(2)     = 100.0_wp ; density(2) = 400.0_wp
    temperature = T0

    call check("pore volume of the top layer is at the collapse threshold", &
               layer_pore_volume(mass(1),density(1),rho_i) .le. TOL_TINY, nfail)

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff)

    call check("collapsed layer is emptied", mass_w(1) .eq. 0.0_wp, nfail)
    call check_rel("all its water is now in the layer below", &
                   real(mass_w(2),wp_acc), 5.0_wp_acc, nfail)
    call check("nothing left the column", runoff .eq. 0.0_wp_acc, nfail)

    ! Same layer, but now it is the lowest one: the water must run off.
    n           = 1
    mass        = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1)     = 100.0_wp ; density(1) = rho_i ; mass_w(1) = 5.0_wp

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff)

    call check_rel("collapsed lowest layer runs off instead", runoff, 5.0_wp_acc, nfail)

    write(*,*)
    write(*,"(a)") "--- percolation: retention is exactly max_lwc ---"

    ! Single layer, heavily over-saturated. Whatever is left must sit exactly
    ! at the retention capacity 0.1*rho_w*phi.
    n          = 1
    mass       = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1)    = 300.0_wp
    density(1) = 400.0_wp
    mass_w(1)  = 100.0_wp

    phi      = layer_pore_volume(mass(1),density(1),rho_i)
    expected = real(PERCOLATION_MAX_LWC_DEFAULT,wp_acc)*phi*real(rho_w,wp_acc)

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff)

    call check_rel("retained water == max_lwc*rho_w*phi", &
                   real(mass_w(1),wp_acc), expected, nfail)
    call check_rel("runoff == the excess", runoff, 100.0_wp_acc - expected, nfail)
    call check_rel("resulting LWC == max_lwc", &
                   layer_lwc(mass(1),mass_w(1),density(1),rho_i,rho_w), &
                   real(PERCOLATION_MAX_LWC_DEFAULT,wp_acc), nfail)

    ! A layer sitting exactly at capacity must not lose anything measurable.
    n          = 1
    mass_w(1)  = real(expected,wp)

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff)

    call check("a layer already at capacity sheds nothing (to sp round-off)", &
               runoff .le. 1.0e-6_wp_acc*expected, nfail)

    ! max_lwc is an independent parameter, NOT c%max_lwc_albedo. Doubling it
    ! must double the retention. See docs/PLAN.md section 5, item 4.
    n          = 1
    mass_w(1)  = 100.0_wp

    call apply_percolation(mass,mass_w,density,n,rho_i,rho_w,runoff,max_lwc=0.2_wp)

    call check_rel("max_lwc is honoured as an independent argument", &
                   real(mass_w(1),wp_acc), 2.0_wp_acc*expected, nfail)

    ! ======================================================================
    ! REFREEZING
    ! ======================================================================
    write(*,*)
    write(*,"(a)") "--- refreezing: partial branch (cold content limits) ---"

    n              = 1
    mass           = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    temperature    = T0
    mass(1)        = 100.0_wp
    density(1)     = 400.0_wp
    mass_w(1)      = 10.0_wp
    temperature(1) = 272.15_wp

    ! Qcold = (T0-T)*ci*m_s, evaluated from the stored sp values in dp.
    expected = (real(T0,wp_acc)-real(temperature(1),wp_acc)) &
               *real(ci,wp_acc)*real(mass(1),wp_acc)/real(Lm,wp_acc)

    call apply_refreezing(mass,mass_w,density,temperature,n,T0,ci,Lm,rho_i,refrozen)

    call check_rel("refrozen == Qcold/Lm", refrozen, expected, nfail)
    call check("temperature is set to T0 exactly", temperature(1) .eq. T0, nfail)
    call check_rel("mass gained the refrozen increment", &
                   real(mass(1),wp_acc), 100.0_wp_acc + expected, nfail)
    call check_rel("liquid water lost the same increment", &
                   real(mass_w(1),wp_acc), 10.0_wp_acc - expected, nfail)
    call check("liquid water remains (partial, not complete)", &
               mass_w(1) .gt. 0.0_wp, nfail)
    call check_rel("density scaled by the mass ratio", &
                   real(density(1),wp_acc), &
                   400.0_wp_acc*(100.0_wp_acc+expected)/100.0_wp_acc, nfail)

    write(*,*)
    write(*,"(a)") "--- refreezing: complete branch (water limits) ---"

    n              = 1
    mass           = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    temperature    = T0
    mass(1)        = 1000.0_wp
    density(1)     = 400.0_wp
    mass_w(1)      = 1.0_wp
    temperature(1) = 272.15_wp

    expected = (real(mass_w(1),wp_acc)*real(Lm,wp_acc)/real(ci,wp_acc)   &
                + real(mass_w(1),wp_acc)*real(T0,wp_acc)                 &
                + real(temperature(1),wp_acc)*real(mass(1),wp_acc))      &
               /(real(mass_w(1),wp_acc) + real(mass(1),wp_acc))

    call apply_refreezing(mass,mass_w,density,temperature,n,T0,ci,Lm,rho_i,refrozen)

    call check_rel("refrozen == all the liquid water", refrozen, 1.0_wp_acc, nfail)
    call check("liquid water is exactly zero", mass_w(1) .eq. 0.0_wp, nfail)
    call check_rel("mass gained the full water mass", &
                   real(mass(1),wp_acc), 1001.0_wp_acc, nfail)
    call check_rel("temperature follows the mixing formula", &
                   real(temperature(1),wp_acc), expected, nfail)
    call check("layer stays below the melting point", temperature(1) .lt. T0, nfail)

    write(*,*)
    write(*,"(a)") "--- refreezing: rho_i cap, and its deliberate volume loss ---"

    ! Trap 7 (docs/PLAN.md section 5): density is capped at rho_i while mass
    ! still gains the FULL increment, so mass is conserved and volume is not.
    n              = 1
    mass           = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    temperature    = T0
    mass(1)        = 100.0_wp
    density(1)     = 900.0_wp
    mass_w(1)      = 10.0_wp
    temperature(1) = 250.0_wp

    call apply_refreezing(mass,mass_w,density,temperature,n,T0,ci,Lm,rho_i,refrozen)

    call check_rel("refrozen == all the water (cold content is ample)", &
                   refrozen, 10.0_wp_acc, nfail)
    call check("density is capped at rho_i", density(1) .eq. rho_i, nfail)
    call check_rel("mass still gained the FULL increment (mass conserved)", &
                   real(mass(1),wp_acc), 110.0_wp_acc, nfail)
    ! Uncapped, the new density would have been 900*110/100 = 990, leaving the
    ! layer thickness unchanged at 0.111111 m. With the cap it becomes
    ! 110/917 = 0.119956 m: volume is NOT conserved. Assert the discrepancy so
    ! that anyone "fixing" it breaks this test on purpose.
    call check("volume is deliberately NOT conserved by the cap", &
               abs(real(mass(1),wp_acc)/real(density(1),wp_acc) - 100.0_wp_acc/900.0_wp_acc) &
               .gt. 1.0e-3_wp_acc, nfail)

    write(*,*)
    write(*,"(a)") "--- refreezing: energy closure over a mixed column ---"

    ! Layers 1 and 3 refreeze partially, layer 2 completely, layer 4 is above
    ! the melting point (no entry), layer 5 has no liquid water.
    n           = 5
    mass        = [ 100.0_wp,  800.0_wp,  150.0_wp,  200.0_wp,  250.0_wp ]
    density     = [ 350.0_wp,  400.0_wp,  450.0_wp,  500.0_wp,  550.0_wp ]
    mass_w      = [  20.0_wp,    1.0_wp,   15.0_wp,    4.0_wp,    0.0_wp ]
    temperature = [ 272.15_wp, 272.15_wp, 271.15_wp, 273.15_wp, 260.0_wp ]

    ! Consumed cold content per layer = min(Qcold, m_w*Lm), zero if the layer
    ! does not satisfy the strict triple entry condition.
    expected = 0.0_wp_acc
    do k = 1, n
        if (mass(k) .gt. 0.0_wp .and. mass_w(k) .gt. 0.0_wp &
            .and. temperature(k) .lt. T0) then
            expected = expected + min( (real(T0,wp_acc)-real(temperature(k),wp_acc)) &
                                       *real(ci,wp_acc)*real(mass(k),wp_acc),        &
                                       real(mass_w(k),wp_acc)*real(Lm,wp_acc) )
        end if
    end do

    before = sum(real(mass(1:n),wp_acc)) + sum(real(mass_w(1:n),wp_acc))

    call apply_refreezing(mass,mass_w,density,temperature,n,T0,ci,Lm,rho_i,refrozen)

    after = sum(real(mass(1:n),wp_acc)) + sum(real(mass_w(1:n),wp_acc))

    ! Cold content is evaluated in dp from sp inputs, so this closes far
    ! tighter than the sp-limited 1e-6 used elsewhere.
    call check_tol("refrozen*Lm == sum of consumed cold content", &
                   refrozen*real(Lm,wp_acc), expected, 1.0e-12_wp_acc, nfail)

    call check_rel("refreezing moves mass but does not create it", after, before, nfail)

    call check("layer above the melting point did not refreeze", &
               mass_w(4) .eq. 4.0_wp, nfail)
    call check("dry cold layer is untouched", &
               temperature(5) .eq. 260.0_wp .and. mass(5) .eq. 250.0_wp, nfail)

    ! ======================================================================
    ! MELT
    ! ======================================================================
    write(*,*)
    write(*,"(a)") "--- melt: fast path ---"

    n           = 3
    mass        = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1:3)   = [ 200.0_wp, 300.0_wp, 300.0_wp ]
    density(1:3)= [ 350.0_wp, 400.0_wp, 450.0_wp ]
    temperature = T0
    runoff      = 0.0_wp_acc
    t_srf       = T0
    albedo_dyn  = c%alpha_dry

    call apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo_dyn, &
                    300.0_wp,100.0_wp,50.0_wp,c,melted)

    call check_rel("fast path melts the full requested mass", melted, 50.0_wp_acc, nfail)
    call check_rel("surface solid mass reduced", real(mass(1),wp_acc), 150.0_wp_acc, nfail)
    call check_rel("surface liquid water increased", real(mass_w(1),wp_acc), 50.0_wp_acc, nfail)
    call check("layer count unchanged", n .eq. 3, nfail)
    call check("no runoff on the fast path", runoff .eq. 0.0_wp_acc, nfail)

    write(*,*)
    write(*,"(a)") "--- melt: mass conservation through layer depletion ---"

    n            = 3
    mass         = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1:3)    = [ 200.0_wp, 300.0_wp, 300.0_wp ]
    mass_w(1:3)  = [   5.0_wp,   3.0_wp,   1.0_wp ]
    density(1:3) = [ 350.0_wp, 400.0_wp, 450.0_wp ]
    temperature  = T0
    runoff       = 0.0_wp_acc
    t_srf        = T0
    albedo_dyn   = c%alpha_dry

    before = sum(real(mass(1:n),wp_acc)) + sum(real(mass_w(1:n),wp_acc))

    call apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo_dyn, &
                    300.0_wp,100.0_wp,450.0_wp,c,melted)

    after = 0.0_wp_acc
    if (n .gt. 0) after = sum(real(mass(1:n),wp_acc)) + sum(real(mass_w(1:n),wp_acc))
    after = after + runoff

    call check_rel("total water is conserved through melt", after, before, nfail)
    call check_rel("melted mass equals the requested melt", melted, 450.0_wp_acc, nfail)
    call check("column still has layers", n .gt. 0, nfail)

    write(*,*)
    write(*,"(a)") "--- melt: column melted away entirely ---"

    n            = 2
    mass         = 0.0_wp ; mass_w = 0.0_wp ; density = 0.0_wp
    mass(1:2)    = [ 100.0_wp, 100.0_wp ]
    mass_w(1:2)  = [   2.0_wp,   2.0_wp ]
    density(1:2) = [ 350.0_wp, 400.0_wp ]
    temperature  = T0
    runoff       = 0.0_wp_acc
    t_srf        = 260.0_wp
    albedo_dyn   = c%alpha_dry

    before = sum(real(mass(1:n),wp_acc)) + sum(real(mass_w(1:n),wp_acc))

    call apply_melt(mass,mass_w,density,temperature,n,runoff,t_srf,albedo_dyn, &
                    300.0_wp,100.0_wp,500.0_wp,c,melted)

    call check("all layers are gone", n .eq. 0, nfail)
    call check_rel("melted only what was there", melted, 200.0_wp_acc, nfail)
    call check_rel("everything left the column as runoff", runoff, before, nfail)
    call check("t_srf reset to T0", t_srf .eq. T0, nfail)
    call check("albedo reset to bare ice", albedo_dyn .eq. c%alpha_ice, nfail)

    ! === Summary ==========================================================
    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP6: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP6: ", nfail, " CHECK(S) FAILED"
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

    subroutine check_rel(label,value,expected,nfail)
        ! Relative check at 1e-6, the sp-limited tolerance mandated by
        ! docs/PLAN.md section 3.1.

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        call check_tol(label,value,expected,1.0e-6_wp_acc,nfail)

        return

    end subroutine check_rel

    subroutine check_tol(label,value,expected,rtol,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        real(wp_acc),     intent(IN)    :: rtol
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp_acc) :: tol

        tol = max(abs(expected)*rtol,1.0e-30_wp_acc)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g16.8)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_tol

end program test_water
