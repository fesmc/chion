program test_pdd
    ! WP9 acceptance test: bulk PDD model.
    !
    ! Covers:
    !   1. mass-balance closure to 1e-10 (all terms wp_acc, so a dp check),
    !      including an explicit regression guard on defect D2 -- the reason
    !      the three-reservoir identity does NOT close.
    !   2. the PISM / Calov-Greve expectation integral against its three limits
    !      and against the Abramowitz-Stegun polynomial it replaces.
    !   3. snowpack_swe never goes negative under extreme melt.
    !   4. ice melt occurs only once the snow reservoir is exhausted within a
    !      step.
    !   5. both pdd_method branches, and the physical difference between them.
    !
    ! See docs/pdd_defects.md for what each defect check is guarding.

    use chion_defs, only : wp, wp_acc, chion_const_class, chion_forcing_class, &
                           chion_step_forcing_class, chion_const_init, &
                           chion_forcing_alloc, chion_forcing_dealloc
    use snow_pdd

    implicit none

    type(chion_const_class) :: c
    type(pdd_par_class)     :: par
    integer                 :: nfail

    nfail = 0

    call chion_const_init(c)
    call pdd_par_init(par)

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP9 acceptance test: snow_pdd"
    write(*,"(a)") "=========================================================="
    write(*,*)

    call test_parameters(par,nfail)
    call test_pism_limits(par,c,nfail)
    call test_normal_cdf(nfail)
    call test_mass_balance(par,c,nfail)
    call test_no_negative_snowpack(par,c,nfail)
    call test_ice_melt_ordering(par,c,nfail)
    call test_method_branches(par,c,nfail)
    call test_multicolumn_step(par,c,nfail)

    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP9: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP9: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    ! =====================================================================

    subroutine test_parameters(par,nfail)

        implicit none

        type(pdd_par_class), intent(IN)    :: par
        integer,             intent(INOUT) :: nfail

        write(*,"(a)") "--- parameters and defaults ---"

        call check("ddf_snow default = 3.0",  par%ddf_snow            .eq. 3.0_wp, nfail)
        call check("ddf_ice default = 8.0",   par%ddf_ice             .eq. 8.0_wp, nfail)
        call check("refreezing_fraction = 0.6", &
                                             par%refreezing_fraction .eq. 0.6_wp, nfail)
        call check("temperature_sigma = 5.0", par%temperature_sigma   .eq. 5.0_wp, nfail)
        call check("default method reproduces Julia daily path (simple)", &
                   par%pdd_method .eq. CHION_PDD_SIMPLE, nfail)

        call check("pdd_method_flag('simple')", &
                   pdd_method_flag("simple") .eq. CHION_PDD_SIMPLE, nfail)
        call check("pdd_method_flag('pism')", &
                   pdd_method_flag("pism") .eq. CHION_PDD_PISM, nfail)

        return

    end subroutine test_parameters

    ! =====================================================================

    subroutine test_pism_limits(par,c,nfail)
        ! The three limits of E[T+] = sigma*phi(z) + Tbar*Phi(z).

        implicit none

        type(pdd_par_class),     intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(pdd_par_class) :: p
        real(wp_acc) :: sigma, teff, expect
        real(wp_acc) :: pdd
        real(wp)     :: dt

        write(*,*)
        write(*,"(a)") "--- PISM / Calov-Greve expectation integral ---"

        sigma = real(par%temperature_sigma,wp_acc)

        ! Limit 1: Tbar -> +inf gives Tbar (all fluctuations are positive).
        teff = pdd_expected_positive_temperature(1000.0_wp_acc,sigma)
        call check_close("Tbar -> +inf : E[T+] -> Tbar", teff, 1000.0_wp_acc, &
                         1.0e-10_wp_acc, nfail)

        ! Limit 2: Tbar = 0 gives sigma/sqrt(2*pi) = sigma*phi(0).
        expect = sigma/sqrt(8.0_wp_acc*atan(1.0_wp_acc))
        teff   = pdd_expected_positive_temperature(0.0_wp_acc,sigma)
        call check_close("Tbar =  0    : E[T+] = sigma/sqrt(2 pi)", teff, expect, &
                         1.0e-14_wp_acc, nfail)

        ! Limit 3: Tbar -> -inf gives 0.
        teff = pdd_expected_positive_temperature(-1000.0_wp_acc,sigma)
        call check_close("Tbar -> -inf : E[T+] -> 0", teff, 0.0_wp_acc, &
                         1.0e-30_wp_acc, nfail)

        ! Monotonicity and strict positivity: E[T+] > max(Tbar,0) always,
        ! which is the whole point of the scheme -- a sub-freezing mean still
        ! produces melt.
        call check("E[T+] > 0 at Tbar = -3 degC", &
                   pdd_expected_positive_temperature(-3.0_wp_acc,sigma) .gt. 0.0_wp_acc, nfail)
        call check("E[T+] > Tbar at Tbar = +3 degC", &
                   pdd_expected_positive_temperature(3.0_wp_acc,sigma) .gt. 3.0_wp_acc, nfail)
        call check("E[T+] monotone increasing in Tbar", &
                   pdd_expected_positive_temperature(-1.0_wp_acc,sigma) .lt. &
                   pdd_expected_positive_temperature( 1.0_wp_acc,sigma), nfail)

        ! And the same limits carried through pdd_degree_days, which is where
        ! the dt_days factor enters. Monthly step, as Chion.jl would use.
        p    = par
        p%pdd_method = CHION_PDD_PISM
        dt   = 30.0_wp

        pdd    = pdd_degree_days(c%T0,dt,p,c)
        expect = real(dt,wp_acc)*sigma/sqrt(8.0_wp_acc*atan(1.0_wp_acc))
        call check_close("pdd at T = T0, dt = 30 d", pdd, expect, 1.0e-6_wp_acc, nfail)

        return

    end subroutine test_pism_limits

    ! =====================================================================

    subroutine test_normal_cdf(nfail)
        ! Document the one deliberate numerical deviation from Chion.jl: the
        ! Abramowitz-Stegun polynomial _normal_cdf (pdd.jl:16-21) is replaced
        ! by 0.5*erfc(-z/sqrt(2)). Assert they agree to the polynomial's own
        ! stated accuracy, and that the exact form is the better one in the
        ! far negative tail, where the polynomial cannot represent the value
        ! at all.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer      :: i
        real(wp_acc) :: z, dmax, d

        write(*,*)
        write(*,"(a)") "--- normal CDF: erfc form vs Abramowitz-Stegun polynomial ---"

        dmax = 0.0_wp_acc

        do i = -600, 600
            z    = real(i,wp_acc)*0.01_wp_acc
            d    = abs(pdd_normal_cdf(z) - normal_cdf_as(z))
            dmax = max(dmax,d)
        end do

        write(*,"(a,g14.6)") "         max |Phi_erfc - Phi_AS| over z in [-6,6] = ", dmax

        call check("agreement is at the polynomial's stated ~1e-7 level", &
                   dmax .lt. 1.0e-6_wp_acc, nfail)
        call check("difference is not zero (the two really are different code)", &
                   dmax .gt. 0.0_wp_acc, nfail)

        call check_close("Phi(0) = 1/2 exactly", pdd_normal_cdf(0.0_wp_acc), &
                         0.5_wp_acc, 1.0e-15_wp_acc, nfail)

        ! Far tail: erfc resolves this, the polynomial saturates at zero.
        call check("Phi(-20) is resolved, not flushed to zero", &
                   pdd_normal_cdf(-20.0_wp_acc) .gt. 0.0_wp_acc, nfail)

        return

    end subroutine test_normal_cdf

    ! =====================================================================

    subroutine test_mass_balance(par,c,nfail)
        ! Cumulative mass-balance closure over a synthetic annual cycle.
        !
        ! WHAT ACTUALLY CLOSES:
        !     d(smb_ice) + d(runoff) == snowfall + rainfall
        ! exactly, every step, because smb_ice and runoff are credited with
        ! complementary halves of every term.
        !
        ! WHAT DOES NOT CLOSE, and is the substance of defect D2:
        !     d(smb_ice) + d(runoff) == snowfall + rainfall - d(snowpack_swe)
        ! The residual of that identity is exactly d(snowpack_swe), because
        ! smb_ice is credited with the snowpack change as well as with the
        ! flux to the ice. smb_ice is documented in Chion.jl's NetCDF metadata
        ! as "Net mass forcing to the ice sheet"; it is not that. Asserting
        ! the residual EQUALS d(snowpack_swe) turns the defect into a
        ! regression guard: the check breaks the moment upstream fixes it.

        implicit none

        type(pdd_par_class),     intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: nstep = 365

        type(chion_step_forcing_class) :: forc
        type(pdd_par_class) :: p
        integer      :: n, imethod
        real(wp)     :: swe, swe0, pi_c
        real(wp_acc) :: smb, runoff, pdd_sum
        real(wp_acc) :: smb0, runoff0
        real(wp_acc) :: input_total, residual, dswe
        character(len=8) :: mname

        write(*,*)
        write(*,"(a)") "--- mass balance over a 365-step annual cycle ---"

        pi_c = 4.0_wp*atan(1.0_wp)

        do imethod = 1, 2

            p = par
            if (imethod .eq. 1) then
                p%pdd_method = CHION_PDD_SIMPLE
                mname = "simple"
            else
                p%pdd_method = CHION_PDD_PISM
                mname = "pism"
            end if

            call step_forcing_init(forc)
            forc%dt_days = 1.0_wp

            ! Start with a pre-existing snowpack, so that d(snowpack_swe) is
            ! not merely the accumulated snowfall.
            swe     = 250.0_wp
            smb     = 0.0_wp_acc
            runoff  = 0.0_wp_acc
            pdd_sum = 0.0_wp_acc

            swe0    = swe
            smb0    = smb
            runoff0 = runoff

            input_total = 0.0_wp_acc

            do n = 1, nstep

                ! Sinusoidal annual cycle, peak +8 degC, trough -18 degC.
                forc%air_temperature = c%T0 - 5.0_wp &
                    + 13.0_wp*sin(2.0_wp*pi_c*(real(n,wp)-100.0_wp)/365.0_wp)

                ! Snowfall in the cold half, rain in the warm half, with a
                ! deliberately negative rate on some steps to exercise the
                ! max(rate,0) clip.
                forc%snowfall_rate = 3.0e-5_wp*(1.0_wp &
                    - sin(2.0_wp*pi_c*(real(n,wp)-100.0_wp)/365.0_wp))
                forc%rainfall_rate = 1.0e-5_wp &
                    * sin(2.0_wp*pi_c*(real(n,wp)-100.0_wp)/365.0_wp)

                input_total = input_total &
                            + pdd_step_mass(forc%snowfall_rate,forc%dt_days,c) &
                            + pdd_step_mass(forc%rainfall_rate,forc%dt_days,c)

                call pdd_column_step(swe,smb,runoff,pdd_sum,forc,p,c)

            end do

            dswe = real(swe,wp_acc) - real(swe0,wp_acc)

            write(*,"(a,a,a)") "  method = ", trim(mname), ":"
            write(*,"(a,g16.8)") "         total input   [kg m-2] = ", input_total
            write(*,"(a,g16.8)") "         d(smb_ice)    [kg m-2] = ", smb-smb0
            write(*,"(a,g16.8)") "         d(runoff)     [kg m-2] = ", runoff-runoff0
            write(*,"(a,g16.8)") "         d(snowpack)   [kg m-2] = ", dswe
            write(*,"(a,g16.8)") "         pdd_sum       [K d]    = ", pdd_sum

            ! --- the closure that holds ---
            residual = (smb-smb0) + (runoff-runoff0) - input_total
            call check_close("closure: d(smb)+d(runoff) == snowfall+rainfall", &
                             residual, 0.0_wp_acc, &
                             1.0e-10_wp_acc*max(input_total,1.0_wp_acc), nfail)

            ! --- the closure that does not, and by exactly how much ---
            residual = (smb-smb0) + (runoff-runoff0) &
                     - (input_total - dswe)

            call check("D2 guard: three-reservoir identity is violated (non-vacuous)", &
                       abs(dswe) .gt. 1.0_wp_acc, nfail)
            call check_close("D2 guard: residual equals d(snowpack_swe) exactly", &
                             residual, dswe, &
                             1.0e-6_wp_acc*max(abs(dswe),1.0_wp_acc), nfail)

            call check("runoff is non-negative", runoff .ge. 0.0_wp_acc, nfail)
            call check("pdd_sum is non-negative", pdd_sum .ge. 0.0_wp_acc, nfail)

        end do

        return

    end subroutine test_mass_balance

    ! =====================================================================

    subroutine test_no_negative_snowpack(par,c,nfail)
        ! snowpack_swe must never go negative, however extreme the melt.
        ! Guaranteed by snow_melt = min(available_snow, ddf_snow*pdd), so the
        ! test is a structural guard on that min() rather than on physics.

        implicit none

        type(pdd_par_class),     intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_step_forcing_class) :: forc
        type(pdd_par_class) :: p
        integer      :: n
        real(wp)     :: swe, swe_min
        real(wp_acc) :: smb, runoff, pdd_sum

        write(*,*)
        write(*,"(a)") "--- snowpack_swe never negative under extreme melt ---"

        ! Case 1: refreezing_fraction = 0, so the reservoir can actually empty.
        p = par
        p%refreezing_fraction = 0.0_wp

        call step_forcing_init(forc)
        forc%dt_days         = 1.0_wp
        forc%air_temperature = c%T0 + 1000.0_wp     ! absurd, on purpose
        forc%snowfall_rate   = 0.0_wp
        forc%rainfall_rate   = 0.0_wp

        swe     = 50.0_wp
        smb     = 0.0_wp_acc
        runoff  = 0.0_wp_acc
        pdd_sum = 0.0_wp_acc

        swe_min = swe

        do n = 1, 20
            call pdd_column_step(swe,smb,runoff,pdd_sum,forc,p,c)
            swe_min = min(swe_min,swe)
        end do

        call check("f_refrz = 0: reservoir empties exactly to zero", &
                   swe .eq. 0.0_wp, nfail)
        call check("f_refrz = 0: never negative", swe_min .ge. 0.0_wp, nfail)

        ! Case 2: default refreezing_fraction = 0.6. The reservoir decays
        ! geometrically by (1-f) each step and is NEVER exhausted -- defect D6.
        p = par

        swe     = 50.0_wp
        smb     = 0.0_wp_acc
        runoff  = 0.0_wp_acc
        pdd_sum = 0.0_wp_acc

        call pdd_column_step(swe,smb,runoff,pdd_sum,forc,p,c)
        call check_val("f_refrz = 0.6: one step of infinite melt leaves 0.6*swe", &
                       swe, 30.0_wp, nfail)

        swe_min = swe
        do n = 1, 20
            call pdd_column_step(swe,smb,runoff,pdd_sum,forc,p,c)
            swe_min = min(swe_min,swe)
        end do

        call check("f_refrz = 0.6: never negative", swe_min .ge. 0.0_wp, nfail)
        call check("D6 guard: reservoir still non-zero after 21 steps of &
                   &infinite melt", swe .gt. 0.0_wp, nfail)

        ! Case 3: negative precipitation rates are clipped, not subtracted.
        swe     = 10.0_wp
        smb     = 0.0_wp_acc
        runoff  = 0.0_wp_acc
        pdd_sum = 0.0_wp_acc

        forc%air_temperature = c%T0 - 10.0_wp
        forc%snowfall_rate   = -1.0e-3_wp
        forc%rainfall_rate   = -1.0e-3_wp

        call pdd_column_step(swe,smb,runoff,pdd_sum,forc,p,c)

        call check_val("negative snowfall rate is clipped to zero", swe, 10.0_wp, nfail)
        call check_close("negative rainfall rate contributes no runoff", &
                         runoff, 0.0_wp_acc, 1.0e-12_wp_acc, nfail)

        return

    end subroutine test_no_negative_snowpack

    ! =====================================================================

    subroutine test_ice_melt_ordering(par,c,nfail)
        ! Ice melt may only draw on degree days LEFT OVER after the snow
        ! reservoir is exhausted within the step.
        !
        ! ice_melt is not stored, but it is recoverable exactly:
        !     d(snowpack) = snowfall - snow_melt + refrozen
        !     d(smb_ice)  = snowfall - snow_melt + refrozen - ice_melt
        ! so  ice_melt   = d(snowpack) - d(smb_ice).

        implicit none

        type(pdd_par_class),     intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_step_forcing_class) :: forc
        real(wp)     :: swe, swe0
        real(wp_acc) :: smb, runoff, pdd_sum, ice_melt

        write(*,*)
        write(*,"(a)") "--- ice melt only after the snow reservoir is exhausted ---"

        call step_forcing_init(forc)
        forc%dt_days       = 1.0_wp
        forc%rainfall_rate = 0.0_wp
        ! 6 kg m-2 of snowfall in the step: 6/86400 kg m-2 s-1.
        forc%snowfall_rate = 6.0_wp/real(c%seconds_per_day,wp)

        ! available_snow = 90 + 6 = 96; ddf_snow = 3, so 32 K d exhausts it.

        ! (a) pdd = 1: deep in the snow-limited regime.
        ! ice_melt is recovered through snowpack_swe, which is stored in sp,
        ! so the recovery carries sp round-off (~1e-5 of the reservoir). The
        ! zero-ice-melt checks therefore use a small absolute tolerance rather
        ! than an exact comparison; 1e-3 kg m-2 is four orders of magnitude
        ! below the ~8 kg m-2 that a single leftover degree day would produce.
        forc%air_temperature = c%T0 + 1.0_wp
        call run_one(swe0,swe,smb,runoff,pdd_sum,ice_melt,90.0_wp,forc,par,c)
        call check_close("pdd =  1 K d (potential  3 << 96): no ice melt", &
                         ice_melt, 0.0_wp_acc, 1.0e-3_wp_acc, nfail)
        call check("pdd =  1 K d: snowpack still present", swe .gt. 0.0_wp, nfail)

        ! (b) pdd = 31: potential snow melt 93 < 96, still snow-limited.
        forc%air_temperature = c%T0 + 31.0_wp
        call run_one(swe0,swe,smb,runoff,pdd_sum,ice_melt,90.0_wp,forc,par,c)
        call check_close("pdd = 31 K d (potential 93 <  96): still no ice melt", &
                         ice_melt, 0.0_wp_acc, 1.0e-3_wp_acc, nfail)

        ! (c) pdd = 33: potential 99 > 96, so 1 K d is left for ice.
        forc%air_temperature = c%T0 + 33.0_wp
        call run_one(swe0,swe,smb,runoff,pdd_sum,ice_melt,90.0_wp,forc,par,c)
        call check("pdd = 33 K d (potential 99 >  96): ice melt occurs", &
                   ice_melt .gt. 0.0_wp_acc, nfail)
        call check_close("ice melt = ddf_ice*(pdd - available/ddf_snow) = 8*1", &
                         ice_melt, 8.0_wp_acc, 0.05_wp_acc, nfail)

        ! (d) bare column, no snowfall: all degree days go to ice.
        forc%snowfall_rate   = 0.0_wp
        forc%air_temperature = c%T0 + 10.0_wp
        call run_one(swe0,swe,smb,runoff,pdd_sum,ice_melt,0.0_wp,forc,par,c)
        call check_close("bare column: ice melt = ddf_ice*pdd = 80", &
                         ice_melt, 80.0_wp_acc, 0.05_wp_acc, nfail)
        call check("bare column: snowpack stays at zero", swe .eq. 0.0_wp, nfail)
        call check_close("bare column: smb_ice = -ice_melt", smb, -ice_melt, &
                         1.0e-12_wp_acc, nfail)
        call check_close("bare column: runoff = +ice_melt", runoff, ice_melt, &
                         1.0e-12_wp_acc, nfail)

        ! (e) rainfall never touches the ice budget: it is pure runoff.
        forc%rainfall_rate   = 20.0_wp/real(c%seconds_per_day,wp)
        forc%air_temperature = c%T0 - 10.0_wp
        call run_one(swe0,swe,smb,runoff,pdd_sum,ice_melt,100.0_wp,forc,par,c)
        call check_close("rainfall goes entirely to runoff", runoff, &
                         20.0_wp_acc, 1.0e-4_wp_acc, nfail)
        call check_close("rainfall does not enter smb_ice", smb, 0.0_wp_acc, &
                         1.0e-4_wp_acc, nfail)
        call check_val("rainfall does not enter the snowpack", swe, 100.0_wp, nfail)

        return

    end subroutine test_ice_melt_ordering

    ! =====================================================================

    subroutine test_method_branches(par,c,nfail)
        ! Both pdd_method branches, and the physical difference between them:
        ! at a sub-freezing mean temperature the simple scheme produces NO
        ! melt while the expectation integral produces a real amount. smbpal
        ! always uses the integral (smbpal.f90:300, :383); Chion.jl uses it
        ! only for 27-32 day steps. Defect D4.

        implicit none

        type(pdd_par_class),     intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_step_forcing_class) :: forc
        type(pdd_par_class) :: p_simple, p_pism
        real(wp)     :: swe_s, swe_p
        real(wp_acc) :: smb_s, runoff_s, pdd_s
        real(wp_acc) :: smb_p, runoff_p, pdd_p

        write(*,*)
        write(*,"(a)") "--- pdd_method branches ---"

        p_simple = par ; p_simple%pdd_method = CHION_PDD_SIMPLE
        p_pism   = par ; p_pism%pdd_method   = CHION_PDD_PISM

        call step_forcing_init(forc)
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 0.0_wp
        forc%rainfall_rate   = 0.0_wp
        forc%air_temperature = c%T0 - 5.0_wp        ! sub-freezing mean

        swe_s = 100.0_wp ; smb_s = 0.0_wp_acc ; runoff_s = 0.0_wp_acc ; pdd_s = 0.0_wp_acc
        swe_p = 100.0_wp ; smb_p = 0.0_wp_acc ; runoff_p = 0.0_wp_acc ; pdd_p = 0.0_wp_acc

        call pdd_column_step(swe_s,smb_s,runoff_s,pdd_s,forc,p_simple,c)
        call pdd_column_step(swe_p,smb_p,runoff_p,pdd_p,forc,p_pism,c)

        write(*,"(a,g16.8)") "         T = T0-5, simple pdd [K d] = ", pdd_s
        write(*,"(a,g16.8)") "         T = T0-5, pism   pdd [K d] = ", pdd_p

        call check("simple: no melt below freezing", pdd_s .eq. 0.0_wp_acc, nfail)
        call check("simple: snowpack unchanged", swe_s .eq. 100.0_wp, nfail)
        call check("pism: sub-freezing mean still melts", pdd_p .gt. 0.0_wp_acc, nfail)
        call check("pism: snowpack decreases", swe_p .lt. 100.0_wp, nfail)

        ! sigma = 5, Tbar = -5 -> z = -1: E[T+] = 5*phi(-1) - 5*Phi(-1)
        !                                = 1.209854 - 0.793274 = 0.416580
        call check_close("pism pdd at Tbar = -5, sigma = 5", pdd_p, &
                         0.41658_wp_acc, 1.0e-4_wp_acc, nfail)

        ! Well above freezing the two converge: E[T+] -> Tbar.
        forc%air_temperature = c%T0 + 40.0_wp
        swe_s = 0.0_wp ; smb_s = 0.0_wp_acc ; runoff_s = 0.0_wp_acc ; pdd_s = 0.0_wp_acc
        swe_p = 0.0_wp ; smb_p = 0.0_wp_acc ; runoff_p = 0.0_wp_acc ; pdd_p = 0.0_wp_acc

        call pdd_column_step(swe_s,smb_s,runoff_s,pdd_s,forc,p_simple,c)
        call pdd_column_step(swe_p,smb_p,runoff_p,pdd_p,forc,p_pism,c)

        call check_close("the two methods converge far above freezing", &
                         pdd_p, pdd_s, 1.0e-6_wp_acc, nfail)

        return

    end subroutine test_method_branches

    ! =====================================================================

    subroutine test_multicolumn_step(par,c,nfail)
        ! pdd_step over a column list, and the active-mask fix (defect D8):
        ! columns absent from active_idx must not advance.

        implicit none

        type(pdd_par_class),     intent(IN)    :: par
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: ncol = 4

        type(pdd_class)           :: pdd
        type(chion_forcing_class) :: forc
        integer :: active_idx(3)
        integer :: i

        write(*,*)
        write(*,"(a)") "--- pdd_step over a column list ---"

        pdd%par = par
        call pdd_alloc(pdd,ncol)

        call check("cold start: snowpack_swe = 0", &
                   all(pdd%now%snowpack_swe .eq. 0.0_wp), nfail)
        call check("cold start: accumulators = 0", &
                   all(pdd%now%smb_ice .eq. 0.0_wp_acc) .and. &
                   all(pdd%now%runoff  .eq. 0.0_wp_acc), nfail)

        call chion_forcing_alloc(forc,ncol)

        forc%air_temperature = c%T0 + 5.0_wp
        forc%snowfall_rate   = 10.0_wp/real(c%seconds_per_day,wp)
        forc%rainfall_rate   = 0.0_wp

        ! Column 3 is deactivated.
        active_idx = [1,2,4]

        do i = 1, 10
            call pdd_step(pdd,forc,1.0_wp,active_idx,c)
        end do

        call check("D8 guard: inactive column did not advance", &
                   pdd%now%pdd_sum(3) .eq. 0.0_wp_acc .and. &
                   pdd%now%snowpack_swe(3) .eq. 0.0_wp, nfail)
        call check("active columns advanced", &
                   pdd%now%pdd_sum(1) .gt. 0.0_wp_acc, nfail)
        call check("identical forcing gives identical columns", &
                   pdd%now%snowpack_swe(1) .eq. pdd%now%snowpack_swe(2) .and. &
                   pdd%now%smb_ice(1)      .eq. pdd%now%smb_ice(4), nfail)

        call pdd_reset_column(pdd,1)
        call check("pdd_reset_column zeroes one column only", &
                   pdd%now%snowpack_swe(1) .eq. 0.0_wp .and. &
                   pdd%now%snowpack_swe(2) .gt. 0.0_wp, nfail)

        call chion_forcing_dealloc(forc)
        call pdd_dealloc(pdd)

        call check("pdd_dealloc releases state", &
                   .not. allocated(pdd%now%snowpack_swe), nfail)

        return

    end subroutine test_multicolumn_step

    ! =====================================================================
    ! Helpers
    ! =====================================================================

    subroutine run_one(swe0,swe,smb,runoff,pdd_sum,ice_melt,swe_init,forc,par,c)
        ! One step from a clean start, returning the recovered ice melt.

        implicit none

        real(wp),                       intent(OUT) :: swe0
        real(wp),                       intent(OUT) :: swe
        real(wp_acc),                   intent(OUT) :: smb
        real(wp_acc),                   intent(OUT) :: runoff
        real(wp_acc),                   intent(OUT) :: pdd_sum
        real(wp_acc),                   intent(OUT) :: ice_melt
        real(wp),                       intent(IN)  :: swe_init
        type(chion_step_forcing_class), intent(IN)  :: forc
        type(pdd_par_class),            intent(IN)  :: par
        type(chion_const_class),        intent(IN)  :: c

        swe0    = swe_init
        swe     = swe_init
        smb     = 0.0_wp_acc
        runoff  = 0.0_wp_acc
        pdd_sum = 0.0_wp_acc

        call pdd_column_step(swe,smb,runoff,pdd_sum,forc,par,c)

        ice_melt = (real(swe,wp_acc) - real(swe0,wp_acc)) - smb

        return

    end subroutine run_one

    subroutine step_forcing_init(forc)
        ! Neutral per-column forcing. PDD reads only air_temperature, dt_days,
        ! snowfall_rate and rainfall_rate, but the whole type is filled so the
        ! test never depends on uninitialized memory.

        implicit none

        type(chion_step_forcing_class), intent(OUT) :: forc

        forc%air_temperature = 273.15_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 0.0_wp
        forc%rainfall_rate   = 0.0_wp
        forc%shortwave_down  = 0.0_wp
        forc%wind_speed      = 0.0_wp

        forc%q_sw_net  = 0.0_wp
        forc%q_lw_down = 0.0_wp
        forc%q_sh      = 0.0_wp
        forc%q_lh      = 0.0_wp

        forc%has_q_sw_net  = .FALSE.
        forc%has_q_lw_down = .FALSE.
        forc%has_q_sh      = .FALSE.
        forc%has_q_lh      = .FALSE.

        forc%relative_humidity     = 0.0_wp
        forc%has_relative_humidity = .FALSE.

        forc%air_pressure          = 101325.0_wp
        forc%prescribed_albedo     = 0.0_wp
        forc%has_prescribed_albedo = .FALSE.

        forc%latitude_deg        = 0.0_wp
        forc%day_of_year         = 1.0_wp
        forc%solar_longitude_deg = 0.0_wp

        return

    end subroutine step_forcing_init

    pure function normal_cdf_as(x) result(cdf)
        ! Chion.jl _normal_cdf (src/processes/pdd.jl:16-21), the
        ! Abramowitz-Stegun rational-polynomial approximation, reproduced here
        ! ONLY so the test can quantify what replacing it costs. It is not
        ! used by snow_pdd.

        implicit none

        real(wp_acc), intent(IN) :: x
        real(wp_acc) :: cdf

        ! Local variables
        real(wp_acc) :: t, poly, pdf

        t    = 1.0_wp_acc/(1.0_wp_acc + 0.2316419_wp_acc*abs(x))
        poly = t*(0.319381530_wp_acc &
             + t*(-0.356563782_wp_acc &
             + t*( 1.781477937_wp_acc &
             + t*(-1.821255978_wp_acc &
             + t*  1.330274429_wp_acc))))
        pdf  = 0.398942280401432678_wp_acc*exp(-0.5_wp_acc*x*x)

        cdf = 1.0_wp_acc - pdf*poly
        if (x .lt. 0.0_wp_acc) cdf = 1.0_wp_acc - cdf

        return

    end function normal_cdf_as

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

    subroutine check_close(label,value,expected,tol,nfail)
        ! Absolute-tolerance dp comparison. The caller states the tolerance
        ! explicitly because WP9 mixes wp_acc closure checks (1e-10) with
        ! sp-limited state checks (1e-6 relative).

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        real(wp_acc),     intent(IN)    :: tol
        integer,          intent(INOUT) :: nfail

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g16.8)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_close

end program test_pdd
