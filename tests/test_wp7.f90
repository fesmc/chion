program test_wp7
    ! WP7 acceptance test: albedo, densification, diurnal shortwave, and
    ! fresh-snow density.
    !
    ! Coverage (docs/PLAN.md WP7 acceptance list):
    !   (i)   fresh-snow density clamps at BOTH ends, BOTH schemes
    !   (ii)  albedo stays within [alpha_wet, alpha_dry] under extreme inputs;
    !         aging is monotone non-brightening; snowfall brightening saturates
    !         at alpha_dry; the constant scheme is memoryless; a single snowfall
    !         event cannot brighten by more than alpha_dry - alpha_wet
    !   (iii) densification is monotone non-decreasing and capped at rho_i, in
    !         all three density regimes and under both low-density schemes
    !   (iv)  diurnal interval averages over a full [-pi,pi] tiling recover the
    !         daily mean to 1e-6 relative; polar day and polar night both give
    !         sensible values; substep_count returns only 1 or max_substeps
    !   (v)   the elevation-dependent diurnal T amplitude clamps at 0 and A_max,
    !         and a missing (NaN) surface height adds no excess (Chion.jl d0146e1)
    !
    ! apply_accumulation itself is exercised by WP4's and WP8's tests, since it
    ! is mostly a driver for the layer-structure routines; what is tested here
    ! is the part WP7 owns, the fresh-snow density and the density mixing it
    ! feeds.

    use chion_defs, only : wp, wp_acc, TOL_TINY, chion_const_class, &
                           chion_const_init, ALBEDO_AGING_BINARY_REFRESH, &
                           CHION_ALBEDO_CONSTANT, CHION_ALBEDO_DYNAMIC, &
                           CHION_ALBEDO_PRESCRIBED, CHION_ALBEDO_AGING, &
                           CHION_FRESH_SNOW_DENSITY_CONSTANT, &
                           CHION_FRESH_SNOW_DENSITY_PARAMETERIZED, &
                           CHION_DENSIFY_BESSI, CHION_DENSIFY_HTESSEL
    use snow_albedo
    use snow_densify
    use snow_diurnal
    use snow_accumulation, only : fresh_snow_density, FRESH_SNOW_DENSITY_MIN

    implicit none

    integer, parameter :: Ntot = 5

    type(chion_const_class) :: c
    integer                 :: nfail

    nfail = 0

    call chion_const_init(c)

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP7 acceptance test: albedo / densify / diurnal /"
    write(*,"(a)") "                            fresh-snow density"
    write(*,"(a)") "=========================================================="
    write(*,*)

    call test_fresh_snow_density(c,nfail)
    call test_albedo(c,nfail)
    call test_albedo_aging(c,nfail)
    call test_albedo_aging_rejuvenate(c,nfail)
    call test_densification(c,nfail)
    call test_diurnal(nfail)

    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP7: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP7: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    ! ======================================================================
    ! (i) fresh-snow density
    ! ======================================================================

    subroutine test_fresh_snow_density(c,nfail)

        implicit none

        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class) :: cc
        real(wp) :: rho

        write(*,"(a)") "--- fresh_snow_density: clamps at both ends, both schemes ---"

        ! === constant scheme ==============================================
        cc = c
        cc%fresh_snow_density_scheme = CHION_FRESH_SNOW_DENSITY_CONSTANT

        cc%rho_s = 315.0_wp
        rho = fresh_snow_density(cc,263.15_wp,5.0_wp)
        call check_val("constant, in range -> rho_s", rho, 315.0_wp, nfail)

        ! The constant scheme clamps too: rho_s is NOT passed through blind.
        cc%rho_s = 1.0_wp
        rho = fresh_snow_density(cc,263.15_wp,5.0_wp)
        call check_val("constant, rho_s below 50 -> clamped to 50", &
                       rho, FRESH_SNOW_DENSITY_MIN, nfail)

        cc%rho_s = 5000.0_wp
        rho = fresh_snow_density(cc,263.15_wp,5.0_wp)
        call check_val("constant, rho_s above rho_i -> clamped to rho_i", &
                       rho, cc%rho_i, nfail)

        ! Constant scheme ignores T and wind entirely.
        cc%rho_s = 315.0_wp
        call check_val("constant, ignores air temperature", &
                       fresh_snow_density(cc,150.0_wp,0.0_wp), 315.0_wp, nfail)
        call check_val("constant, ignores wind speed", &
                       fresh_snow_density(cc,263.15_wp,80.0_wp), 315.0_wp, nfail)

        ! === parameterized scheme =========================================
        cc%fresh_snow_density_scheme = CHION_FRESH_SNOW_DENSITY_PARAMETERIZED

        ! a + b*(T-T0) + c*sqrt(V) = 109 + 6*0 + 26*2 = 161  at T = T0, V = 4
        rho = fresh_snow_density(cc,cc%T0,4.0_wp)
        call check_val("parameterized, in range", rho, 161.0_wp, nfail)

        ! Very cold air drives the expression far negative -> lower clamp.
        ! 109 + 6*(-100) + 0 = -491
        rho = fresh_snow_density(cc,cc%T0-100.0_wp,0.0_wp)
        call check_val("parameterized, very cold -> clamped to 50", &
                       rho, FRESH_SNOW_DENSITY_MIN, nfail)

        ! Very warm + very windy drives it above rho_i -> upper clamp.
        ! 109 + 6*200 + 26*sqrt(2500) = 109 + 1200 + 1300 = 2609
        rho = fresh_snow_density(cc,cc%T0+200.0_wp,2500.0_wp)
        call check_val("parameterized, warm and windy -> clamped to rho_i", &
                       rho, cc%rho_i, nfail)

        ! Negative wind is floored at zero before the sqrt (no NaN).
        call check_val("parameterized, negative wind treated as zero", &
                       fresh_snow_density(cc,cc%T0,-9.0_wp), 109.0_wp, nfail)

        ! Monotone in wind speed, over the unclamped range.
        call check("parameterized, increasing in wind speed", &
                   fresh_snow_density(cc,cc%T0,9.0_wp) .gt. &
                   fresh_snow_density(cc,cc%T0,1.0_wp), nfail)

        write(*,*)

        return

    end subroutine test_fresh_snow_density

    ! ======================================================================
    ! (ii) albedo
    ! ======================================================================

    subroutine test_albedo(c,nfail)

        implicit none

        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class) :: cc
        real(wp) :: mass(Ntot), mass_w(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: alb, alb_prev, alb_a, alb_b, alb_exp, span, wet_r
        real(wp) :: age   ! snow age; untouched outside the aging scheme
        integer  :: k, i
        logical  :: ok
        real(wp), parameter :: WET_RS(6) = [0.0_wp, 0.1_wp, 0.5_wp, 0.99_wp, 1.0_wp, 2.0_wp]

        write(*,"(a)") "--- albedo: bounds, monotonicity, saturation, memorylessness ---"

        age = 0.0_wp

        cc = c
        span = cc%alpha_dry - cc%alpha_wet

        ! A plain dry snow column.
        mass        = 0.0_wp
        mass_w      = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        mass(1)        = 200.0_wp
        density(1)     = 350.0_wp
        temperature(1) = 260.0_wp

        ! === bare surface -> alpha_ice ====================================
        cc%albedo_scheme = CHION_ALBEDO_DYNAMIC
        alb = 0.5_wp
        call albedo_update(mass,mass_w,density,temperature,0,cc,1.0_wp,alb)
        call check_val("dynamic, n=0 -> alpha_ice", alb, cc%alpha_ice, nfail)

        mass(1) = 0.0_wp
        alb = 0.5_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
        call check_val("dynamic, empty surface layer -> alpha_ice", alb, cc%alpha_ice, nfail)
        mass(1) = 200.0_wp

        ! === bounds under extreme inputs ==================================
        ! Sweep absurd starting albedos, temperatures and wetnesses; the result
        ! must always land inside [alpha_wet, alpha_dry].
        ok = .TRUE.
        do k = 1, 7
            temperature(1) = 150.0_wp + real(k-1,wp)*30.0_wp      ! 150 .. 330 K
            mass_w(1)      = real(k-1,wp)*20.0_wp                 ! 0 .. 120 kg m-2

            alb = -5.0_wp
            call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
            if (alb .lt. cc%alpha_wet .or. alb .gt. cc%alpha_dry) ok = .FALSE.

            alb = 5.0_wp
            call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
            if (alb .lt. cc%alpha_wet .or. alb .gt. cc%alpha_dry) ok = .FALSE.
        end do
        call check("dynamic, always within [alpha_wet, alpha_dry] under extremes", ok, nfail)

        ! === aging is monotone non-brightening ============================
        ! Repeated daily steps with no snowfall and no liquid water must never
        ! raise the albedo.
        mass_w(1)      = 0.0_wp
        temperature(1) = 265.0_wp
        alb            = cc%alpha_dry

        ok = .TRUE.
        do k = 1, 30
            alb_prev = alb
            call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
            if (alb .gt. alb_prev) ok = .FALSE.
        end do
        call check("dynamic, aging never brightens over 30 calls", ok, nfail)
        call check_val("dynamic, aging floors at alpha_wet", alb, cc%alpha_wet, nfail)

        ! Aging scales with dt (Chion.jl 6d077c5): two half-day steps decay
        ! as much as one daily step, to round-off, away from the floor.
        alb_a = cc%alpha_dry
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_a)
        alb_b = cc%alpha_dry
        call albedo_update(mass,mass_w,density,temperature,1,cc,0.5_wp,alb_b)
        call albedo_update(mass,mass_w,density,temperature,1,cc,0.5_wp,alb_b)
        call check("dynamic, aging decays below alpha_dry in one day", &
                   alb_a .lt. cc%alpha_dry .and. alb_a .gt. cc%alpha_wet, nfail)
        call check("dynamic, two half-day steps == one daily step", &
                   abs(alb_b - alb_a) .le. 4.0_wp*epsilon(1.0_wp), nfail)

        ! Very cold surface: the aging bracket turns negative, and the min()
        ! must then hold the albedo at its previous value rather than brighten.
        temperature(1) = cc%T0 - 100.0_wp
        alb            = 0.75_wp
        alb_prev       = alb
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
        call check("dynamic, very cold surface cannot brighten via aging", &
                   alb .le. alb_prev, nfail)
        temperature(1) = 265.0_wp

        ! === wetness pulls towards alpha_wet ==============================
        mass_w(1) = 0.0_wp
        alb_a     = 0.80_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_a)

        mass_w(1) = 30.0_wp
        alb_b     = 0.80_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_b)

        call check("dynamic, liquid water darkens the surface", alb_b .lt. alb_a, nfail)
        call check("dynamic, wet result still >= alpha_wet", alb_b .ge. cc%alpha_wet, nfail)

        call check("surface_liquid_water_content > 0 for a wet, porous surface", &
                   surface_liquid_water_content(mass,mass_w,density,1,cc) .gt. 0.0_wp_acc, nfail)

        ! Ice-dense surface has no pore space -> lwc must be exactly zero.
        density(1) = cc%rho_i
        call check_val_acc("surface LWC = 0 when the surface is ice-dense", &
                           surface_liquid_water_content(mass,mass_w,density,1,cc), &
                           0.0_wp_acc, nfail)
        density(1) = 350.0_wp
        mass_w(1)  = 0.0_wp

        ! === wetness relaxation composes across substeps ===================
        ! Port of Chion.jl test/test_albedo_timestep.jl (03bb445): at fixed
        ! wetness r, (1-r)**dt makes one daily step, 24 hourly steps, three
        ! uneven steps and two steps of one day vs one of two days agree.
        ! T0 - 30 K disables the temperature aging (its bracket is negative).
        temperature(1) = cc%T0 - 30.0_wp
        density(1)     = 300.0_wp
        mass(1)        = 100.0_wp
        ok = .TRUE.
        do k = 1, 6
            wet_r = WET_RS(k)
            mass_w(1) = wet_r*cc%max_lwc_albedo &
                        *(mass(1)/density(1) - mass(1)/cc%rho_i)*cc%rho_w

            alb_a = cc%alpha_dry
            call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_a)
            alb_exp = cc%alpha_dry - (cc%alpha_dry - cc%alpha_wet)*min(max(wet_r,0.0_wp),1.0_wp)
            if (abs(alb_a - alb_exp) .gt. 64.0_wp*epsilon(1.0_wp)) ok = .FALSE.

            alb_b = cc%alpha_dry
            do i = 1, 24
                call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp/24.0_wp,alb_b)
            end do
            if (abs(alb_b - alb_a) .gt. 64.0_wp*epsilon(1.0_wp)) ok = .FALSE.

            alb_b = cc%alpha_dry
            call albedo_update(mass,mass_w,density,temperature,1,cc,0.1_wp,alb_b)
            call albedo_update(mass,mass_w,density,temperature,1,cc,0.2_wp,alb_b)
            call albedo_update(mass,mass_w,density,temperature,1,cc,0.7_wp,alb_b)
            if (abs(alb_b - alb_a) .gt. 64.0_wp*epsilon(1.0_wp)) ok = .FALSE.

            alb_a = cc%alpha_dry
            call albedo_update(mass,mass_w,density,temperature,1,cc,2.0_wp,alb_a)
            alb_b = cc%alpha_dry
            call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_b)
            call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_b)
            if (abs(alb_b - alb_a) .gt. 64.0_wp*epsilon(1.0_wp)) ok = .FALSE.

            alb_b = cc%alpha_dry
            call albedo_update(mass,mass_w,density,temperature,1,cc,0.0_wp,alb_b)
            if (alb_b .ne. cc%alpha_dry) ok = .FALSE.
        end do
        call check("dynamic, wetness relaxation is substep-invariant ((1-r)**dt)", ok, nfail)
        temperature(1) = 265.0_wp
        density(1)     = 350.0_wp
        mass(1)        = 200.0_wp
        mass_w(1)      = 0.0_wp

        ! === snowfall brightening =========================================
        ! Saturates at alpha_dry, however large the event.
        alb = cc%alpha_wet
        call albedo_refresh_from_snowfall(alb,age,cc,1.0e6_wp)
        call check_val("snowfall brightening saturates at alpha_dry", alb, cc%alpha_dry, nfail)

        ! A single event cannot brighten by more than alpha_dry - alpha_wet.
        ok = .TRUE.
        do k = 1, 8
            alb_prev = cc%alpha_wet + real(k-1,wp)*span/8.0_wp
            alb      = alb_prev
            call albedo_refresh_from_snowfall(alb,age,cc,1.0e9_wp)
            if (alb - alb_prev .gt. span + 8.0_wp*epsilon(1.0_wp)) ok = .FALSE.
            if (alb .gt. cc%alpha_dry) ok = .FALSE.
        end do
        call check("one snowfall event brightens by at most alpha_dry-alpha_wet", ok, nfail)

        ! e-folding is 3 kg m-2: dm = 3 gives exactly (1-1/e) of the span.
        alb = cc%alpha_wet
        call albedo_refresh_from_snowfall(alb,age,cc,ALBEDO_SNOWFALL_EFOLD_MASS)
        call check_val("snowfall e-folding mass is 3 kg m-2", &
                       alb, cc%alpha_wet + span*(1.0_wp - exp(-1.0_wp)), nfail)

        ! Below TOL_TINY of added mass it is a no-op.
        alb = 0.73_wp
        call albedo_refresh_from_snowfall(alb,age,cc,0.0_wp)
        call check_val("zero snowfall mass -> no-op", alb, 0.73_wp, nfail)

        ! Monotone non-darkening in the added mass.
        alb_a = cc%alpha_wet
        call albedo_refresh_from_snowfall(alb_a,age,cc,1.0_wp)
        alb_b = cc%alpha_wet
        call albedo_refresh_from_snowfall(alb_b,age,cc,10.0_wp)
        call check("snowfall brightening increases with added mass", alb_b .gt. alb_a, nfail)

        ! === constant scheme is memoryless ================================
        cc%albedo_scheme = CHION_ALBEDO_CONSTANT

        temperature(1) = 260.0_wp
        alb_a = 0.05_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_a)
        alb_b = 0.95_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_b)
        call check_val("constant, result independent of previous albedo", alb_a, alb_b, nfail)
        call check_val("constant, cold surface -> alpha_dry", alb_a, cc%alpha_dry, nfail)

        temperature(1) = cc%T0
        alb = 0.05_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
        call check_val("constant, surface at T0 -> alpha_wet", alb, cc%alpha_wet, nfail)

        temperature(1) = cc%T0 + 5.0_wp
        alb = 0.05_wp
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
        call check_val("constant, surface above T0 -> alpha_wet", alb, cc%alpha_wet, nfail)

        ! Repeated calls do not drift: memoryless means idempotent.
        alb_prev = alb
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb)
        call check_val("constant, repeated calls are idempotent", alb, alb_prev, nfail)

        ! Constant-scheme snowfall refresh resets to alpha_dry outright.
        alb = cc%alpha_wet
        call albedo_refresh_from_snowfall(alb,age,cc,0.001_wp)
        call check_val("constant, snowfall refresh -> alpha_dry", alb, cc%alpha_dry, nfail)

        ! === prescribed takes the DYNAMIC path (trap 9) ===================
        cc%albedo_scheme = CHION_ALBEDO_PRESCRIBED

        temperature(1) = 265.0_wp
        alb_a          = cc%alpha_dry
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_a)

        cc%albedo_scheme = CHION_ALBEDO_DYNAMIC
        alb_b            = cc%alpha_dry
        call albedo_update(mass,mass_w,density,temperature,1,cc,1.0_wp,alb_b)

        call check_val("prescribed scheme follows the DYNAMIC path (trap 9)", &
                       alb_a, alb_b, nfail)

        write(*,*)

        return

    end subroutine test_albedo

    subroutine test_albedo_aging(c,nfail)
        ! The snowfall-age scheme (Chion.jl 6d06af6): reset on snowfall,
        ! exponential relaxation towards alpha_wet otherwise, with the melting
        ! timescale at T0 and the cold one below.

        implicit none

        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class) :: cc
        real(wp) :: mass(Ntot), temperature(Ntot)
        real(wp) :: alb, age, alb_a, age_a, alb_b, age_b, expected

        write(*,"(a)") "--- albedo: aging scheme (snowfall reset, exponential relaxation) ---"

        cc = c
        cc%albedo_scheme = CHION_ALBEDO_AGING

        mass           = 0.0_wp
        temperature    = 0.0_wp
        mass(1)        = 200.0_wp
        temperature(1) = 260.0_wp

        ! === bare surface -> alpha_ice, age 0 =============================
        alb = 0.75_wp
        age = 4.0_wp
        call albedo_update_aging(mass,temperature,0,cc,0.0_wp,1.0_wp,alb,age)
        call check_val("aging, n=0 -> alpha_ice", alb, cc%alpha_ice, nfail)
        call check_val("aging, n=0 -> snow age 0", age, 0.0_wp, nfail)

        ! === snowfall: Chion.jl's binary reset under legacy_chion only ====
        ! Otherwise the refresh happens in accumulation (D30) and the update
        ! ages the snow regardless of the snowfall rate.
        alb = 0.72_wp
        age = 5.0_wp
        call albedo_update_aging(mass,temperature,1,cc,1.0e-9_wp,1.0_wp,alb,age)
        if (ALBEDO_AGING_BINARY_REFRESH) then
            call check_val("aging (legacy), any snowfall -> alpha_dry", alb, cc%alpha_dry, nfail)
            call check_val("aging (legacy), any snowfall -> snow age 0", age, 0.0_wp, nfail)
        else
            call check_val("aging, snowfall rate does not reset in the update", alb, &
                           cc%alpha_wet + (0.72_wp - cc%alpha_wet) &
                                         *exp(-1.0_wp/cc%aging_cold_timescale_days), nfail)
            call check_val("aging, snowfall rate does not reset the age", age, 6.0_wp, nfail)
        end if

        ! === relaxation, cold surface: tau = aging_cold_timescale_days ====
        alb = cc%alpha_dry
        age = 0.0_wp
        call albedo_update_aging(mass,temperature,1,cc,0.0_wp,1.0_wp,alb,age)
        expected = cc%alpha_wet + (cc%alpha_dry - cc%alpha_wet) &
                                 *exp(-1.0_wp/cc%aging_cold_timescale_days)
        call check_val("aging, cold relaxation alpha_wet + (a-alpha_wet)*exp(-dt/tau_cold)", &
                       alb, expected, nfail)
        call check_val("aging, snow age advances by dt", age, 1.0_wp, nfail)

        ! === relaxation, melting surface: tau = aging_melting_timescale_days
        temperature(1) = cc%T0
        alb = cc%alpha_dry
        age = 3.0_wp
        call albedo_update_aging(mass,temperature,1,cc,0.0_wp,0.5_wp,alb,age)
        expected = cc%alpha_wet + (cc%alpha_dry - cc%alpha_wet) &
                                 *exp(-0.5_wp/cc%aging_melting_timescale_days)
        call check_val("aging, melting relaxation uses tau_melt at Ts = T0", alb, expected, nfail)
        call check_val("aging, snow age advances by a partial dt", age, 3.5_wp, nfail)

        ! === dt composition: two half-day steps == one daily step =========
        temperature(1) = 260.0_wp
        alb_a = cc%alpha_dry
        age_a = 0.0_wp
        call albedo_update_aging(mass,temperature,1,cc,0.0_wp,1.0_wp,alb_a,age_a)
        alb_b = cc%alpha_dry
        age_b = 0.0_wp
        call albedo_update_aging(mass,temperature,1,cc,0.0_wp,0.5_wp,alb_b,age_b)
        call albedo_update_aging(mass,temperature,1,cc,0.0_wp,0.5_wp,alb_b,age_b)
        call check("aging, two half-day steps == one daily step", &
                   abs(alb_b - alb_a) .le. 4.0_wp*epsilon(1.0_wp), nfail)
        call check_val("aging, two half-day steps age one day", age_b, 1.0_wp, nfail)

        ! === out-of-range previous albedo is clamped first ================
        alb = cc%alpha_ice
        age = -2.0_wp
        call albedo_update_aging(mass,temperature,1,cc,0.0_wp,1.0_wp,alb,age)
        call check_val("aging, previous alpha_ice is clamped up to alpha_wet", &
                       alb, cc%alpha_wet, nfail)
        call check_val("aging, negative snow age is clamped to 0 first", age, 1.0_wp, nfail)

        ! === snowfall refresh in the accumulation slot ====================
        ! Legacy: Chion.jl's reset. Otherwise D30's partial refresh, which a
        ! trace of snow on fully aged snow leaves at alpha_wet.
        alb = cc%alpha_wet
        age = 4.0_wp
        call albedo_refresh_from_snowfall(alb,age,cc,0.001_wp)
        if (ALBEDO_AGING_BINARY_REFRESH) then
            call check_val("aging (legacy), snowfall refresh -> alpha_dry", alb, cc%alpha_dry, nfail)
        else
            call check_val("aging, trace snowfall on alpha_wet stays alpha_wet", alb, cc%alpha_wet, nfail)
            call check_val("aging, trace snowfall scales the age by (1-f)", age, &
                           4.0_wp*(1.0_wp - 0.001_wp/cc%aging_snowfall_ref), nfail)
        end if

        write(*,*)

        return

    end subroutine test_albedo_aging

    subroutine test_albedo_aging_rejuvenate(c,nfail)
        ! D30: partial rejuvenation f = min(1, S/S_ref) of the aging albedo,
        ! E = -ln((a-alpha_wet)/(alpha_dry-alpha_wet)) <- (1-f)*E. Tested on
        ! the routine itself, which legacy_chion builds do not call.

        implicit none

        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class) :: cc
        real(wp) :: alb, age, alb_a, age_a, alb_prev, span, e0, e2, s1, s2
        integer  :: k
        logical  :: ok

        write(*,"(a)") "--- albedo: aging rejuvenation by snowfall (D30) ---"

        cc = c
        cc%albedo_scheme = CHION_ALBEDO_AGING
        span = cc%alpha_dry - cc%alpha_wet

        ! S = 0 leaves albedo and age unchanged.
        alb = 0.75_wp
        age = 6.0_wp
        call albedo_aging_rejuvenate(alb,age,cc,0.0_wp)
        call check_val("rejuvenate, S = 0 leaves albedo unchanged", alb, 0.75_wp, nfail)
        call check_val("rejuvenate, S = 0 leaves snow age unchanged", age, 6.0_wp, nfail)

        ! S >= S_ref: Chion.jl's reset, from anywhere, including alpha_wet.
        alb = cc%alpha_wet
        age = 30.0_wp
        call albedo_aging_rejuvenate(alb,age,cc,cc%aging_snowfall_ref)
        call check_val("rejuvenate, S = S_ref -> alpha_dry", alb, cc%alpha_dry, nfail)
        call check_val("rejuvenate, S = S_ref -> snow age 0", age, 0.0_wp, nfail)

        alb = 0.74_wp
        age = 3.0_wp
        call albedo_aging_rejuvenate(alb,age,cc,5.0_wp*cc%aging_snowfall_ref)
        call check_val("rejuvenate, S > S_ref -> alpha_dry", alb, cc%alpha_dry, nfail)
        call check_val("rejuvenate, S > S_ref -> snow age 0", age, 0.0_wp, nfail)

        ! Half the reference snowfall halves the aging progress E.
        alb = 0.75_wp
        age = 8.0_wp
        e0  = -log((alb - cc%alpha_wet)/span)
        call albedo_aging_rejuvenate(alb,age,cc,0.5_wp*cc%aging_snowfall_ref)
        call check_val("rejuvenate, S = S_ref/2 halves E", &
                       -log((alb - cc%alpha_wet)/span), 0.5_wp*e0, nfail)
        call check_val("rejuvenate, S = S_ref/2 halves the snow age", age, 4.0_wp, nfail)

        ! Successive refreshes compose: E*(1-f1)*(1-f2).
        s1  = 0.3_wp*cc%aging_snowfall_ref
        s2  = 0.6_wp*cc%aging_snowfall_ref
        alb = 0.72_wp
        age = 10.0_wp
        e0  = -log((alb - cc%alpha_wet)/span)
        call albedo_aging_rejuvenate(alb,age,cc,s1)
        call albedo_aging_rejuvenate(alb,age,cc,s2)
        e2  = -log((alb - cc%alpha_wet)/span)
        call check_val("rejuvenate, two refreshes compose as E*(1-f1)*(1-f2)", &
                       e2, e0*(1.0_wp - 0.3_wp)*(1.0_wp - 0.6_wp), nfail)
        call check_val("rejuvenate, two refreshes compose on the snow age", &
                       age, 10.0_wp*(1.0_wp - 0.3_wp)*(1.0_wp - 0.6_wp), nfail)

        ! Monotone in S: more snow never darkens, from below alpha_dry.
        ok = .TRUE.
        alb_prev = 0.0_wp
        do k = 0, 12
            alb = 0.73_wp
            age = 5.0_wp
            call albedo_aging_rejuvenate(alb,age,cc,real(k,wp)*0.1_wp*cc%aging_snowfall_ref)
            if (alb .lt. alb_prev) ok = .FALSE.
            if (alb .lt. 0.73_wp .or. alb .gt. cc%alpha_dry) ok = .FALSE.
            alb_prev = alb
        end do
        call check("rejuvenate, albedo monotone non-decreasing in S, within bounds", ok, nfail)

        ! alpha_wet exactly (E infinite): only a full refresh restores it, and
        ! alpha_dry = alpha_wet is well defined (no 0/0).
        alb = cc%alpha_wet
        age = 2.0_wp
        call albedo_aging_rejuvenate(alb,age,cc,0.9_wp*cc%aging_snowfall_ref)
        call check_val("rejuvenate, at alpha_wet a partial refresh stays alpha_wet", &
                       alb, cc%alpha_wet, nfail)

        cc%alpha_wet = cc%alpha_dry
        alb_a = 0.5_wp
        age_a = 2.0_wp
        call albedo_aging_rejuvenate(alb_a,age_a,cc,0.5_wp*cc%aging_snowfall_ref)
        call check_val("rejuvenate, alpha_dry = alpha_wet -> that albedo", &
                       alb_a, cc%alpha_dry, nfail)

        write(*,*)

        return

    end subroutine test_albedo_aging_rejuvenate

    ! ======================================================================
    ! (iii) densification
    ! ======================================================================

    subroutine test_densification(c,nfail)

        implicit none

        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class) :: cc
        real(wp) :: mass(Ntot), mass_w(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: rho_before(Ntot), mass_w_prev(Ntot)
        real(wp) :: dt_seconds, acc_rate
        integer  :: k, istep, ischeme
        logical  :: ok_monotone, ok_capped, ok_moved

        write(*,"(a)") "--- densification: monotone, capped at rho_i, all regimes, both schemes ---"

        dt_seconds = 86400.0_wp
        acc_rate   = 1.0e-5_wp      ! [kg m-2 s-1] ~0.86 kg m-2 d-1

        ! Layer mass matters here. The mid and high branches scale with the
        ! CUBE of the overburden pressure, so a thin test column produces
        ! tendencies far below sp resolution and looks (wrongly) inert. Real
        ! firn at 600-850 kg m-3 sits under tens of metres of overburden, so
        ! the test column uses 2e4 kg m-2 per layer (~25 m of firn).
        do ischeme = 1, 2

            cc = c
            if (ischeme .eq. 1) then
                cc%low_density_densification = CHION_DENSIFY_BESSI
            else
                cc%low_density_densification = CHION_DENSIFY_HTESSEL
            end if

            ! One layer in each of the three regimes, plus an at-cap layer and
            ! a massless layer.
            mass        = 20000.0_wp
            mass_w      = 0.0_wp
            temperature = 263.15_wp

            density(1) = 300.0_wp     ! low  regime, rho < 550
            density(2) = 600.0_wp     ! mid  regime, 550 <= rho < 800
            density(3) = 850.0_wp     ! high regime, rho >= 800
            density(4) = cc%rho_i     ! already at the cap
            density(5) = 400.0_wp     ! low regime but massless
            mass(5)    = 0.0_wp

            ok_monotone = .TRUE.
            ok_capped   = .TRUE.

            do istep = 1, 200
                rho_before = density
                call densify_column(mass,density,temperature,Ntot,cc,acc_rate,dt_seconds)
                do k = 1, Ntot
                    if (density(k) .lt. rho_before(k)) ok_monotone = .FALSE.
                    if (density(k) .gt. cc%rho_i)      ok_capped   = .FALSE.
                end do
            end do

            if (ischeme .eq. 1) then
                call check("BESSI: density is monotone non-decreasing", ok_monotone, nfail)
                call check("BESSI: density never exceeds rho_i", ok_capped, nfail)
                call check("BESSI: low-density layer densified", &
                           density(1) .gt. 300.0_wp, nfail)
            else
                call check("HTESSEL: density is monotone non-decreasing", ok_monotone, nfail)
                call check("HTESSEL: density never exceeds rho_i", ok_capped, nfail)
                call check("HTESSEL: low-density layer densified", &
                           density(1) .gt. 300.0_wp, nfail)
            end if

            call check("mid regime densified (550 <= rho < 800)", &
                       density(2) .gt. 600.0_wp, nfail)
            call check("high regime densified (rho >= 800)", &
                       density(3) .gt. 850.0_wp, nfail)
            call check_val("layer already at rho_i is snapped and held", &
                           density(4), cc%rho_i, nfail)
            call check_val("massless layer is left untouched", &
                           density(5), 400.0_wp, nfail)

        end do

        ! === the scheme flag selects ONLY the low-density branch ===========
        mass        = 20000.0_wp
        temperature = 263.15_wp
        density     = 600.0_wp

        cc = c
        cc%low_density_densification = CHION_DENSIFY_BESSI
        call densify_column(mass,density,temperature,2,cc,acc_rate,dt_seconds)
        rho_before = density

        density = 600.0_wp
        cc%low_density_densification = CHION_DENSIFY_HTESSEL
        call densify_column(mass,density,temperature,2,cc,acc_rate,dt_seconds)

        call check_val("mid-regime tendency is scheme-independent", &
                       density(1), rho_before(1), nfail)

        ! === the two low-density schemes actually differ ===================
        density = 300.0_wp
        cc%low_density_densification = CHION_DENSIFY_BESSI
        call densify_column(mass,density,temperature,2,cc,acc_rate,dt_seconds)
        rho_before = density

        density = 300.0_wp
        cc%low_density_densification = CHION_DENSIFY_HTESSEL
        call densify_column(mass,density,temperature,2,cc,acc_rate,dt_seconds)

        call check("low-density branch differs between the two schemes", &
                   abs(density(1) - rho_before(1)) .gt. 0.0_wp, nfail)

        ! === BESSI rate vanishes for non-positive accumulation ============
        density = 300.0_wp
        cc%low_density_densification = CHION_DENSIFY_BESSI
        call densify_column(mass,density,temperature,2,cc,-1.0_wp,dt_seconds)
        call check_val("BESSI: negative accumulation gives zero tendency", &
                       density(1), 300.0_wp, nfail)

        ! === overburden increases with depth ==============================
        ! Two identical mid-regime layers: the deeper one carries more overburden
        ! and must therefore densify more.
        mass        = 20000.0_wp
        density     = 700.0_wp
        temperature = 263.15_wp
        call densify_column(mass,density,temperature,3,cc,acc_rate,dt_seconds)
        call check("deeper layers densify more (overburden increases downward)", &
                   density(3) .gt. density(2) .and. density(2) .gt. density(1), nfail)

        ! === bubble pressure ==============================================
        call check_val_acc("bubble pressure is zero below rho_e", &
                           bubble_pressure_mpa(700.0_wp,c%rho_i), 0.0_wp_acc, nfail)
        call check_val_acc("bubble pressure is zero exactly at rho_e", &
                           bubble_pressure_mpa(real(DENSIFY_RHO_E,wp),c%rho_i), &
                           0.0_wp_acc, nfail)
        call check("bubble pressure is positive above rho_e", &
                   bubble_pressure_mpa(880.0_wp,c%rho_i) .gt. 0.0_wp_acc, nfail)

        ! === HTESSEL liquid-water compaction ==============================
        mass        = 100.0_wp
        mass_w      = 0.0_wp
        mass_w_prev = 0.0_wp
        temperature = 265.0_wp

        density(1) = 300.0_wp    ! low  -> eligible
        density(2) = 600.0_wp    ! high -> not eligible
        density(3) = 300.0_wp    ! low, but loses water -> no change
        density(4) = 900.0_wp
        density(5) = 300.0_wp    ! low, but no water gain

        mass_w_prev(1) = 0.0_wp  ; mass_w(1) = 10.0_wp    ! gains 10
        mass_w_prev(2) = 0.0_wp  ; mass_w(2) = 10.0_wp    ! gains 10 but too dense
        mass_w_prev(3) = 10.0_wp ; mass_w(3) = 2.0_wp     ! loses water
        mass_w_prev(4) = 0.0_wp  ; mass_w(4) = 10.0_wp
        mass_w_prev(5) = 5.0_wp  ; mass_w(5) = 5.0_wp     ! no change

        rho_before = density
        call apply_htessel_liquid_water_compaction(mass,mass_w,density,Ntot,mass_w_prev,c)

        ! rho*(1 + dmw/m) = 300*(1 + 10/100) = 330
        call check_val("HTESSEL compaction: rho -> rho*(1 + dmw/m)", &
                       density(1), 330.0_wp, nfail)
        call check_val("HTESSEL compaction: skips rho >= 550", &
                       density(2), rho_before(2), nfail)
        call check_val("HTESSEL compaction: no-op when water is lost", &
                       density(3), rho_before(3), nfail)
        call check_val("HTESSEL compaction: no-op when water is unchanged", &
                       density(5), rho_before(5), nfail)

        ! Cap at rho_i even for an absurd water gain.
        density(1)     = 500.0_wp
        mass_w_prev(1) = 0.0_wp
        mass_w(1)      = 1.0e6_wp
        call apply_htessel_liquid_water_compaction(mass,mass_w,density,1,mass_w_prev,c)
        call check_val("HTESSEL compaction: capped at rho_i", density(1), c%rho_i, nfail)

        ! Never decreases density.
        ok_moved = .TRUE.
        density  = 300.0_wp
        mass_w   = 1.0_wp
        mass_w_prev = 0.0_wp
        rho_before  = density
        call apply_htessel_liquid_water_compaction(mass,mass_w,density,Ntot,mass_w_prev,c)
        do k = 1, Ntot
            if (density(k) .lt. rho_before(k)) ok_moved = .FALSE.
        end do
        call check("HTESSEL compaction is monotone non-decreasing", ok_moved, nfail)

        write(*,*)

        return

    end subroutine test_densification

    ! ======================================================================
    ! (iv) diurnal shortwave
    ! ======================================================================

    subroutine test_diurnal(nfail)

        use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer,  parameter :: nsub_max = 24
        real(wp), parameter :: PI_WP = 3.14159265358979_wp

        real(wp)     :: lat, lon, qbar, h_a, h_b, q, tbar, amp, t_int, dec, h0, nan_wp
        real(wp_acc) :: total, weight, w_total, I_day, A, B
        integer      :: i, n, ilat, ilon, nsub
        logical      :: ok_nonneg, ok_tiling

        write(*,"(a)") "--- diurnal: tiling recovers the daily mean, polar limits, substeps ---"

        qbar = 250.0_wp

        ! === geometry sanity ==============================================
        call check_val("declination is zero at solar longitude 0", &
                       solar_declination_deg(0.0_wp), 0.0_wp, nfail)
        call check_val("declination = +obliquity at solar longitude 90", &
                       solar_declination_deg(90.0_wp), DIURNAL_OBLIQUITY_DEG, nfail)
        call check_val("declination = -obliquity at solar longitude 270", &
                       solar_declination_deg(270.0_wp), -DIURNAL_OBLIQUITY_DEG, nfail)

        call check_val("equator at equinox: sunset hour angle = pi/2", &
                       sunset_hour_angle(0.0_wp,0.0_wp), PI_WP/2.0_wp, nfail)

        ! === full [-pi,pi] tiling recovers the daily mean =================
        ! sum_i (w_i/2pi) * Qbar_i  ==  Qbar, for any tiling.
        ok_tiling = .TRUE.

        do ilat = 1, 5
            lat = -80.0_wp + real(ilat-1,wp)*40.0_wp     ! -80 .. 80
            do ilon = 1, 4
                lon = real(ilon-1,wp)*90.0_wp            ! 0, 90, 180, 270

                ! Skip polar night/day-degenerate configurations, which are
                ! covered explicitly below.
                dec = solar_declination_deg(lon)
                h0  = sunset_hour_angle(lat,dec)
                if (h0 .le. 0.0_wp) cycle

                do n = 2, 12

                    total   = 0.0_wp_acc
                    w_total = 0.0_wp_acc

                    do i = 1, n
                        call diurnal_substep_bounds(i,n,h_a,h_b)
                        q = diurnal_shortwave_interval_average(qbar,lat,lon,h_a,h_b)
                        weight  = real(h_b,wp_acc) - real(h_a,wp_acc)
                        total   = total   + weight*real(q,wp_acc)
                        w_total = w_total + weight
                    end do

                    ! Weights must tile the full day exactly.
                    if (abs(w_total - 2.0_wp_acc*real(PI_WP,wp_acc)) .gt. 1.0e-6_wp_acc) &
                        ok_tiling = .FALSE.

                    ! Energy-conserving reconstruction.
                    if (abs(total/w_total - real(qbar,wp_acc))/real(qbar,wp_acc) &
                        .gt. 1.0e-6_wp_acc) ok_tiling = .FALSE.

                end do
            end do
        end do

        call check("tiling [-pi,pi] recovers the daily mean to 1e-6 relative &
                   &(5 latitudes x 4 solar longitudes x n=2..12)", ok_tiling, nfail)

        ! === nocturnal intervals dilute, never go negative ================
        lat = 65.0_wp
        lon = 90.0_wp

        ok_nonneg = .TRUE.
        do i = 1, 24
            call diurnal_substep_bounds(i,24,h_a,h_b)
            q = diurnal_shortwave_interval_average(qbar,lat,lon,h_a,h_b)
            if (q .lt. 0.0_wp) ok_nonneg = .FALSE.
        end do
        call check("interval average is never negative", ok_nonneg, nfail)

        ! A fully nocturnal interval (just past sunset) returns exactly zero.
        dec = solar_declination_deg(lon)
        h0  = sunset_hour_angle(lat,dec)
        call check_val("fully nocturnal interval returns zero", &
                       diurnal_shortwave_interval_average(qbar,lat,lon, &
                                                          h0+0.1_wp,h0+0.2_wp), &
                       0.0_wp, nfail)

        ! Zero-width and reversed intervals return zero.
        call check_val("zero-width interval returns zero", &
                       diurnal_shortwave_interval_average(qbar,lat,lon,0.5_wp,0.5_wp), &
                       0.0_wp, nfail)
        call check_val("reversed interval returns zero", &
                       diurnal_shortwave_interval_average(qbar,lat,lon,0.5_wp,0.1_wp), &
                       0.0_wp, nfail)
        call check_val("non-positive daily mean returns zero", &
                       diurnal_shortwave_interval_average(0.0_wp,lat,lon,-1.0_wp,1.0_wp), &
                       0.0_wp, nfail)

        ! === polar night ==================================================
        ! Northern winter (lon = 270 -> declination = -23.44) at 80 N.
        lat = 80.0_wp
        lon = 270.0_wp
        call diurnal_daylight_integral(lat,lon,dec,h0,I_day,A,B)
        call check_val("polar night: sunset hour angle = 0", h0, 0.0_wp, nfail)
        call check_val("polar night: interval average = 0 over the whole day", &
                       diurnal_shortwave_interval_average(qbar,lat,lon,-PI_WP,PI_WP), &
                       0.0_wp, nfail)
        call check_val("polar night: peak flux = 0", &
                       diurnal_shortwave_peak_flux(qbar,lat,lon), 0.0_wp, nfail)

        ! === polar day ====================================================
        ! Northern summer (lon = 90 -> declination = +23.44) at 80 N.
        lat = 80.0_wp
        lon = 90.0_wp
        call diurnal_daylight_integral(lat,lon,dec,h0,I_day,A,B)
        call check_val("polar day: sunset hour angle = pi", h0, PI_WP, nfail)

        q = diurnal_shortwave_interval_average(qbar,lat,lon,-PI_WP,PI_WP)
        call check_val("polar day: whole-day average returns the daily mean", q, qbar, nfail)

        call check("polar day: peak flux exceeds the daily mean", &
                   diurnal_shortwave_peak_flux(qbar,lat,lon) .gt. qbar, nfail)
        call check("polar day: sun is up in every interval", &
                   diurnal_shortwave_interval_average(qbar,lat,lon,PI_WP-0.2_wp,PI_WP) &
                   .gt. 0.0_wp, nfail)

        ! === temperature interval average =================================
        tbar = 270.0_wp
        amp  = 5.0_wp

        call check_val("temperature: whole-day average returns the daily mean", &
                       diurnal_temperature_interval_average(tbar,amp,-PI_WP,PI_WP), &
                       tbar, nfail)
        call check_val("temperature: zero amplitude returns the daily mean", &
                       diurnal_temperature_interval_average(tbar,0.0_wp,-1.0_wp,1.0_wp), &
                       tbar, nfail)
        call check_val("temperature: non-positive width returns the daily mean", &
                       diurnal_temperature_interval_average(tbar,amp,1.0_wp,1.0_wp), &
                       tbar, nfail)

        ! Warmest around solar noon (h = 0), coldest around midnight.
        call check("temperature: warmest interval brackets solar noon", &
                   diurnal_temperature_interval_average(tbar,amp,-0.5_wp,0.5_wp) .gt. &
                   diurnal_temperature_interval_average(tbar,amp,PI_WP-0.5_wp,PI_WP), nfail)

        ! The tiling preserves the daily mean temperature too.
        total   = 0.0_wp_acc
        w_total = 0.0_wp_acc
        do i = 1, 6
            call diurnal_substep_bounds(i,6,h_a,h_b)
            t_int   = diurnal_temperature_interval_average(tbar,amp,h_a,h_b)
            weight  = real(h_b,wp_acc) - real(h_a,wp_acc)
            total   = total   + weight*real(t_int,wp_acc)
            w_total = w_total + weight
        end do
        call check("temperature tiling preserves the daily mean", &
                   abs(total/w_total - real(tbar,wp_acc)) .le. 1.0e-4_wp_acc, nfail)

        ! === elevation-dependent amplitude (Chion.jl d0146e1) ==============
        ! A = clamp(A0 + gamma*max(z - z_ref, 0), 0, A_max), gamma in K/km.
        call check("amplitude: neutral defaults return A0 exactly", &
                   diurnal_temperature_amplitude(5.0_wp,0.0_wp,0.0_wp,1.0e30_wp,3000.0_wp) &
                   .eq. 5.0_wp, nfail)
        call check_val("amplitude: gradient above z_ref", &
                       diurnal_temperature_amplitude(1.0_wp,2.0_wp,1000.0_wp,10.0_wp,2500.0_wp), &
                       4.0_wp, nfail)
        call check_val("amplitude: no excess below z_ref", &
                       diurnal_temperature_amplitude(1.0_wp,2.0_wp,1000.0_wp,10.0_wp,500.0_wp), &
                       1.0_wp, nfail)
        call check_val("amplitude: clamped at A_max", &
                       diurnal_temperature_amplitude(1.0_wp,2.0_wp,0.0_wp,3.0_wp,4000.0_wp), &
                       3.0_wp, nfail)
        call check_val("amplitude: negative gradient clamped at zero", &
                       diurnal_temperature_amplitude(1.0_wp,-2.0_wp,0.0_wp,3.0_wp,4000.0_wp), &
                       0.0_wp, nfail)
        nan_wp = ieee_value(nan_wp,ieee_quiet_nan)
        amp = diurnal_temperature_amplitude(1.0_wp,2.0_wp,0.0_wp,3.0_wp,nan_wp)
        call check("amplitude: missing (NaN) surface height gives A0, not NaN", &
                   ieee_is_finite(amp) .and. amp .eq. 1.0_wp, nfail)

        ! === substep_count returns only 1 or max_substeps ==================
        lat = 65.0_wp
        lon = 90.0_wp

        nsub = diurnal_substep_count(1.0_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("substepping activates in polar summer", nsub, 3, nfail)

        nsub = diurnal_substep_count(1.0_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,1)
        call check_int("max_substeps = 1 -> 1", nsub, 1, nfail)

        nsub = diurnal_substep_count(0.5_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("dt_days below 0.75 -> 1", nsub, 1, nfail)

        nsub = diurnal_substep_count(2.0_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("dt_days above 1.25 -> 1", nsub, 1, nfail)

        nsub = diurnal_substep_count(0.75_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("dt_days exactly 0.75 is inside the window", nsub, 3, nfail)

        nsub = diurnal_substep_count(1.25_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("dt_days exactly 1.25 is inside the window", nsub, 3, nfail)

        nsub = diurnal_substep_count(1.0_wp,0.0_wp,270.0_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("zero daily mean shortwave -> 1", nsub, 1, nfail)

        nsub = diurnal_substep_count(1.0_wp,qbar,265.15_wp,265.15_wp,lat,lon,0.0_wp,3)
        call check_int("air temperature EQUAL to the minimum -> 1 (strict test)", &
                       nsub, 1, nfail)

        nsub = diurnal_substep_count(1.0_wp,qbar,270.0_wp,265.15_wp,lat,lon,1.0e6_wp,3)
        call check_int("peak-minus-mean below threshold -> 1", nsub, 1, nfail)

        nsub = diurnal_substep_count(1.0_wp,qbar,270.0_wp,265.15_wp,80.0_wp,270.0_wp,0.0_wp,3)
        call check_int("polar night -> 1", nsub, 1, nfail)

        ! Only ever 1 or max_substeps, never an intermediate value.
        ok_tiling = .TRUE.
        do i = 1, nsub_max
            nsub = diurnal_substep_count(1.0_wp,qbar,270.0_wp,265.15_wp,lat,lon,0.0_wp,i)
            if (nsub .ne. 1 .and. nsub .ne. i) ok_tiling = .FALSE.
        end do
        call check("substep_count returns only 1 or max_substeps (1..24)", ok_tiling, nfail)

        write(*,*)

        return

    end subroutine test_diurnal

    ! ======================================================================
    ! Check helpers (same style as tests/test_column_utils.f90)
    ! ======================================================================

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

        tol = max(abs(expected)*32.0_wp*epsilon(1.0_wp),1.0e-6_wp)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_val

    subroutine check_val_acc(label,value,expected,nfail)

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

    end subroutine check_val_acc

    subroutine check_int(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        integer,          intent(IN)    :: value
        integer,          intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        if (value .eq. expected) then
            write(*,"(a,a,a,i0)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,i0,a,i0)") "  FAIL : ", trim(label), &
                                       " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_int

end program test_wp7
