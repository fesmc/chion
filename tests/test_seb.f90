program test_seb
    ! Acceptance test for snow_seb_semix: the CLIMBER-X SEMIX aerodynamic
    ! turbulent exchange (docs/semix_port_scope.md, rung 2).
    !
    ! PRECISION: the default build is wp = sp (single). Hand-computed values
    ! are therefore checked to ~1e-5 relative, and the two num/denom identities
    ! at the end to 1e-4: they subtract two quantities near 1100 and 2100 W m-2
    ! to leave 6 and 40 W m-2, so single-precision cancellation dominates.
    !
    ! Covers:
    !   * semix_resistance against hand-computed values, in all three
    !     stability regimes plus the dead-calm guard
    !   * semix_air_density against the ideal gas law
    !   * both saturation-humidity variants (c%semix_qsat) against
    !     hand-computed values, and their derivatives against finite
    !     differences of the humidity itself
    !   * the exchange coefficients f_sh / f_lh, the no-humidity case and the
    !     dew-inhibition switch l_dew
    !   * the sign convention (positive INTO the surface, opposite to SEMIX's
    !     own diagnostics)
    !   * THE COUPLING IDENTITY: ebal's num_lh/denom_lh decomposition, which is
    !     what the energy solve will fold into q_const/q_lin, reproduces the
    !     exact latent flux at the linearization point
    !
    ! No call site uses this module yet -- the flux sites are wired in the
    ! following commit. This test pins the physics on its own first.

    use chion_defs,   only : wp, chion_const_class, chion_const_init, &
                             DEF_SEA_LEVEL_AIR_PRESSURE, &
                             SEMIX_QSAT_SEMIX, SEMIX_QSAT_BESSI, &
                             CHION_SEB_BESSI, CHION_SEB_SEMIX, &
                             chion_seb_scheme_flag, chion_semix_qsat_flag
    use snow_seb_semix

    implicit none

    type(chion_const_class)    :: c
    type(semix_exchange_class) :: x, x_warm

    real(wp), parameter :: P0 = DEF_SEA_LEVEL_AIR_PRESSURE
    real(wp), parameter :: T_AIR  = 263.15_wp     ! [K] -10 C
    real(wp), parameter :: T_COLD = 258.15_wp     ! [K] surface 5 K below air
    real(wp), parameter :: T_WARM = 268.15_wp     ! [K] surface 5 K above air
    real(wp), parameter :: WIND   = 5.0_wp        ! [m s-1]
    real(wp), parameter :: H_DEEP = 1.0_wp        ! [m] snow depth

    real(wp) :: r_stable, r_unstable, r_neutral, r_calm
    real(wp) :: rhoa, q_a, dq_fd, num_lh, denom_lh
    integer  :: nfail

    nfail = 0

    call chion_const_init(c)

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion acceptance test: snow_seb_semix (SEMIX SEB, rung 2)"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! === Scheme flags ====================================================
    write(*,"(a)") "--- scheme flags ---"

    call check("seb_scheme default is bessi", c%seb_scheme .eq. CHION_SEB_BESSI, nfail)
    call check("semix_qsat default is semix", c%semix_qsat .eq. SEMIX_QSAT_SEMIX, nfail)
    call check("flag('bessi') = CHION_SEB_BESSI", &
               chion_seb_scheme_flag("bessi") .eq. CHION_SEB_BESSI, nfail)
    call check("flag('semix') = CHION_SEB_SEMIX", &
               chion_seb_scheme_flag("semix") .eq. CHION_SEB_SEMIX, nfail)
    call check("qsat flag('chion') = SEMIX_QSAT_BESSI", &
               chion_semix_qsat_flag("chion") .eq. SEMIX_QSAT_BESSI, nfail)
    call check("qsat flag('climberx') = SEMIX_QSAT_SEMIX", &
               chion_semix_qsat_flag("climberx") .eq. SEMIX_QSAT_SEMIX, nfail)

    ! === Aerodynamic resistance ==========================================
    ! Hand-computed for h_snow = 1 m, wind = 5 m s-1, z0m_ice = 0.002,
    ! z0m_snow = 0.0024, zm_to_zh = exp(-2), z_sfl = 100:
    !   fsnow      = 1/(1 + 10*0.002)       = 0.9803921569
    !   rough_m    = 0.9803921569*0.0024 + 0.0196078431*0.002
    !              = 0.0023921569
    !   rough_h    = rough_m*exp(-2)        = 3.2373e-4
    !   Ch_neutral = (0.4/ln(100/rough_m))*(0.4/ln(100/rough_h))
    !              = 1.1895327409e-3
    !   Ri         = 9.81*100*(1 - T_s/263.15)/25
    write(*,*)
    write(*,"(a)") "--- semix_resistance ---"

    ! Stable (surface colder than air): Ri = +0.745582, Ch = Ch_neutral, so
    ! r_a = 1/(1.1895327409e-3 * 5) = 168.1332452
    r_stable = semix_resistance(H_DEEP,T_AIR,T_COLD,WIND,c)
    call check_close("r_a stable  (T_s = T_air - 5)", r_stable, 168.1332452_wp, 1.0e-5_wp, nfail)

    ! Unstable (surface warmer than air): Ri = -0.745582, so the exchange is
    ! enhanced by (1 - 2*Ri) = 2.491164 and r_a = 168.1332452/2.491164
    r_unstable = semix_resistance(H_DEEP,T_AIR,T_WARM,WIND,c)
    call check_close("r_a unstable (T_s = T_air + 5)", r_unstable, 67.4918213_wp, 1.0e-5_wp, nfail)

    ! The stable branch is deliberately NOT damped in SEMIX -- it just uses the
    ! neutral coefficient -- so forcing neutral must reproduce it exactly.
    c%l_neutral = .TRUE.
    r_neutral = semix_resistance(H_DEEP,T_AIR,T_COLD,WIND,c)
    call check("l_neutral leaves the stable branch unchanged", &
               r_neutral .eq. r_stable, nfail)
    r_neutral = semix_resistance(H_DEEP,T_AIR,T_WARM,WIND,c)
    call check("l_neutral suppresses the unstable enhancement", &
               r_neutral .gt. r_unstable, nfail)
    c%l_neutral = .FALSE.

    ! Dead calm: Ch*wind = 0, so SEMIX returns its 1e20 sentinel rather than
    ! dividing by zero.
    r_calm = semix_resistance(H_DEEP,T_AIR,T_COLD,0.0_wp,c)
    call check_close("r_a at zero wind = 1e20", r_calm, 1.0e20_wp, 1.0e-12_wp, nfail)

    ! Bare ground uses the ice roughness alone; deep snow the snow roughness.
    ! z0m_snow > z0m_ice, so deep snow is the ROUGHER surface and exchanges
    ! more freely.
    call check("deep snow has lower r_a than bare ice", &
               semix_resistance(H_DEEP,T_AIR,T_COLD,WIND,c) .lt. &
               semix_resistance(0.0_wp,T_AIR,T_COLD,WIND,c), nfail)

    ! === Air density =====================================================
    write(*,*)
    write(*,"(a)") "--- semix_air_density ---"

    rhoa = semix_air_density(T_AIR,P0,c)
    call check_close("rho_a(263.15 K, 101325 Pa)", rhoa, 1.3413545395_wp, 1.0e-6_wp, nfail)

    ! === Saturation humidity =============================================
    ! At T0 the SEMIX exponent vanishes, leaving 380.1726/p; the BESSI variant
    ! leaves 0.622*611.2/p. The two forms differ by their fit coefficients.
    write(*,*)
    write(*,"(a)") "--- semix_q_sat / semix_dqsat_dT ---"

    c%semix_qsat = SEMIX_QSAT_SEMIX
    call check_close("q_sat semix at T0", semix_q_sat(c%T0,P0,c), &
                     3.7520118431e-3_wp, 1.0e-6_wp, nfail)
    call check_close("q_sat semix at 263.15 K", semix_q_sat(263.15_wp,P0,c), &
                     1.5940374185e-3_wp, 1.0e-6_wp, nfail)

    c%semix_qsat = SEMIX_QSAT_BESSI
    call check_close("q_sat bessi at T0", semix_q_sat(c%T0,P0,c), &
                     3.7519506538e-3_wp, 1.0e-6_wp, nfail)
    call check_close("q_sat bessi at 263.15 K", semix_q_sat(263.15_wp,P0,c), &
                     1.5952776445e-3_wp, 1.0e-6_wp, nfail)

    ! The two variants must stay close: the option exists to test sensitivity,
    ! not to change the answer.
    c%semix_qsat = SEMIX_QSAT_SEMIX
    q_a = semix_q_sat(263.15_wp,P0,c)
    c%semix_qsat = SEMIX_QSAT_BESSI
    call check("the two q_sat variants agree within 1%", &
               abs(semix_q_sat(263.15_wp,P0,c) - q_a) .lt. 0.01_wp*q_a, nfail)

    ! Derivatives against a central difference of the humidity itself. The
    ! step is 0.5 K: q_sat is sp and varies by only ~2%/K, so a smaller step
    ! is swamped by rounding, while the O(h^2) truncation error at 0.5 K is
    ! only a few tenths of a percent.
    !
    ! The BESSI pair is analytically consistent, so it matches to truncation
    ! error. The SEMIX pair is NOT quite: its dqsat_dT_i is built from the
    ! 273.86-offset e_sat_i while q_sat_i uses the T+0.71 approximation. That
    ! ~0.2% mismatch is CLIMBER-X's, reproduced deliberately, and the looser
    ! tolerance below is what pins it.
    c%semix_qsat = SEMIX_QSAT_BESSI
    dq_fd = (semix_q_sat(263.65_wp,P0,c) - semix_q_sat(262.65_wp,P0,c))/1.0_wp
    call check_close("dq/dT bessi matches a central difference", &
                     semix_dqsat_dT(263.15_wp,P0,c), dq_fd, 5.0e-3_wp, nfail)

    c%semix_qsat = SEMIX_QSAT_SEMIX
    dq_fd = (semix_q_sat(263.65_wp,P0,c) - semix_q_sat(262.65_wp,P0,c))/1.0_wp
    call check_close("dq/dT semix matches a central difference to 1%", &
                     semix_dqsat_dT(263.15_wp,P0,c), dq_fd, 1.0e-2_wp, nfail)

    ! === Air humidity ====================================================
    write(*,*)
    write(*,"(a)") "--- semix_air_humidity ---"

    call check_close("rh = 1 gives saturation at t2m", &
                     semix_air_humidity(T_AIR,1.0_wp,P0,c), &
                     semix_q_sat(T_AIR,P0,c), 1.0e-6_wp, nfail)
    call check_close("rh = 80 (percent) equals rh = 0.8", &
                     semix_air_humidity(T_AIR,80.0_wp,P0,c), &
                     semix_air_humidity(T_AIR,0.8_wp,P0,c), 1.0e-6_wp, nfail)

    ! === Exchange coefficients ===========================================
    write(*,*)
    write(*,"(a)") "--- semix_turbulent_exchange ---"

    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_COLD,WIND,P0,0.8_wp,.TRUE.)

    call check_close("r_a matches semix_resistance", x%r_a, r_stable, 1.0e-6_wp, nfail)
    call check_close("f_sh = rho_a*cp_air/r_a", x%f_sh, rhoa*c%cp_air/r_stable, &
                     1.0e-6_wp, nfail)
    call check_close("f_lh = (Lv+Lm)*rho_a/r_a", x%f_lh, (c%Lv + c%Lm)/r_stable*rhoa, &
                     1.0e-6_wp, nfail)
    call check_close("qsat is taken at the surface temperature", x%qsat, &
                     semix_q_sat(T_COLD,P0,c), 1.0e-6_wp, nfail)
    call check_close("q_air is taken at the air temperature", x%q_air, &
                     semix_air_humidity(T_AIR,0.8_wp,P0,c), 1.0e-6_wp, nfail)

    ! Stable f_sh sits well below BESSI's fixed D_sh; unstable well above it.
    ! This is the whole point of the scheme, so it is pinned rather than left
    ! as a remark in a comment.
    call check("stable f_sh < D_sh", x%f_sh .lt. c%D_sh, nfail)
    x_warm = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_WARM,WIND,P0,0.8_wp,.TRUE.)
    call check("unstable f_sh > D_sh", x_warm%f_sh .gt. c%D_sh, nfail)

    ! No humidity forcing: the latent flux vanishes entirely, exactly as the
    ! BESSI three-way selection does.
    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_COLD,WIND,P0,0.8_wp,.FALSE.)
    call check("no humidity forcing zeroes f_lh", x%f_lh .eq. 0.0_wp, nfail)
    call check("no humidity forcing zeroes q_air", x%q_air .eq. 0.0_wp, nfail)
    call check("no humidity forcing zeroes the latent flux", &
               semix_latent_heat_flux(x) .eq. 0.0_wp, nfail)
    call check("no humidity forcing leaves f_sh alone", x%f_sh .gt. 0.0_wp, nfail)

    ! Dew inhibition: with a surface colder than saturated air, q_air > qsat,
    ! which is the deposition condition l_dew switches off.
    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_COLD,WIND,P0,1.0_wp,.TRUE.)
    call check("saturated air over a cold surface would deposit", &
               x%q_air .gt. x%qsat, nfail)
    call check("l_dew = T allows deposition", semix_latent_heat_flux(x) .gt. 0.0_wp, nfail)

    c%l_dew = .FALSE.
    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_COLD,WIND,P0,1.0_wp,.TRUE.)
    call check("l_dew = F zeroes f_lh when q_air > qsat", x%f_lh .eq. 0.0_wp, nfail)
    call check("l_dew = F zeroes the deposition flux", &
               semix_latent_heat_flux(x) .eq. 0.0_wp, nfail)

    ! ... but must NOT suppress sublimation, which is the other sign.
    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_WARM,WIND,P0,0.4_wp,.TRUE.)
    call check("l_dew = F leaves sublimation untouched", &
               semix_latent_heat_flux(x) .lt. 0.0_wp, nfail)
    c%l_dew = .TRUE.

    ! === Sign convention =================================================
    ! chion counts fluxes positive INTO the surface; SEMIX's own diagnostics
    ! count sensible and latent positive away from it.
    write(*,*)
    write(*,"(a)") "--- sign convention (positive into the surface) ---"

    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_COLD,WIND,P0,0.8_wp,.TRUE.)
    call check("warm air over a cold surface heats it", &
               semix_sensible_heat_flux(x,T_AIR,T_COLD) .gt. 0.0_wp, nfail)
    call check_close("Q_sh = f_sh*(T_air - T_s)", &
                     semix_sensible_heat_flux(x,T_AIR,T_COLD), &
                     x%f_sh*(T_AIR - T_COLD), 1.0e-6_wp, nfail)

    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_WARM,WIND,P0,0.4_wp,.TRUE.)
    call check("dry air over a warm surface sublimates (cools it)", &
               semix_latent_heat_flux(x) .lt. 0.0_wp, nfail)

    ! === The coupling identity ===========================================
    ! smb_ebal.f90:122-123 splits the latent flux into
    !     num_lh   = -f_lh*(qsat - dqsatdT*T_s - q_air)
    !     denom_lh =  f_lh*dqsatdT
    ! so that the flux is num_lh - denom_lh*T. Those two are exactly what
    ! snow_energy will add to q_const and q_lin (coupling decision alpha), and
    ! at the linearization point the split must reproduce the exact flux.
    write(*,*)
    write(*,"(a)") "--- ebal num_lh/denom_lh -> chion q_const/q_lin ---"

    x = semix_turbulent_exchange(c,H_DEEP,T_AIR,T_COLD,WIND,P0,0.8_wp,.TRUE.)

    num_lh   = -x%f_lh*(x%qsat - x%dqsatdT*T_COLD - x%q_air)
    denom_lh =  x%f_lh*x%dqsatdT

    ! 1e-4 rather than machine tolerance: num_lh and denom_lh*T_s are both
    ! near 1100 W m-2 and cancel down to about 6, which costs four digits in
    ! single precision. The identity itself is exact.
    call check_close("num_lh - denom_lh*T_s = exact latent flux", &
                     num_lh - denom_lh*T_COLD, semix_latent_heat_flux(x), &
                     1.0e-4_wp, nfail)
    call check("denom_lh > 0 (a warmer surface loses more vapour)", &
               denom_lh .gt. 0.0_wp, nfail)

    ! The sensible split is trivial but is pinned for the same reason.
    call check_close("num_sh - denom_sh*T_s = exact sensible flux", &
                     x%f_sh*T_AIR - x%f_sh*T_COLD, &
                     semix_sensible_heat_flux(x,T_AIR,T_COLD), 1.0e-4_wp, nfail)

    ! === Summary =========================================================
    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " snow_seb_semix: ALL CHECKS PASSED"
    else
        write(*,"(a,i0,a)") " snow_seb_semix: ", nfail, " CHECK(S) FAILED"
    end if
    write(*,"(a)") "=========================================================="

    if (nfail .gt. 0) stop 1

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

    subroutine check_close(label,value,expected,rtol,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value
        real(wp),         intent(IN)    :: expected
        real(wp),         intent(IN)    :: rtol
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp) :: tol

        tol = max(abs(expected)*rtol,1.0e-20_wp)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_close

end program test_seb
