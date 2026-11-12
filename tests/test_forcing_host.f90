program test_forcing_host
    ! Acceptance test: chion_forcing_host.
    !
    !   1. forcing_shortwave_transmissivity: TOA*(trans0 + trans_per_km*z[km]),
    !      clamped to [0, 1]; non-finite or negative height = sea level.
    !   2. forcing_wind_speed: elevation term; a uniform slope gives u_s*slope
    !      inside the grid, the same on a grid twice as coarse; at least 0.
    !   3. forcing_tvar_sigma: the line and its bounds.
    !   4. forcing_red_noise_monthly: zero mean in every month, unit variance,
    !      lag-one correlation near exp(-dt/tau), reproducible from the seed,
    !      and the Park-Miller start (seed 1: x = 33614, 564950498).

    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

    use chion_defs, only : wp
    use chion_forcing_host

    implicit none

    integer :: nfail

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion acceptance test: chion_forcing_host"
    write(*,"(a)") "=========================================================="

    call test_shortwave(nfail)
    call test_wind(nfail)
    call test_sigma(nfail)
    call test_noise(nfail)

    write(*,*)
    if (nfail .eq. 0) then
        write(*,"(a)") " chion_forcing_host: ALL CHECKS PASSED"
    else
        write(*,"(a,i0,a)") " chion_forcing_host: ", nfail, " CHECK(S) FAILED"
        error stop 1
    end if

contains

    subroutine test_shortwave(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: nan_wp

        write(*,*); write(*,"(a)") "--- shortwave through a transmissivity ---"
        nan_wp = ieee_value(nan_wp, ieee_quiet_nan)
        call check_close("Greenland: SWdn = TOA (0.587 + 0.052*1.5)", &
                         forcing_shortwave_transmissivity(500.0_wp, 1500.0_wp, 0.587_wp, 0.052_wp), &
                         500.0_wp*(0.587_wp + 0.052_wp*1.5_wp), 1.0e-6_wp, nfail)
        call check_close("negative height = sea level", &
                         forcing_shortwave_transmissivity(500.0_wp, -50.0_wp, 0.587_wp, 0.052_wp), &
                         500.0_wp*0.587_wp, 1.0e-6_wp, nfail)
        call check_close("NaN height = sea level", &
                         forcing_shortwave_transmissivity(500.0_wp, nan_wp, 0.587_wp, 0.052_wp), &
                         500.0_wp*0.587_wp, 1.0e-6_wp, nfail)
        call check_close("transmissivity capped at 1", &
                         forcing_shortwave_transmissivity(500.0_wp, 9000.0_wp, 0.68_wp, 0.041_wp), &
                         500.0_wp, 1.0e-6_wp, nfail)
        call check_close("transmissivity floored at 0", &
                         forcing_shortwave_transmissivity(500.0_wp, 1000.0_wp, -0.2_wp, 0.1_wp), &
                         0.0_wp, 1.0e-6_wp, nfail)
    end subroutine test_shortwave

    subroutine test_wind(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: z(40,30), z2(20,15), u(40,30), u2(20,15)
        integer  :: i

        write(*,*); write(*,"(a)") "--- wind from elevation and slope ---"
        z = 2000.0_wp
        u = forcing_wind_speed(z, 16.0e3_wp, 2.5_wp, 1.27_wp, 0.0_wp, 50.0_wp)
        call check_close("flat: u0 + u_z z[km]", u(10,10), 2.5_wp + 1.27_wp*2.0_wp, 1.0e-9_wp, nfail)
        do i = 1, 40
            z(i,:) = 20.0_wp*16.0_wp*real(i-1, wp)          ! 20 m/km along x
        end do
        do i = 1, 20
            z2(i,:) = 20.0_wp*32.0_wp*real(i-1, wp)
        end do
        u  = forcing_wind_speed(z,  16.0e3_wp, 3.0_wp, 0.0_wp, 0.2_wp, 50.0_wp)
        u2 = forcing_wind_speed(z2, 32.0e3_wp, 3.0_wp, 0.0_wp, 0.2_wp, 50.0_wp)
        call check_close("uniform slope 20 m/km: u0 + u_s*20 inside", u(20,15), 3.0_wp + 0.2_wp*20.0_wp, 1.0e-9_wp, nfail)
        call check_close("the same on a grid twice as coarse", u2(10,8), u(20,15), 1.0e-9_wp, nfail)
        call check_close("one-sided at the edge", u(1,1), 3.0_wp + 0.2_wp*20.0_wp, 1.0e-9_wp, nfail)
        u = forcing_wind_speed(0.0_wp*z, 16.0e3_wp, -1.0_wp, 0.0_wp, 0.0_wp, 50.0_wp)
        call check("never negative", all(u .eq. 0.0_wp), nfail)
    end subroutine test_wind

    subroutine test_sigma(nfail)
        integer, intent(inout) :: nfail
        write(*,*); write(*,"(a)") "--- sigma of the daily temperature ---"
        call check_close("sigma(T0 - 10) = 2 + 1.4", forcing_tvar_sigma(263.15_wp, 273.15_wp, 2.0_wp, -0.14_wp, 2.0_wp, 6.0_wp), &
                         3.4_wp, 1.0e-9_wp, nfail)
        call check_close("bounded below at T > T0", forcing_tvar_sigma(278.15_wp, 273.15_wp, 2.0_wp, -0.14_wp, 2.0_wp, 6.0_wp), &
                         2.0_wp, 1.0e-12_wp, nfail)
        call check_close("bounded above in the cold", forcing_tvar_sigma(223.15_wp, 273.15_wp, 2.0_wp, -0.14_wp, 2.0_wp, 6.0_wp), &
                         6.0_wp, 1.0e-12_wp, nfail)
    end subroutine test_sigma

    subroutine test_noise(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: xi(360), xb(360), r1
        integer  :: m

        write(*,*); write(*,"(a)") "--- red noise, zero monthly mean ---"
        xi = forcing_red_noise_monthly(360, 1, 30, 4.0_wp, 1)
        call check("zero mean in every month", all([(abs(sum(xi(30*m+1:30*m+30)))/30.0_wp .lt. 1.0e-5_wp, m = 0, 11)]), nfail)
        call check_close("unit variance", sum(xi*xi)/360.0_wp, 1.0_wp, 1.0e-5_wp, nfail)
        r1 = sum(xi(1:359)*xi(2:360))/sum(xi*xi)
        write(*,"(a,f6.3,a,f6.3)") "        lag-1 correlation ", r1, ", exp(-1/4) = ", exp(-0.25_wp)
        call check("lag-1 correlation within 0.15 of exp(-dt/tau)", abs(r1 - exp(-0.25_wp)) .lt. 0.15_wp, nfail)
        xb = forcing_red_noise_monthly(360, 1, 30, 4.0_wp, 1)
        call check("the same seed, the same sequence", all(xb .eq. xi), nfail)
        xb = forcing_red_noise_monthly(360, 1, 30, 4.0_wp, 2)
        call check("another seed, another sequence", any(xb .ne. xi), nfail)
        xb = forcing_red_noise_monthly(360, 5, 30, 4.0_wp, 1)
        call check("dt = 5 d: 72 values of zero monthly mean", &
                   all([(abs(sum(xb(6*m+1:6*m+6)))/6.0_wp .lt. 1.0e-5_wp, m = 0, 11)]), nfail)
        write(*,"(a,3f12.7)") "        seed 1, tau 4 d, dt 1 d: xi(1:3) = ", xi(1:3)
    end subroutine test_noise

    subroutine check(label, ok, nfail)
        character(len=*), intent(in) :: label
        logical, intent(in) :: ok
        integer, intent(inout) :: nfail
        if (ok) then
            write(*,"(a,a)") "  ok   : ", trim(label)
        else
            write(*,"(a,a)") "  FAIL : ", trim(label); nfail = nfail + 1
        end if
    end subroutine check

    subroutine check_close(label, val, expect, tol, nfail)
        character(len=*), intent(in) :: label
        real(wp), intent(in) :: val, expect, tol
        integer, intent(inout) :: nfail
        if (abs(val - expect) .le. tol*max(abs(expect), 1.0_wp)) then
            write(*,"(a,a,a,g14.7)") "  ok   : ", trim(label), " = ", val
        else
            write(*,"(a,a,2g16.8)") "  FAIL : ", trim(label), val, expect; nfail = nfail + 1
        end if
    end subroutine check_close

end program test_forcing_host
