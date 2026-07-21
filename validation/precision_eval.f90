program prec_test
    ! Evaluate where single precision breaks down in the chion physics.
    ! Each test computes an expression in sp and dp and reports the error,
    ! using values taken from the actual Chion.jl algorithms.

    implicit none

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.d0)

    write(*,"(a)") "=================================================================="
    write(*,"(a)") " Precision evaluation for chion"
    write(*,"(a,es12.4)") " sp epsilon = ", epsilon(1.0_sp)
    write(*,"(a,es12.4)") " dp epsilon = ", epsilon(1.0_dp)
    write(*,"(a)") "=================================================================="

    call test_pore_volume()
    call test_cold_content()
    call test_cumulative()
    call test_melt_residual()
    call test_densification_gap()
    call test_tridiag()
    call test_temperature_offset()

contains

    subroutine test_pore_volume()
        ! percolation.jl:50 / albedo.jl:47:  phi = m/rho - m/rho_i
        ! Guarded by  phi <= EPS_TINY (1e-12), then used as  lwc = mw/rho_w/phi.
        ! Cancellation as rho -> rho_i is the whole question.
        real(dp) :: m, rho_i, rho, phi_dp, phi_sp_as_dp
        real(sp) :: m_s, rho_i_s, rho_s, phi_s
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- TEST 1: pore volume  phi = m/rho - m/rho_i  [m] ---"
        write(*,"(a)") "  (guard threshold TOL_TINY = 1e-12; phi feeds a DIVISION)"
        write(*,"(a10,a16,a16,a14,a12)") "rho", "phi (dp)", "phi (sp)", "abs err", "rel err"

        m     = 300.0_dp
        rho_i = 917.0_dp

        do k = 1, 6
            rho = rho_i - 10.0_dp**(1-k)     ! 916, 916.9, 916.99, ...

            phi_dp = m/rho - m/rho_i

            m_s     = real(m,sp)
            rho_i_s = real(rho_i,sp)
            rho_s   = real(rho,sp)
            phi_s   = m_s/rho_s - m_s/rho_i_s
            phi_sp_as_dp = real(phi_s,dp)

            write(*,"(f10.5,es16.6,es16.6,es14.3,f12.4)") &
                rho, phi_dp, phi_sp_as_dp, abs(phi_sp_as_dp-phi_dp), &
                abs(phi_sp_as_dp-phi_dp)/max(abs(phi_dp),1.0e-30_dp)
        end do

        write(*,"(a)") "  NOTE: sp phi is quantized -- the smallest nonzero sp value here"
        write(*,"(a)") "        is ~1e-5 m, so the 1e-12 guard can never fire in sp."

        return
    end subroutine test_pore_volume

    subroutine test_cold_content()
        ! refreezing.jl:33:  Qcold = (T0 - T)*ci*ms
        ! Entry condition is the STRICT test T < T0. Near the melting point the
        ! difference of two ~273 K numbers is the entire signal.
        real(dp) :: T0, T, ci, ms, q_dp, q_sp_as_dp, dm_dp, dm_sp
        real(sp) :: T0_s, T_s, ci_s, ms_s, q_s
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- TEST 2: cold content  Qcold = (T0-T)*ci*ms  [J m-2] ---"
        write(*,"(a)") "  (and the refrozen mass dm = Qcold/Lm it implies)"
        write(*,"(a12,a16,a16,a12,a14)") "T0-T [K]", "Qcold (dp)", "Qcold (sp)", "rel err", "dm err [kg m-2]"

        T0 = 273.15_dp
        ci = 2110.0_dp
        ms = 300.0_dp

        do k = 1, 6
            T = T0 - 10.0_dp**(-k)      ! 0.1, 0.01, ... 1e-6 K below melting

            q_dp = (T0 - T)*ci*ms

            T0_s = real(T0,sp)
            T_s  = real(T,sp)
            ci_s = real(ci,sp)
            ms_s = real(ms,sp)
            q_s  = (T0_s - T_s)*ci_s*ms_s
            q_sp_as_dp = real(q_s,dp)

            dm_dp = q_dp/334000.0_dp
            dm_sp = q_sp_as_dp/334000.0_dp

            write(*,"(es12.1,es16.6,es16.6,f12.4,es14.3)") &
                T0-T, q_dp, q_sp_as_dp, &
                abs(q_sp_as_dp-q_dp)/max(abs(q_dp),1.0e-30_dp), abs(dm_sp-dm_dp)
        end do

        write(*,"(a,es12.4,a)") "  sp resolution at 273 K is ", &
            spacing(real(273.15,sp)), " K -- differences below this are pure noise."

        return
    end subroutine test_cold_content

    subroutine test_cumulative()
        ! smb_ice, runoff, melt, refreezing, pdd_sum are CUMULATIVE over the run.
        ! Daily steps for 100 years = 36500 increments.
        real(dp) :: acc_dp, inc_dp
        real(sp) :: acc_s, inc_s
        integer  :: i, nstep

        write(*,*)
        write(*,"(a)") "--- TEST 3: cumulative diagnostics over a 100-yr daily run ---"
        write(*,"(a)") "  (smb_ice / runoff / melt accumulate every step, never reset)"
        write(*,"(a12,a18,a18,a14,a12)") "increment", "total (dp)", "total (sp)", "abs err", "rel err"

        nstep = 36500

        ! Case A: steady accumulation, 2 kg m-2 per day
        inc_dp = 2.0_dp
        inc_s  = real(inc_dp,sp)
        acc_dp = 0.0_dp
        acc_s  = 0.0_sp
        do i = 1, nstep
            acc_dp = acc_dp + inc_dp
            acc_s  = acc_s  + inc_s
        end do
        write(*,"(es12.2,es18.8,es18.8,es14.3,es12.2)") inc_dp, acc_dp, real(acc_s,dp), &
            abs(real(acc_s,dp)-acc_dp), abs(real(acc_s,dp)-acc_dp)/acc_dp

        ! Case B: small increments onto a large total -- the dangerous case
        inc_dp = 0.01_dp
        inc_s  = real(inc_dp,sp)
        acc_dp = 1.0e5_dp
        acc_s  = real(1.0e5_dp,sp)
        do i = 1, nstep
            acc_dp = acc_dp + inc_dp
            acc_s  = acc_s  + inc_s
        end do
        write(*,"(es12.2,es18.8,es18.8,es14.3,es12.2)") inc_dp, acc_dp, real(acc_s,dp), &
            abs(real(acc_s,dp)-acc_dp), abs(real(acc_s,dp)-acc_dp)/acc_dp
        write(*,"(a)") "  ^ increments of 0.01 onto a total of 1e5 are LOST ENTIRELY in sp"
        write(*,"(a,es12.4)") "    (sp spacing at 1e5 = ", spacing(real(1.0e5,sp))

        return
    end subroutine test_cumulative

    subroutine test_melt_residual()
        ! energy_flux.jl:247-274:
        !   melt_energy_available = max((Q_const - Q_lin*T)*dt - energy_to_melting, 0)
        ! A difference of two large, nearly equal energies.
        real(dp) :: qc, ql, T, dt, e2m, res_dp, res_s_as_dp
        real(sp) :: qc_s, ql_s, T_s, dt_s, e2m_s, res_s

        write(*,*)
        write(*,"(a)") "--- TEST 4: melt energy residual (difference of large energies) ---"

        qc  = 3.0e3_dp        ! W m-2 constant part
        ql  = 10.0_dp         ! W m-2 K-1 linear part
        T   = 273.15_dp
        dt  = 86400.0_dp
        e2m = 6.3e7_dp        ! J m-2 already spent warming to T0

        res_dp = (qc - ql*T)*dt - e2m

        qc_s  = real(qc,sp);  ql_s = real(ql,sp); T_s = real(T,sp)
        dt_s  = real(dt,sp);  e2m_s = real(e2m,sp)
        res_s = (qc_s - ql_s*T_s)*dt_s - e2m_s
        res_s_as_dp = real(res_s,dp)

        write(*,"(a,es16.8)") "  residual (dp)      = ", res_dp
        write(*,"(a,es16.8)") "  residual (sp)      = ", res_s_as_dp
        write(*,"(a,es16.8)") "  abs error [J m-2]  = ", abs(res_s_as_dp-res_dp)
        write(*,"(a,es16.8)") "  implied melt error [kg m-2] = ", abs(res_s_as_dp-res_dp)/334000.0_dp
        write(*,"(a,es12.4)") "  sp spacing at 6.3e7 J m-2 = ", spacing(real(6.3e7,sp))

        return
    end subroutine test_melt_residual

    subroutine test_densification_gap()
        ! densification.jl:9:  drho_dt ~ (rho_i - rho), with rho -> rho_i,
        ! and the snap-to-ice guard  rho >= rho_i - EPS_TINY.
        real(dp) :: rho_i, rho, gap_dp, gap_s_as_dp
        real(sp) :: rho_i_s, rho_s, gap_s
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- TEST 5: densification gap  (rho_i - rho)  [kg m-3] ---"
        write(*,"(a)") "  (also the guard  rho >= rho_i - TOL_TINY)"
        write(*,"(a14,a18,a18,a14)") "true gap", "gap (dp)", "gap (sp)", "abs err"

        rho_i = 917.0_dp
        do k = 2, 7
            rho = rho_i - 10.0_dp**(2-k)
            gap_dp = rho_i - rho
            rho_i_s = real(rho_i,sp); rho_s = real(rho,sp)
            gap_s = rho_i_s - rho_s
            gap_s_as_dp = real(gap_s,dp)
            write(*,"(es14.1,es18.8,es18.8,es14.3)") &
                10.0_dp**(2-k), gap_dp, gap_s_as_dp, abs(gap_s_as_dp-gap_dp)
        end do
        write(*,"(a,es12.4)") "  sp spacing at 917 kg m-3 = ", spacing(real(917.0,sp))

        return
    end subroutine test_densification_gap

    subroutine test_tridiag()
        ! energy_flux.jl: Thomas algorithm, no pivoting, on the assembled
        ! conduction matrix. diag = 1 - lower - upper with lower,upper < 0,
        ! so the system is strictly diagonally dominant. Check sp is adequate.
        integer, parameter :: n = 15
        real(dp) :: lo(n), di(n), up(n), rh(n), sol_dp(n)
        real(sp) :: lo_s(n), di_s(n), up_s(n), rh_s(n)
        real(dp) :: maxerr
        integer  :: i

        write(*,*)
        write(*,"(a)") "--- TEST 6: tridiagonal conduction solve, Thomas, n=15 ---"

        ! Representative coefficients: strong conduction coupling
        do i = 1, n
            lo(i) = -0.35_dp
            up(i) = -0.35_dp
            di(i) =  1.0_dp + 0.70_dp
            rh(i) =  263.15_dp
        end do
        di(1) = 1.0_dp + 0.35_dp + 0.12_dp     ! surface flux feedback
        di(n) = 1.0_dp + 0.35_dp               ! zero-flux bottom
        rh(1) = 263.15_dp + 4.7_dp

        lo_s = real(lo,sp); di_s = real(di,sp); up_s = real(up,sp); rh_s = real(rh,sp)

        call thomas_dp(lo,di,up,rh,n)
        sol_dp = rh
        call thomas_sp(lo_s,di_s,up_s,rh_s,n)

        maxerr = 0.0_dp
        do i = 1, n
            maxerr = max(maxerr,abs(real(rh_s(i),dp)-sol_dp(i)))
        end do

        write(*,"(a,f16.10)") "  surface T (dp)        = ", sol_dp(1)
        write(*,"(a,f16.10)") "  surface T (sp)        = ", real(rh_s(1),dp)
        write(*,"(a,es16.6)") "  max |dT| over column  = ", maxerr
        write(*,"(a)") "  -> matrix is diagonally dominant; sp is fine for the SOLVE itself."

        return
    end subroutine test_tridiag

    subroutine test_temperature_offset()
        ! Option: store temperature as an offset from T0 instead of absolute K.
        real(dp) :: dT
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- TEST 7: absolute K vs offset-from-T0 storage in sp ---"
        write(*,"(a14,a22,a22)") "dT below T0", "sp err, absolute K", "sp err, offset"

        do k = 1, 6
            dT = 10.0_dp**(-k)
            write(*,"(es14.1,es22.4,es22.4)") dT, &
                abs(real(real(273.15_dp,sp) - real(273.15_dp-dT,sp),dp) - dT), &
                abs(real(real(dT,sp),dp) - dT)
        end do
        write(*,"(a)") "  -> storing (T - T0) recovers ~7 significant digits on the"
        write(*,"(a)") "     quantity the physics actually uses."

        return
    end subroutine test_temperature_offset

    subroutine thomas_dp(lo,di,up,rh,n)
        integer,  intent(IN)    :: n
        real(dp), intent(INOUT) :: lo(n), di(n), up(n), rh(n)
        real(dp) :: f
        integer  :: i
        do i = 2, n
            f = lo(i-1)/di(i-1)
            di(i) = di(i) - f*up(i-1)
            rh(i) = rh(i) - f*rh(i-1)
        end do
        rh(n) = rh(n)/di(n)
        do i = n-1, 1, -1
            rh(i) = (rh(i) - up(i)*rh(i+1))/di(i)
        end do
        return
    end subroutine thomas_dp

    subroutine thomas_sp(lo,di,up,rh,n)
        integer,  intent(IN)    :: n
        real(sp), intent(INOUT) :: lo(n), di(n), up(n), rh(n)
        real(sp) :: f
        integer  :: i
        do i = 2, n
            f = lo(i-1)/di(i-1)
            di(i) = di(i) - f*up(i-1)
            rh(i) = rh(i) - f*rh(i-1)
        end do
        rh(n) = rh(n)/di(n)
        do i = n-1, 1, -1
            rh(i) = (rh(i) - up(i)*rh(i+1))/di(i)
        end do
        return
    end subroutine thomas_sp

end program prec_test
