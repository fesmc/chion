program test_monthly
    ! Acceptance test: chion_forcing_monthly.
    !
    ! Covers:
    !   1. Interpolation stencil sanity: weights are a partition of unity, the
    !      mid-month day returns the month's control value exactly, and the
    !      series is continuous across month boundaries.
    !   2. Mean preservation: for arbitrary monthly means, the daily series
    !      built from the mean-preserving controls averages back to the input
    !      means per month (and hence the annual total is conserved) to
    !      machine precision.
    !   3. Exactness on a constant field: a flat input yields a flat daily
    !      series with controls equal to the input.
    !   4. Column-list path (2D controls) agrees with the per-column path.

    use chion_defs, only : wp, dp
    use chion_forcing_monthly

    implicit none

    type(monthly_to_daily_class) :: md
    integer, parameter :: nmon = 12
    integer, parameter :: ndm  = 30
    integer, parameter :: nday = nmon*ndm
    integer            :: nfail

    nfail = 0

    call monthly_to_daily_init(md, nmon, ndm)

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion acceptance test: chion_forcing_monthly"
    write(*,"(a)") "=========================================================="
    write(*,*)

    call test_stencil(nfail)
    call test_mean_preservation(nfail)
    call test_constant_field(nfail)
    call test_column_path(nfail)

    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " monthly: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " monthly: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    subroutine report(name, ok, val, tol, nfail)
        character(len=*), intent(IN)    :: name
        logical,          intent(IN)    :: ok
        real(wp),         intent(IN)    :: val, tol
        integer,          intent(INOUT) :: nfail
        character(len=4) :: tag
        tag = "PASS"
        if (.not. ok) then
            tag = "FAIL"
            nfail = nfail + 1
        end if
        write(*,"(a,a,a,es12.4,a,es10.2,a)") "   [", tag, "] "//name// &
            "  (max ", val, " vs tol ", tol, ")"
        return
    end subroutine report

    subroutine test_stencil(nfail)
        integer, intent(INOUT) :: nfail
        real(wp) :: ctrl(nmon), vday, vprev, err, bound
        integer  :: d, q, mid, qn

        write(*,"(a)") " 1. Interpolation stencil"

        ! Partition of unity: wt0 + wt1 == 1 everywhere.
        err = 0.0_wp
        do d = 1, nday
            err = max(err, abs((md%wt0(d) + md%wt1(d)) - 1.0_wp))
        end do
        call report("weights sum to one", err .le. 1.0e-6_wp, err, 1.0e-6_wp, nfail)

        ! Mid-month day returns the month's own control value. Use a distinct
        ! value per month so any leakage from a neighbour shows up.
        do q = 1, nmon
            ctrl(q) = real(q, wp)
        end do
        mid = ndm/2
        err = 0.0_wp
        do q = 1, nmon
            d = (q-1)*ndm + mid
            call interp_monthly_to_day(md, ctrl, d, vday)
            err = max(err, abs(vday - real(q, wp)))
        end do
        call report("mid-month day is exact", err .le. 1.0e-6_wp, err, 1.0e-6_wp, nfail)

        ! Continuity: no day-to-day jump exceeds the sharpest ramp the stencil
        ! can produce, which is the largest adjacent control step (including
        ! the Dec->Jan wrap, an 11-unit step for this ctrl(q)=q sawtooth) spread
        ! over one month. Derive the bound rather than hardcode it.
        bound = 0.0_wp
        do q = 1, nmon
            qn = q + 1
            if (qn .eq. nmon+1) qn = 1
            bound = max(bound, abs(ctrl(qn) - ctrl(q)))
        end do
        bound = bound/real(ndm, wp) + 1.0e-4_wp

        call interp_monthly_to_day(md, ctrl, nday, vprev)
        err = 0.0_wp
        do d = 1, nday
            call interp_monthly_to_day(md, ctrl, d, vday)
            err  = max(err, abs(vday - vprev))
            vprev = vday
        end do
        call report("series is continuous", err .le. bound, err, bound, nfail)

        write(*,*)
        return
    end subroutine test_stencil

    subroutine test_mean_preservation(nfail)
        integer, intent(INOUT) :: nfail
        real(wp) :: means(nmon), ctrl(nmon), rebuilt(nmon)
        real(wp) :: vday, err, tot_in, tot_out
        integer  :: q, k, d

        write(*,"(a)") " 2. Mean preservation"

        ! A deliberately non-smooth monthly series (a spring precip peak) so
        ! the off-diagonal correction is genuinely exercised.
        means = [ 0.5_wp, 0.6_wp, 1.2_wp, 3.8_wp, 6.1_wp, 2.0_wp, &
                  1.1_wp, 0.9_wp, 1.5_wp, 2.7_wp, 1.0_wp, 0.4_wp ]

        call monthly_to_daily_controls(md, means, ctrl)

        ! Average the daily series back over each month.
        rebuilt = 0.0_wp
        do q = 1, nmon
            do k = 1, ndm
                d = (q-1)*ndm + k
                call interp_monthly_to_day(md, ctrl, d, vday)
                rebuilt(q) = rebuilt(q) + vday
            end do
            rebuilt(q) = rebuilt(q) / real(ndm, wp)
        end do

        err = maxval(abs(rebuilt - means))
        call report("monthly means reproduced", err .le. 1.0e-4_wp, err, 1.0e-4_wp, nfail)

        tot_in  = sum(means)
        tot_out = sum(rebuilt)
        call report("annual total conserved", abs(tot_out-tot_in) .le. 1.0e-4_wp, &
                    abs(tot_out-tot_in), 1.0e-4_wp, nfail)

        write(*,*)
        return
    end subroutine test_mean_preservation

    subroutine test_constant_field(nfail)
        integer, intent(INOUT) :: nfail
        real(wp) :: means(nmon), ctrl(nmon), vday, err
        integer  :: d

        write(*,"(a)") " 3. Constant field"

        means = 4.2_wp
        call monthly_to_daily_controls(md, means, ctrl)

        err = maxval(abs(ctrl - 4.2_wp))
        call report("controls equal input", err .le. 1.0e-4_wp, err, 1.0e-4_wp, nfail)

        err = 0.0_wp
        do d = 1, nday
            call interp_monthly_to_day(md, ctrl, d, vday)
            err = max(err, abs(vday - 4.2_wp))
        end do
        call report("daily series is flat", err .le. 1.0e-4_wp, err, 1.0e-4_wp, nfail)

        write(*,*)
        return
    end subroutine test_constant_field

    subroutine test_column_path(nfail)
        integer, intent(INOUT) :: nfail
        integer, parameter :: ncol = 5
        real(wp) :: means2(ncol,nmon), ctrl2(ncol,nmon)
        real(wp) :: means1(nmon), ctrl1(nmon)
        real(wp) :: vday2(ncol), vday1, err
        integer  :: i, q, d

        write(*,"(a)") " 4. Column-list path matches per-column"

        do i = 1, ncol
            do q = 1, nmon
                means2(i,q) = real(i, wp) + 0.3_wp*real(q, wp) - 0.02_wp*real(q*q, wp)
            end do
        end do

        call monthly_to_daily_controls(md, means2, ctrl2)

        err = 0.0_wp
        do i = 1, ncol
            means1 = means2(i,:)
            call monthly_to_daily_controls(md, means1, ctrl1)
            err = max(err, maxval(abs(ctrl1 - ctrl2(i,:))))
        end do
        call report("2D controls match 1D", err .le. 1.0e-5_wp, err, 1.0e-5_wp, nfail)

        err = 0.0_wp
        do d = 1, nday, 7
            call interp_monthly_to_day(md, ctrl2, d, vday2)
            do i = 1, ncol
                call interp_monthly_to_day(md, ctrl2(i,:), d, vday1)
                err = max(err, abs(vday2(i) - vday1))
            end do
        end do
        call report("2D interp matches 1D", err .le. 1.0e-5_wp, err, 1.0e-5_wp, nfail)

        write(*,*)
        return
    end subroutine test_column_path

end program test_monthly
