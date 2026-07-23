module chion_forcing_monthly
    ! Monthly -> daily forcing interface (library level, model-neutral).
    !
    ! A host commonly has climatological forcing as 12 monthly means but must
    ! drive chion on a daily (or shorter) step over a repeating annual cycle.
    ! This module turns the 12 means into a smooth daily series in a way that
    ! is MEAN-PRESERVING: the daily values, averaged back over each month,
    ! reproduce the input monthly mean exactly. For precipitation that means
    ! the monthly (and hence annual) total is conserved to machine precision,
    ! which a naive linear interpolation between mid-month values is not.
    !
    ! Discretisation. The year is nmon months of nday_mon days each (the chion
    ! default is a 360-day year, 12 x 30). The per-day interpolation stencil is
    ! the one from rembo (~/models/rembo/libs/monthlydaily.f90): each day is a
    ! linear blend of two adjacent monthly CONTROL values, ramping across the
    ! month so that the mid-month day takes the current month alone. That
    ! stencil is NOT mean-preserving on its own; mean preservation is recovered
    ! by first solving for control values whose stencil-average matches the
    ! target means.
    !
    ! The mean-preserving solve. Let w be the stencil: for day d,
    !     value(d) = C(m0(d))*wt0(d) + C(m1(d))*wt1(d)
    ! and let the monthly-mean operator be
    !     Mbar(q) = (1/nday_mon) * sum_{d in month q} value(d).
    ! Mbar is linear in the 12 controls, Mbar = A C, and because the stencil of
    ! any day in month q references only months q-1, q, q+1, A is CYCLIC
    ! TRIDIAGONAL. Given target means M we want controls C with A C = M, i.e.
    ! C = A^{-1} M. A depends only on the stencil, not on the data, so A^{-1} is
    ! built ONCE at init and every column/field is a single mat-vec.
    !
    ! A is assembled directly from the stencil weights rather than from a
    ! hand-derived formula, so it stays correct if nday_mon changes and it is
    ! self-consistent with the interpolation actually used at runtime.
    !
    ! Caveat: the mean-preserving control values can overshoot for a sharply
    ! peaked field (a control may fall below the smallest input mean), so a
    ! non-negative input does not guarantee a non-negative daily series. For
    ! smooth climatological precipitation this is rare; a host that needs a
    ! hard floor should clamp the daily output. Conservation is exact; the
    ! shape is not monotone.

    use chion_defs, only : wp, dp

    implicit none

    private

    type monthly_to_daily_class
        integer :: nmon                          ! months per year   (12)
        integer :: nday_mon                      ! days per month     (30)
        integer :: nday_year                     ! days per year     (360)

        ! rembo interpolation stencil, indexed by day-of-year (1..nday_year)
        integer,  allocatable :: m0(:)           ! first  month referenced
        integer,  allocatable :: m1(:)           ! second month referenced
        real(wp), allocatable :: wt0(:)          ! weight on m0
        real(wp), allocatable :: wt1(:)          ! weight on m1

        ! mean-preserving map: controls = Ainv * monthly_means
        real(wp), allocatable :: Ainv(:,:)       ! (nmon,nmon)
    end type

    interface interp_monthly_to_day
        module procedure interp_monthly_to_day_pt
        module procedure interp_monthly_to_day_1D
    end interface

    interface monthly_to_daily_controls
        module procedure monthly_to_daily_controls_1D
        module procedure monthly_to_daily_controls_2D
    end interface

    public :: monthly_to_daily_class
    public :: monthly_to_daily_init
    public :: monthly_to_daily_controls
    public :: interp_monthly_to_day
    public :: monthly_to_daily_month_of_day

contains

    subroutine monthly_to_daily_init(md, nmon, nday_mon)
        ! Build the interpolation stencil and the mean-preserving map A^{-1}.

        implicit none

        type(monthly_to_daily_class), intent(OUT) :: md
        integer, intent(IN) :: nmon
        integer, intent(IN) :: nday_mon

        ! Local variables
        real(dp), allocatable :: A(:,:)
        integer :: d, q

        md%nmon      = nmon
        md%nday_mon  = nday_mon
        md%nday_year = nmon*nday_mon

        if (allocated(md%m0))   deallocate(md%m0)
        if (allocated(md%m1))   deallocate(md%m1)
        if (allocated(md%wt0))  deallocate(md%wt0)
        if (allocated(md%wt1))  deallocate(md%wt1)
        if (allocated(md%Ainv)) deallocate(md%Ainv)
        allocate(md%m0(md%nday_year), md%m1(md%nday_year))
        allocate(md%wt0(md%nday_year), md%wt1(md%nday_year))
        allocate(md%Ainv(nmon,nmon))

        ! --- Interpolation stencil (rembo kernel) -------------------------
        call interp_weights_monthly_to_daily(md%m0, md%m1, md%wt0, md%wt1, &
                                              nday_mon, nmon)

        ! --- Mean-preserving map ------------------------------------------
        ! A(q,:) = monthly-mean operator applied to the unit controls, i.e.
        ! the stencil weights of the days in month q, accumulated by month and
        ! divided by nday_mon. Done in dp: A is close to the identity plus a
        ! small off-diagonal, and the inverse is reused across every column.
        allocate(A(nmon,nmon))
        A = 0.0_dp
        do d = 1, md%nday_year
            q = monthly_to_daily_month_of_day(md, d)
            A(q, md%m0(d)) = A(q, md%m0(d)) + real(md%wt0(d), dp)
            A(q, md%m1(d)) = A(q, md%m1(d)) + real(md%wt1(d), dp)
        end do
        A = A / real(nday_mon, dp)

        md%Ainv = real(invert_matrix(A), wp)

        deallocate(A)

        return

    end subroutine monthly_to_daily_init

    elemental function monthly_to_daily_month_of_day(md, day) result(q)
        ! Month (1..nmon) containing day-of-year `day` (1..nday_year).

        implicit none

        type(monthly_to_daily_class), intent(IN) :: md
        integer, intent(IN) :: day
        integer :: q

        q = (day-1)/md%nday_mon + 1

        return

    end function monthly_to_daily_month_of_day

    subroutine monthly_to_daily_controls_1D(md, var_mon, var_ctrl)
        ! Mean-preserving control values for one column: var_ctrl = Ainv*var_mon.
        ! Feed var_ctrl (NOT var_mon) to interp_monthly_to_day to conserve.

        implicit none

        type(monthly_to_daily_class), intent(IN)  :: md
        real(wp), intent(IN)  :: var_mon(:)        ! (nmon)
        real(wp), intent(OUT) :: var_ctrl(:)       ! (nmon)

        var_ctrl = matmul(md%Ainv, var_mon)

        return

    end subroutine monthly_to_daily_controls_1D

    subroutine monthly_to_daily_controls_2D(md, var_mon, var_ctrl)
        ! Mean-preserving control values for a column list, layout (ncol,nmon).

        implicit none

        type(monthly_to_daily_class), intent(IN)  :: md
        real(wp), intent(IN)  :: var_mon(:,:)      ! (ncol,nmon)
        real(wp), intent(OUT) :: var_ctrl(:,:)     ! (ncol,nmon)

        ! var_ctrl(:,q) = sum_p var_mon(:,p) * Ainv(q,p)
        var_ctrl = matmul(var_mon, transpose(md%Ainv))

        return

    end subroutine monthly_to_daily_controls_2D

    subroutine interp_monthly_to_day_pt(md, var_ctrl, day, var_day)
        ! One column, one day.

        implicit none

        type(monthly_to_daily_class), intent(IN)  :: md
        real(wp), intent(IN)  :: var_ctrl(:)       ! (nmon)
        integer,  intent(IN)  :: day               ! (1..nday_year)
        real(wp), intent(OUT) :: var_day

        var_day = var_ctrl(md%m0(day))*md%wt0(day) + var_ctrl(md%m1(day))*md%wt1(day)

        return

    end subroutine interp_monthly_to_day_pt

    subroutine interp_monthly_to_day_1D(md, var_ctrl, day, var_day)
        ! Column list, one day. Layout (ncol,nmon) -> (ncol).

        implicit none

        type(monthly_to_daily_class), intent(IN)  :: md
        real(wp), intent(IN)  :: var_ctrl(:,:)     ! (ncol,nmon)
        integer,  intent(IN)  :: day               ! (1..nday_year)
        real(wp), intent(OUT) :: var_day(:)        ! (ncol)

        var_day = var_ctrl(:,md%m0(day))*md%wt0(day) + var_ctrl(:,md%m1(day))*md%wt1(day)

        return

    end subroutine interp_monthly_to_day_1D

    ! ==================================================================
    ! Internal
    ! ==================================================================

    subroutine interp_weights_monthly_to_daily(m0, m1, wt0, wt1, nday_mon, nmon)
        ! rembo's stencil (monthlydaily.f90). Days before the mid-month day
        ! blend the previous and current month; days after blend current and
        ! next; the mid-month day is the current month alone. Months wrap.

        implicit none

        integer,  intent(OUT) :: m0(:)
        integer,  intent(OUT) :: m1(:)
        real(wp), intent(OUT) :: wt0(:)
        real(wp), intent(OUT) :: wt1(:)
        integer,  intent(IN)  :: nday_mon
        integer,  intent(IN)  :: nmon

        ! Local variables
        integer :: dnow, q, k
        integer :: q0, q1, q2
        integer :: mid

        dnow = 0
        mid  = nday_mon / 2

        do q = 1, nmon
            do k = 1, nday_mon

                dnow = dnow + 1

                q1 = q                          ! current month
                q0 = q1 - 1                     ! previous
                if (q0 .eq. 0)    q0 = nmon
                q2 = q1 + 1                     ! next
                if (q2 .eq. nmon+1) q2 = 1

                if (k .lt. mid) then
                    m0(dnow)  = q0
                    m1(dnow)  = q1
                    wt0(dnow) = real(mid-k, wp) / real(nday_mon, wp)
                    wt1(dnow) = 1.0_wp - wt0(dnow)
                else if (k .gt. mid) then
                    m0(dnow)  = q1
                    m1(dnow)  = q2
                    wt0(dnow) = 1.0_wp - real(k-mid, wp) / real(nday_mon, wp)
                    wt1(dnow) = 1.0_wp - wt0(dnow)
                else
                    m0(dnow)  = q1
                    m1(dnow)  = q1
                    wt0(dnow) = 1.0_wp
                    wt1(dnow) = 0.0_wp
                end if

            end do
        end do

        return

    end subroutine interp_weights_monthly_to_daily

    function invert_matrix(A) result(Ainv)
        ! Dense inverse by Gauss-Jordan with partial pivoting. Used once at
        ! init on the small (nmon x nmon) mean operator, which is well
        ! conditioned (diagonally dominant), so this is ample.

        implicit none

        real(dp), intent(IN) :: A(:,:)
        real(dp) :: Ainv(size(A,1),size(A,2))

        ! Local variables
        real(dp), allocatable :: M(:,:)
        integer  :: n, i, j, p
        integer  :: ipiv
        real(dp) :: pmax, fac

        n = size(A,1)
        allocate(M(n,2*n))

        M = 0.0_dp
        M(:,1:n) = A
        do i = 1, n
            M(i,n+i) = 1.0_dp
        end do

        do p = 1, n
            ! Partial pivot
            ipiv = p
            pmax = abs(M(p,p))
            do i = p+1, n
                if (abs(M(i,p)) .gt. pmax) then
                    pmax = abs(M(i,p))
                    ipiv = i
                end if
            end do
            if (ipiv .ne. p) then
                call swap_rows(M, p, ipiv)
            end if

            ! Normalise pivot row
            fac = M(p,p)
            M(p,:) = M(p,:) / fac

            ! Eliminate other rows
            do i = 1, n
                if (i .ne. p) then
                    fac = M(i,p)
                    M(i,:) = M(i,:) - fac*M(p,:)
                end if
            end do
        end do

        Ainv = M(:,n+1:2*n)

        deallocate(M)

        return

    end function invert_matrix

    subroutine swap_rows(M, i, j)

        implicit none

        real(dp), intent(INOUT) :: M(:,:)
        integer,  intent(IN)    :: i, j

        ! Local variables
        real(dp), allocatable :: tmp(:)

        allocate(tmp(size(M,2)))
        tmp    = M(i,:)
        M(i,:) = M(j,:)
        M(j,:) = tmp
        deallocate(tmp)

        return

    end subroutine swap_rows

end module chion_forcing_monthly
