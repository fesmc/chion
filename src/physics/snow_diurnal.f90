module snow_diurnal
    ! Energy-conserving adaptive diurnal shortwave substepping: solar geometry,
    ! interval averages, and the substep-count decision.
    !
    ! Port of Chion.jl/src/processes/diurnal_shortwave.jl.
    !
    ! This module is PURE: it holds no state, touches no column arrays, and
    ! depends only on chion_defs. Everything is a function of latitude, solar
    ! longitude and a pair of hour angles.
    !
    ! PRECISION: the geometry is evaluated in wp_acc internally and returned in
    ! wp. Two of the expressions are differences of nearly equal numbers --
    ! sin(h_b) - sin(h_a) for a narrow interval, and the daylight clipping
    ! max/min against +-h0 -- and the acceptance test requires that a full
    ! [-pi,pi] tiling recovers the daily mean to 1e-6 relative, which is at the
    ! edge of sp. See docs/PLAN.md section 3.1.
    !
    ! NOTE FOR THE CALLER (WP8): the substep loop accumulates surface energy
    ! across substeps, and that accumulation MUST be real(wp_acc) -- it is the
    ! third of the three named dp-local expressions in docs/PLAN.md section 3.1.
    ! This module deliberately does not own the loop; it only supplies the
    ! bounds (diurnal_substep_bounds) and the interval averages.
    !
    ! PRESERVED QUIRKS:
    !   * Interval averages divide by the FULL interval width h_b - h_a, not by
    !     the daylight-clipped width. That is what makes a partly- or fully-
    !     nocturnal interval dilute correctly, and what makes the tiling sum
    !     back to the daily mean.
    !   * diurnal_substep_count returns ONLY 1 or max_substeps -- never an
    !     intermediate value -- gated by eight conditions.
    !   * The dt_days window [0.75, 1.25] is a bare pair of literals upstream;
    !     it means the whole scheme silently disables itself for any timestep
    !     that is not roughly one day (docs/PLAN.md section 5, item 3).

    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    use chion_defs, only : wp, wp_acc, io_unit_err

    implicit none

    private

    real(wp_acc), parameter :: PI_ACC  = 3.14159265358979323846_wp_acc
    real(wp_acc), parameter :: DEG2RAD = PI_ACC/180.0_wp_acc

    ! Obliquity of the ecliptic (Chion.jl diurnal_shortwave.jl:5).
    real(wp), parameter, public :: DIURNAL_OBLIQUITY_DEG = 23.439291_wp   ! [deg]

    ! Timestep window within which substepping is permitted
    ! (Chion.jl diurnal_shortwave.jl:111-112).
    real(wp), parameter, public :: DIURNAL_DT_DAYS_MIN = 0.75_wp   ! [d]
    real(wp), parameter, public :: DIURNAL_DT_DAYS_MAX = 1.25_wp   ! [d]

    public :: solar_declination_deg
    public :: sunset_hour_angle
    public :: diurnal_daylight_integral
    public :: diurnal_shortwave_interval_average
    public :: diurnal_shortwave_peak_flux
    public :: diurnal_temperature_interval_average
    public :: diurnal_substep_count
    public :: diurnal_substep_bounds

contains

    pure function solar_declination_deg(solar_longitude_deg) result(dec_deg)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:6-7:
        !     delta = asind(sind(obliquity)*sind(lambda))

        implicit none

        real(wp), intent(IN) :: solar_longitude_deg   ! [deg] lambda
        real(wp) :: dec_deg                           ! [deg] delta

        ! Local variables
        real(wp_acc) :: s

        s = sin(real(DIURNAL_OBLIQUITY_DEG,wp_acc)*DEG2RAD) &
            *sin(real(solar_longitude_deg,wp_acc)*DEG2RAD)

        ! Argument of asin is a product of two sines, so |s| <= 1 exactly in
        ! exact arithmetic; clamp anyway so round-off cannot raise invalid.
        s = min(max(s,-1.0_wp_acc),1.0_wp_acc)

        dec_deg = real(asin(s)/DEG2RAD,wp)

        return

    end function solar_declination_deg

    pure function sunset_hour_angle(latitude_deg,declination_deg) result(h0)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:9-17:
        !     cos_h0 = -tand(phi)*tand(delta)
        !     >= 1  -> 0   (polar night: the sun never rises)
        !     <= -1 -> pi  (polar day:   the sun never sets)
        !     else  -> acos(cos_h0)
        ! The explicit branches are what keep acos in range; do not replace
        ! them with a clamp, because the returned 0 also flags polar night to
        ! the callers below.

        implicit none

        real(wp), intent(IN) :: latitude_deg      ! [deg N] phi
        real(wp), intent(IN) :: declination_deg   ! [deg]   delta
        real(wp) :: h0                            ! [rad]

        ! Local variables
        real(wp_acc) :: cos_h0

        cos_h0 = -tan(real(latitude_deg,wp_acc)*DEG2RAD) &
                 *tan(real(declination_deg,wp_acc)*DEG2RAD)

        if (cos_h0 .ge. 1.0_wp_acc) then
            h0 = 0.0_wp
        else if (cos_h0 .le. -1.0_wp_acc) then
            h0 = real(PI_ACC,wp)
        else
            h0 = real(acos(cos_h0),wp)
        end if

        return

    end function sunset_hour_angle

    pure subroutine diurnal_daylight_integral(latitude_deg,solar_longitude_deg, &
                                              declination_deg,h0,I_day,A,B)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:19-35. The shared geometry
        ! terms:
        !     A     = sin(phi)*sin(delta)
        !     B     = cos(phi)*cos(delta)
        !     I_day = 2*(h0*A + B*sin(h0))
        ! I_day is the integral of the solar shape mu(h) = A + B*cos(h) over the
        ! daylight interval [-h0, h0].

        implicit none

        real(wp),     intent(IN)  :: latitude_deg          ! [deg N]
        real(wp),     intent(IN)  :: solar_longitude_deg   ! [deg]
        real(wp),     intent(OUT) :: declination_deg       ! [deg]
        real(wp),     intent(OUT) :: h0                    ! [rad] sunset hour angle
        real(wp_acc), intent(OUT) :: I_day                 ! [rad] daylight integral
        real(wp_acc), intent(OUT) :: A                     ! [1] sin(phi)*sin(delta)
        real(wp_acc), intent(OUT) :: B                     ! [1] cos(phi)*cos(delta)

        ! Local variables
        real(wp_acc) :: lat_rad, dec_rad, h0_acc

        declination_deg = solar_declination_deg(solar_longitude_deg)
        h0              = sunset_hour_angle(latitude_deg,declination_deg)

        lat_rad = real(latitude_deg,wp_acc)*DEG2RAD
        dec_rad = real(declination_deg,wp_acc)*DEG2RAD

        A = sin(lat_rad)*sin(dec_rad)
        B = cos(lat_rad)*cos(dec_rad)

        h0_acc = real(h0,wp_acc)
        I_day  = 2.0_wp_acc*(h0_acc*A + B*sin(h0_acc))

        return

    end subroutine diurnal_daylight_integral

    pure function diurnal_shortwave_interval_average(shortwave_daily_mean, &
                                                     latitude_deg,solar_longitude_deg, &
                                                     hour_angle_start,hour_angle_end) result(q_sw)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:37-66.
        !
        !     d_a  = max(h_a, -h0)
        !     d_b  = min(h_b,  h0)
        !     I_ab = (d_b - d_a)*A + B*(sin(d_b) - sin(d_a))
        !     S    = Qbar*2*pi/I_day
        !     ->     max(S*I_ab/(h_b - h_a), 0)
        !
        ! The divisor is the FULL interval width, so a fully nocturnal interval
        ! returns 0 and a partly nocturnal one is diluted. Summing w_i*q_i over
        ! a full [-pi,pi] tiling therefore returns 2*pi*Qbar.

        implicit none

        real(wp), intent(IN) :: shortwave_daily_mean   ! [W m-2] Qbar
        real(wp), intent(IN) :: latitude_deg           ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg    ! [deg]
        real(wp), intent(IN) :: hour_angle_start       ! [rad] h_a
        real(wp), intent(IN) :: hour_angle_end         ! [rad] h_b
        real(wp) :: q_sw                               ! [W m-2] interval mean

        ! Local variables
        real(wp)     :: dec_deg, h0
        real(wp_acc) :: width, I_day, A, B, d_a, d_b, I_ab, scale

        q_sw = 0.0_wp

        width = real(hour_angle_end,wp_acc) - real(hour_angle_start,wp_acc)

        if (shortwave_daily_mean .le. 0.0_wp)      return
        if (width .le. 0.0_wp_acc)                 return
        if (.not. ieee_is_finite(latitude_deg))        return
        if (.not. ieee_is_finite(solar_longitude_deg)) return

        call diurnal_daylight_integral(latitude_deg,solar_longitude_deg, &
                                       dec_deg,h0,I_day,A,B)

        ! Julia tests I_day against eps(Float64); the equivalent guard for the
        ! wp_acc arithmetic used here is epsilon(1.0_wp_acc). Polar night gives
        ! h0 = 0 exactly and is caught by the second test.
        if (I_day .le. epsilon(1.0_wp_acc)) return
        if (h0    .le. 0.0_wp)              return

        d_a = max(real(hour_angle_start,wp_acc),-real(h0,wp_acc))
        d_b = min(real(hour_angle_end,wp_acc),   real(h0,wp_acc))

        if (d_b .le. d_a) return

        I_ab  = (d_b - d_a)*A + B*(sin(d_b) - sin(d_a))
        scale = real(shortwave_daily_mean,wp_acc)*2.0_wp_acc*PI_ACC/I_day

        q_sw = real(max(scale*I_ab/width,0.0_wp_acc),wp)

        return

    end function diurnal_shortwave_interval_average

    pure function diurnal_shortwave_peak_flux(shortwave_daily_mean, &
                                              latitude_deg,solar_longitude_deg) result(q_peak)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:68-84.
        !     Q_peak = S*(A + B)      the reconstructed local-noon flux

        implicit none

        real(wp), intent(IN) :: shortwave_daily_mean   ! [W m-2]
        real(wp), intent(IN) :: latitude_deg           ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg    ! [deg]
        real(wp) :: q_peak                             ! [W m-2]

        ! Local variables
        real(wp)     :: dec_deg, h0
        real(wp_acc) :: I_day, A, B, scale

        q_peak = 0.0_wp

        if (shortwave_daily_mean .le. 0.0_wp)          return
        if (.not. ieee_is_finite(latitude_deg))        return
        if (.not. ieee_is_finite(solar_longitude_deg)) return

        call diurnal_daylight_integral(latitude_deg,solar_longitude_deg, &
                                       dec_deg,h0,I_day,A,B)

        if (I_day .le. epsilon(1.0_wp_acc)) return
        if (h0    .le. 0.0_wp)              return

        scale = real(shortwave_daily_mean,wp_acc)*2.0_wp_acc*PI_ACC/I_day

        q_peak = real(max(scale*(A + B),0.0_wp_acc),wp)

        return

    end function diurnal_shortwave_peak_flux

    pure function diurnal_temperature_interval_average(air_temperature_daily_mean,amplitude, &
                                                       hour_angle_start,hour_angle_end) result(t_air)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:86-98.
        !     T_ab = Tbar + A_T*(sin(h_b) - sin(h_a))/(h_b - h_a)
        ! i.e. the interval mean of T(h) = Tbar + A_T*cos(h), warmest at solar
        ! noon (h = 0) and preserving the daily mean over [-pi,pi].
        !
        ! A non-positive amplitude or a non-positive width returns the daily
        ! mean unchanged -- note this returns Tbar, NOT zero, unlike the
        ! shortwave routines above.

        implicit none

        real(wp), intent(IN) :: air_temperature_daily_mean   ! [K] Tbar
        real(wp), intent(IN) :: amplitude                    ! [K] A_T (half-amplitude)
        real(wp), intent(IN) :: hour_angle_start             ! [rad]
        real(wp), intent(IN) :: hour_angle_end               ! [rad]
        real(wp) :: t_air                                    ! [K]

        ! Local variables
        real(wp_acc) :: width, dsin

        t_air = air_temperature_daily_mean

        width = real(hour_angle_end,wp_acc) - real(hour_angle_start,wp_acc)

        if (amplitude .le. 0.0_wp)  return
        if (width .le. 0.0_wp_acc)  return

        ! Difference of nearly equal sines for a narrow interval -> wp_acc.
        dsin = sin(real(hour_angle_end,wp_acc)) - sin(real(hour_angle_start,wp_acc))

        t_air = real(real(air_temperature_daily_mean,wp_acc) &
                     + real(amplitude,wp_acc)*dsin/width, wp)

        return

    end function diurnal_temperature_interval_average

    pure function diurnal_substep_count(dt_days,shortwave_daily_mean,air_temperature, &
                                        min_air_temperature,latitude_deg,solar_longitude_deg, &
                                        threshold,max_substeps) result(n_substeps)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:100-127.
        !
        ! Returns ONLY 1 or max_substeps. The eight gating conditions, in order:
        !   1. max_substeps > 1
        !   2. dt_days >= 0.75
        !   3. dt_days <= 1.25
        !   4. shortwave_daily_mean > 0
        !   5. air_temperature > min_air_temperature   (STRICT)
        !   6. air_temperature, latitude_deg, solar_longitude_deg all finite
        !   7. peak - mean > threshold                 (STRICT; the Julia form
        !      is "<= threshold && return 1", so equality disables substepping)
        !   8. sunset hour angle > 0                   (not polar night)

        implicit none

        real(wp), intent(IN) :: dt_days                ! [d]
        real(wp), intent(IN) :: shortwave_daily_mean   ! [W m-2]
        real(wp), intent(IN) :: air_temperature        ! [K]
        real(wp), intent(IN) :: min_air_temperature    ! [K]
        real(wp), intent(IN) :: latitude_deg           ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg    ! [deg]
        real(wp), intent(IN) :: threshold              ! [W m-2] peak-minus-mean excess
        integer,  intent(IN) :: max_substeps           ! [1]
        integer :: n_substeps

        ! Local variables
        real(wp)     :: q_peak, dec_deg, h0
        real(wp_acc) :: I_day, A, B

        n_substeps = 1

        if (max_substeps .le. 1) return

        if (dt_days .lt. DIURNAL_DT_DAYS_MIN)          return
        if (dt_days .gt. DIURNAL_DT_DAYS_MAX)          return
        if (shortwave_daily_mean .le. 0.0_wp)          return
        if (air_temperature .le. min_air_temperature)  return
        if (.not. ieee_is_finite(air_temperature))     return
        if (.not. ieee_is_finite(latitude_deg))        return
        if (.not. ieee_is_finite(solar_longitude_deg)) return

        q_peak = diurnal_shortwave_peak_flux(shortwave_daily_mean,latitude_deg,solar_longitude_deg)

        if (max(q_peak - shortwave_daily_mean,0.0_wp) .le. threshold) return

        call diurnal_daylight_integral(latitude_deg,solar_longitude_deg, &
                                       dec_deg,h0,I_day,A,B)

        if (h0 .le. 0.0_wp) return

        n_substeps = max_substeps

        return

    end function diurnal_substep_count

    subroutine diurnal_substep_bounds(substep_index,n_substeps,hour_angle_start,hour_angle_end)
        ! Chion.jl/src/step.jl:131-141. The day is tiled uniformly in hour angle
        ! from -pi to pi. The last substep's upper bound is set to exactly +pi
        ! rather than accumulated, so the tiling closes without round-off gaps.

        implicit none

        integer,  intent(IN)  :: substep_index    ! 1 .. n_substeps
        integer,  intent(IN)  :: n_substeps
        real(wp), intent(OUT) :: hour_angle_start ! [rad]
        real(wp), intent(OUT) :: hour_angle_end   ! [rad]

        ! Local variables
        real(wp_acc) :: day_start, day_end, width

        if (n_substeps .lt. 1) then
            write(io_unit_err,*) "diurnal_substep_bounds:: Error: n_substeps must be positive."
            write(io_unit_err,*) "n_substeps = ", n_substeps
            stop "Program stopped."
        end if

        if (substep_index .lt. 1 .or. substep_index .gt. n_substeps) then
            write(io_unit_err,*) "diurnal_substep_bounds:: Error: substep_index out of range."
            write(io_unit_err,*) "substep_index, n_substeps = ", substep_index, n_substeps
            stop "Program stopped."
        end if

        day_start = -PI_ACC
        day_end   =  PI_ACC
        width     = (day_end - day_start)/real(n_substeps,wp_acc)

        hour_angle_start = real(day_start + real(substep_index-1,wp_acc)*width,wp)

        if (substep_index .eq. n_substeps) then
            hour_angle_end = real(day_end,wp)
        else
            hour_angle_end = real(day_start + real(substep_index,wp_acc)*width,wp)
        end if

        return

    end subroutine diurnal_substep_bounds

end module snow_diurnal
