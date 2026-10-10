module snow_diurnal
    ! Energy-conserving adaptive diurnal shortwave substepping: solar geometry,
    ! interval averages, and the substep-count decision.
    !
    ! Port of Chion.jl/src/processes/diurnal_shortwave.jl.
    !
    ! This module holds no state, touches no column arrays, and depends only
    ! on chion_defs. Everything is a function of latitude, solar longitude
    ! and a pair of hour angles.
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

    ! Solar constant of the daily-mean top-of-atmosphere shortwave
    ! (Chion.jl 03bb445 surface_fluxes.jl SOLAR_CONSTANT).
    real(wp_acc), parameter, public :: DIURNAL_SOLAR_CONSTANT = 1361.0_wp_acc   ! [W m-2]

    ! Timestep window within which substepping is permitted
    ! (Chion.jl diurnal_shortwave.jl:111-112).
    real(wp), parameter, public :: DIURNAL_DT_DAYS_MIN = 0.75_wp   ! [d]
    real(wp), parameter, public :: DIURNAL_DT_DAYS_MAX = 1.25_wp   ! [d]

    ! Solar geometry of one column and day: what the interval averages, the
    ! peak flux, the substep criterion and the cloud-proxy TOA take from the
    ! latitude and the solar longitude. A column-day evaluates it once
    ! (diurnal_geometry) and passes it to the *_geom forms, instead of every
    ! call -- three per day plus one per substep -- repeating the trigonometry
    ! of diurnal_daylight_integral. The latitude/solar-longitude forms are
    ! wrappers over the same code, so both give identical results.
    ! The latitude-independent part of it: the solar declination of the day
    ! and its sine, cosine and tangent, the same for every column. A host step
    ! evaluates it once (solar_declination) and builds each column's geometry
    ! from it (diurnal_geometry(latitude_deg,decl)), instead of every column
    ! repeating the asin and three trigonometric functions of the day.
    type solar_declination_class
        logical      :: defined = .FALSE.        ! solar longitude finite
        real(wp)     :: declination_deg = 0.0_wp ! [deg]
        real(wp_acc) :: sin_dec = 0.0_wp_acc     ! [1]
        real(wp_acc) :: cos_dec = 0.0_wp_acc     ! [1]
        real(wp_acc) :: tan_dec = 0.0_wp_acc     ! [1]
    end type solar_declination_class

    type diurnal_geometry_class
        logical      :: defined = .FALSE.        ! latitude and solar longitude finite
        real(wp)     :: declination_deg = 0.0_wp ! [deg]
        real(wp)     :: h0 = 0.0_wp              ! [rad] sunset hour angle
        real(wp_acc) :: I_day = 0.0_wp_acc       ! [rad] daylight integral
        real(wp_acc) :: A = 0.0_wp_acc           ! [1] sin(phi)*sin(delta)
        real(wp_acc) :: B = 0.0_wp_acc           ! [1] cos(phi)*cos(delta)
    end type diurnal_geometry_class

    ! The substep tiling of one column-day: the n_substeps + 1 hour-angle
    ! bounds (diurnal_substep_bounds) and their sines, with the sines of the
    ! sunset/sunrise hour angles +-h0 that clip the shortwave intervals. The
    ! interval averages need these sines only; a substepped day forms each
    ! once (diurnal_substep_tiling) instead of four per substep.
    integer, parameter :: DIURNAL_MAX_SUBSTEPS = 24
    type diurnal_tiling_class
        integer      :: n_substeps = 0
        real(wp)     :: bounds(0:DIURNAL_MAX_SUBSTEPS) = 0.0_wp      ! [rad]
        real(wp_acc) :: sin_bounds(0:DIURNAL_MAX_SUBSTEPS) = 0.0_wp_acc
        real(wp_acc) :: sin_h0 = 0.0_wp_acc       ! sin(h0)
        real(wp_acc) :: sin_neg_h0 = 0.0_wp_acc   ! sin(-h0)
    end type diurnal_tiling_class

    interface diurnal_geometry
        module procedure diurnal_geometry_lon
        module procedure diurnal_geometry_decl
    end interface diurnal_geometry

    interface daily_toa_shortwave
        module procedure daily_toa_shortwave_lat
        module procedure daily_toa_shortwave_geom
    end interface daily_toa_shortwave

    interface diurnal_shortwave_interval_average
        module procedure diurnal_shortwave_interval_average_lat
        module procedure diurnal_shortwave_interval_average_geom
    end interface diurnal_shortwave_interval_average

    interface diurnal_shortwave_peak_flux
        module procedure diurnal_shortwave_peak_flux_lat
        module procedure diurnal_shortwave_peak_flux_geom
    end interface diurnal_shortwave_peak_flux

    interface diurnal_substep_count
        module procedure diurnal_substep_count_lat
        module procedure diurnal_substep_count_geom
    end interface diurnal_substep_count

    public :: solar_declination_class
    public :: solar_declination
    public :: diurnal_geometry_class
    public :: diurnal_geometry
    public :: solar_declination_deg
    public :: sunset_hour_angle
    public :: diurnal_daylight_integral
    public :: daily_toa_shortwave
    public :: calendar_solar_longitude_deg
    public :: diurnal_shortwave_interval_average
    public :: diurnal_shortwave_peak_flux
    public :: diurnal_temperature_amplitude
    public :: diurnal_temperature_interval_average
    public :: diurnal_substep_count
    public :: diurnal_substep_bounds
    public :: diurnal_tiling_class
    public :: diurnal_substep_tiling
    public :: diurnal_shortwave_substep_average
    public :: diurnal_temperature_substep_average

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

        h0 = sunset_hour_angle_tan(latitude_deg,tan(real(declination_deg,wp_acc)*DEG2RAD))

        return

    end function sunset_hour_angle

    pure function sunset_hour_angle_tan(latitude_deg,tan_dec) result(h0)
        ! sunset_hour_angle from the tangent of the declination.

        implicit none

        real(wp),     intent(IN) :: latitude_deg  ! [deg N] phi
        real(wp_acc), intent(IN) :: tan_dec       ! [1] tan(delta)
        real(wp) :: h0                            ! [rad]

        ! Local variables
        real(wp_acc) :: cos_h0

        cos_h0 = -tan(real(latitude_deg,wp_acc)*DEG2RAD)*tan_dec

        if (cos_h0 .ge. 1.0_wp_acc) then
            h0 = 0.0_wp
        else if (cos_h0 .le. -1.0_wp_acc) then
            h0 = real(PI_ACC,wp)
        else
            h0 = real(acos(cos_h0),wp)
        end if

        return

    end function sunset_hour_angle_tan

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
        type(solar_declination_class) :: decl

        decl = declination_terms(solar_declination_deg(solar_longitude_deg))

        declination_deg = decl%declination_deg

        call daylight_integral_decl(latitude_deg,decl,h0,I_day,A,B)

        return

    end subroutine diurnal_daylight_integral

    pure function declination_terms(declination_deg) result(decl)
        ! The declination and its sine, cosine and tangent (in radians, from
        ! the wp declination as diurnal_daylight_integral always took it).

        implicit none

        real(wp), intent(IN) :: declination_deg   ! [deg]
        type(solar_declination_class) :: decl

        ! Local variables
        real(wp_acc) :: dec_rad

        decl%defined         = .TRUE.
        decl%declination_deg = declination_deg

        dec_rad      = real(declination_deg,wp_acc)*DEG2RAD
        decl%sin_dec = sin(dec_rad)
        decl%cos_dec = cos(dec_rad)
        decl%tan_dec = tan(dec_rad)

        return

    end function declination_terms

    pure subroutine daylight_integral_decl(latitude_deg,decl,h0,I_day,A,B)
        ! diurnal_daylight_integral from the day's declination terms.

        implicit none

        real(wp),                      intent(IN)  :: latitude_deg   ! [deg N]
        type(solar_declination_class), intent(IN)  :: decl
        real(wp),                      intent(OUT) :: h0             ! [rad] sunset hour angle
        real(wp_acc),                  intent(OUT) :: I_day          ! [rad] daylight integral
        real(wp_acc),                  intent(OUT) :: A              ! [1] sin(phi)*sin(delta)
        real(wp_acc),                  intent(OUT) :: B              ! [1] cos(phi)*cos(delta)

        ! Local variables
        real(wp_acc) :: lat_rad, h0_acc

        h0 = sunset_hour_angle_tan(latitude_deg,decl%tan_dec)

        lat_rad = real(latitude_deg,wp_acc)*DEG2RAD

        A = sin(lat_rad)*decl%sin_dec
        B = cos(lat_rad)*decl%cos_dec

        h0_acc = real(h0,wp_acc)
        I_day  = 2.0_wp_acc*(h0_acc*A + B*sin(h0_acc))

        return

    end subroutine daylight_integral_decl

    pure function solar_declination(solar_longitude_deg) result(decl)
        ! The day's declination terms (solar_declination_class). Undefined
        ! (all zero) when the solar longitude is not finite.

        implicit none

        real(wp), intent(IN) :: solar_longitude_deg   ! [deg]
        type(solar_declination_class) :: decl

        if (.not. ieee_is_finite(solar_longitude_deg)) return

        decl = declination_terms(solar_declination_deg(solar_longitude_deg))

        return

    end function solar_declination

    pure function diurnal_geometry_lon(latitude_deg,solar_longitude_deg) result(geom)
        ! The column-day's solar geometry (diurnal_geometry_class). Undefined
        ! (all zero) when the latitude or the solar longitude is not finite,
        ! which every *_geom consumer treats as "no daylight geometry", as the
        ! latitude/solar-longitude forms do with their finiteness tests.

        implicit none

        real(wp), intent(IN) :: latitude_deg          ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg   ! [deg]
        type(diurnal_geometry_class) :: geom

        geom = diurnal_geometry_decl(latitude_deg,solar_declination(solar_longitude_deg))

        return

    end function diurnal_geometry_lon

    pure function diurnal_geometry_decl(latitude_deg,decl) result(geom)
        ! diurnal_geometry from the day's declination terms (solar_declination),
        ! identical to the latitude/solar-longitude form.

        implicit none

        real(wp),                      intent(IN) :: latitude_deg   ! [deg N]
        type(solar_declination_class), intent(IN) :: decl
        type(diurnal_geometry_class) :: geom

        if (.not. ieee_is_finite(latitude_deg)) return
        if (.not. decl%defined)                 return

        geom%defined         = .TRUE.
        geom%declination_deg = decl%declination_deg

        call daylight_integral_decl(latitude_deg,decl,geom%h0,geom%I_day,geom%A,geom%B)

        return

    end function diurnal_geometry_decl

    pure function daily_toa_shortwave_lat(latitude_deg,solar_longitude_deg,day_of_year) result(toa)
        ! daily_toa_shortwave from latitude and solar longitude.

        implicit none

        real(wp), intent(IN) :: latitude_deg          ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg   ! [deg]
        real(wp), intent(IN) :: day_of_year           ! [d] fractional, 1-based
        real(wp) :: toa                               ! [W m-2]

        toa = daily_toa_shortwave_geom(diurnal_geometry(latitude_deg,solar_longitude_deg), &
                                       day_of_year)

        return

    end function daily_toa_shortwave_lat

    pure function daily_toa_shortwave_geom(geom,day_of_year) result(toa)
        ! Chion.jl/src/processes/surface_fluxes.jl:_daily_toa_shortwave
        ! (03bb445), the daily-mean top-of-atmosphere shortwave of a fixed,
        ! modern orbit:
        !     TOA = S0*(1 + 0.033*cos(2*pi*doy/365))*max(I_day,0)/(2*pi)
        ! with I_day the daylight integral of diurnal_daylight_integral. Only
        ! the cloud-proxy longwave uses it, and only when the host supplies no
        ! TOA of its own (docs/porting_notes.md D33).
        !
        ! An undefined geometry (non-finite latitude or solar longitude) has
        ! I_day = 0 and gives 0; the cloud proxy does not ask for it then.

        implicit none

        type(diurnal_geometry_class), intent(IN) :: geom
        real(wp),                     intent(IN) :: day_of_year   ! [d] fractional, 1-based
        real(wp) :: toa                                           ! [W m-2]

        ! Local variables
        real(wp_acc) :: eccentricity

        eccentricity = 1.0_wp_acc + 0.033_wp_acc &
                       *cos(2.0_wp_acc*PI_ACC*real(day_of_year,wp_acc)/365.0_wp_acc)

        toa = real(DIURNAL_SOLAR_CONSTANT*eccentricity*max(geom%I_day,0.0_wp_acc) &
                   /(2.0_wp_acc*PI_ACC),wp)

        return

    end function daily_toa_shortwave_geom

    pure function calendar_solar_longitude_deg(day_of_year) result(lon)
        ! Chion.jl/src/forcing.jl:_solar_longitude_deg_from_calendar_day: the
        ! solar longitude of a calendar day of year, from the low-precision
        ! solar-position series (mean longitude, mean anomaly, equation of
        ! centre), with day 1 = 1 January:
        !     d = doy - 1
        !     L = 280.46646 + 0.98564736 d,   g = 357.52911 + 0.98560028 d
        !     lambda = mod(L + 1.914602 sin(g) + 0.019993 sin(2g), 360)
        ! A host or driver that has no orbital solar longitude of its own
        ! derives it here, so that chion and Chion.jl see the same season for
        ! the same calendar day (the diurnal substeps, the SEMIX coszm, the
        ! cloud-proxy TOA).

        implicit none

        real(wp), intent(IN) :: day_of_year   ! [d] fractional, 1-based
        real(wp) :: lon                       ! [deg]

        ! Local variables
        real(wp_acc) :: d, mean_longitude, mean_anomaly, true_longitude

        d              = real(day_of_year,wp_acc) - 1.0_wp_acc
        mean_longitude = 280.46646_wp_acc + 0.98564736_wp_acc*d
        mean_anomaly   = 357.52911_wp_acc + 0.98560028_wp_acc*d
        true_longitude = mean_longitude &
                         + 1.914602_wp_acc*sin(mean_anomaly*DEG2RAD) &
                         + 0.019993_wp_acc*sin(2.0_wp_acc*mean_anomaly*DEG2RAD)

        lon = real(modulo(true_longitude,360.0_wp_acc),wp)

        return

    end function calendar_solar_longitude_deg

    pure function diurnal_shortwave_interval_average_lat(shortwave_daily_mean, &
                                                         latitude_deg,solar_longitude_deg, &
                                                         hour_angle_start,hour_angle_end) result(q_sw)
        ! diurnal_shortwave_interval_average from latitude and solar longitude.

        implicit none

        real(wp), intent(IN) :: shortwave_daily_mean   ! [W m-2] Qbar
        real(wp), intent(IN) :: latitude_deg           ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg    ! [deg]
        real(wp), intent(IN) :: hour_angle_start       ! [rad] h_a
        real(wp), intent(IN) :: hour_angle_end         ! [rad] h_b
        real(wp) :: q_sw                               ! [W m-2] interval mean

        q_sw = diurnal_shortwave_interval_average_geom(shortwave_daily_mean, &
                   diurnal_geometry(latitude_deg,solar_longitude_deg), &
                   hour_angle_start,hour_angle_end)

        return

    end function diurnal_shortwave_interval_average_lat

    pure function diurnal_shortwave_interval_average_geom(shortwave_daily_mean,geom, &
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

        real(wp),                     intent(IN) :: shortwave_daily_mean   ! [W m-2] Qbar
        type(diurnal_geometry_class), intent(IN) :: geom
        real(wp),                     intent(IN) :: hour_angle_start       ! [rad] h_a
        real(wp),                     intent(IN) :: hour_angle_end         ! [rad] h_b
        real(wp) :: q_sw                                                   ! [W m-2] interval mean

        q_sw = shortwave_interval_average_sines(shortwave_daily_mean,geom, &
                   hour_angle_start,hour_angle_end, &
                   sin(real(hour_angle_start,wp_acc)),sin(real(hour_angle_end,wp_acc)), &
                   sin(real(geom%h0,wp_acc)),sin(-real(geom%h0,wp_acc)))

        return

    end function diurnal_shortwave_interval_average_geom

    pure function shortwave_interval_average_sines(shortwave_daily_mean,geom, &
                                                   hour_angle_start,hour_angle_end, &
                                                   sin_start,sin_end,sin_h0,sin_neg_h0) result(q_sw)
        ! diurnal_shortwave_interval_average_geom given the sines of the bounds
        ! and of +-h0: sin(d_a) and sin(d_b) are those of the clipped bounds.

        implicit none

        real(wp),                     intent(IN) :: shortwave_daily_mean   ! [W m-2] Qbar
        type(diurnal_geometry_class), intent(IN) :: geom
        real(wp),                     intent(IN) :: hour_angle_start       ! [rad] h_a
        real(wp),                     intent(IN) :: hour_angle_end         ! [rad] h_b
        real(wp_acc),                 intent(IN) :: sin_start              ! sin(h_a)
        real(wp_acc),                 intent(IN) :: sin_end                ! sin(h_b)
        real(wp_acc),                 intent(IN) :: sin_h0                 ! sin(h0)
        real(wp_acc),                 intent(IN) :: sin_neg_h0             ! sin(-h0)
        real(wp) :: q_sw                                                   ! [W m-2] interval mean

        ! Local variables
        real(wp_acc) :: width, d_a, d_b, I_ab, scale, sin_a, sin_b

        q_sw = 0.0_wp

        width = real(hour_angle_end,wp_acc) - real(hour_angle_start,wp_acc)

        if (shortwave_daily_mean .le. 0.0_wp) return
        if (width .le. 0.0_wp_acc)            return
        if (.not. geom%defined)               return

        ! Julia tests I_day against eps(Float64); the equivalent guard for the
        ! wp_acc arithmetic used here is epsilon(1.0_wp_acc). Polar night gives
        ! h0 = 0 exactly and is caught by the second test.
        if (geom%I_day .le. epsilon(1.0_wp_acc)) return
        if (geom%h0    .le. 0.0_wp)              return

        d_a = max(real(hour_angle_start,wp_acc),-real(geom%h0,wp_acc))
        d_b = min(real(hour_angle_end,wp_acc),   real(geom%h0,wp_acc))

        if (d_b .le. d_a) return

        sin_a = sin_start
        if (real(hour_angle_start,wp_acc) .lt. -real(geom%h0,wp_acc)) sin_a = sin_neg_h0
        sin_b = sin_end
        if (real(hour_angle_end,wp_acc) .gt. real(geom%h0,wp_acc)) sin_b = sin_h0

        I_ab  = (d_b - d_a)*geom%A + geom%B*(sin_b - sin_a)
        scale = real(shortwave_daily_mean,wp_acc)*2.0_wp_acc*PI_ACC/geom%I_day

        q_sw = real(max(scale*I_ab/width,0.0_wp_acc),wp)

        return

    end function shortwave_interval_average_sines

    pure function diurnal_shortwave_peak_flux_lat(shortwave_daily_mean, &
                                                  latitude_deg,solar_longitude_deg) result(q_peak)
        ! diurnal_shortwave_peak_flux from latitude and solar longitude.

        implicit none

        real(wp), intent(IN) :: shortwave_daily_mean   ! [W m-2]
        real(wp), intent(IN) :: latitude_deg           ! [deg N]
        real(wp), intent(IN) :: solar_longitude_deg    ! [deg]
        real(wp) :: q_peak                             ! [W m-2]

        q_peak = diurnal_shortwave_peak_flux_geom(shortwave_daily_mean, &
                     diurnal_geometry(latitude_deg,solar_longitude_deg))

        return

    end function diurnal_shortwave_peak_flux_lat

    pure function diurnal_shortwave_peak_flux_geom(shortwave_daily_mean,geom) result(q_peak)
        ! Chion.jl/src/processes/diurnal_shortwave.jl:68-84.
        !     Q_peak = S*(A + B)      the reconstructed local-noon flux

        implicit none

        real(wp),                     intent(IN) :: shortwave_daily_mean   ! [W m-2]
        type(diurnal_geometry_class), intent(IN) :: geom
        real(wp) :: q_peak                                                 ! [W m-2]

        ! Local variables
        real(wp_acc) :: scale

        q_peak = 0.0_wp

        if (shortwave_daily_mean .le. 0.0_wp) return
        if (.not. geom%defined)               return

        if (geom%I_day .le. epsilon(1.0_wp_acc)) return
        if (geom%h0    .le. 0.0_wp)              return

        scale = real(shortwave_daily_mean,wp_acc)*2.0_wp_acc*PI_ACC/geom%I_day

        q_peak = real(max(scale*(geom%A + geom%B),0.0_wp_acc),wp)

        return

    end function diurnal_shortwave_peak_flux_geom

    pure function diurnal_temperature_amplitude(amplitude_base,gradient_per_km, &
                                                reference_height,amplitude_max, &
                                                surface_height) result(amplitude)
        ! Chion.jl/src/step.jl:_diurnal_air_temperature (d0146e1, 9ec6cc7):
        !     A = clamp(A0 + gamma*max(z_s - z_ref, 0), 0, A_max)
        ! with gamma converted K/km -> K/m as in the BESSIModel constructor
        ! (src/models.jl: gradient_c_per_km/1000).
        !
        ! A non-finite surface height (a host that does not supply one) has no
        ! height excess, so the amplitude stays A0 clamped -- never NaN
        ! (Chion.jl 03bb445; our 1c.6). The default parameters
        ! (gamma = 0, A_max large) give A = A0 exactly.

        implicit none

        real(wp), intent(IN) :: amplitude_base     ! [K] A0
        real(wp), intent(IN) :: gradient_per_km    ! [K km-1] gamma
        real(wp), intent(IN) :: reference_height   ! [m] z_ref
        real(wp), intent(IN) :: amplitude_max      ! [K] A_max
        real(wp), intent(IN) :: surface_height     ! [m] z_s
        real(wp) :: amplitude                      ! [K]

        ! Local variables
        real(wp) :: height_excess

        if (ieee_is_finite(surface_height)) then
            height_excess = max(surface_height - reference_height, 0.0_wp)
        else
            height_excess = 0.0_wp
        end if

        amplitude = min(max(amplitude_base + (gradient_per_km/1000.0_wp)*height_excess, &
                            0.0_wp), amplitude_max)

        return

    end function diurnal_temperature_amplitude

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

        t_air = temperature_interval_average_sines(air_temperature_daily_mean,amplitude, &
                    hour_angle_start,hour_angle_end, &
                    sin(real(hour_angle_start,wp_acc)),sin(real(hour_angle_end,wp_acc)))

        return

    end function diurnal_temperature_interval_average

    pure function temperature_interval_average_sines(air_temperature_daily_mean,amplitude, &
                                                     hour_angle_start,hour_angle_end, &
                                                     sin_start,sin_end) result(t_air)
        ! diurnal_temperature_interval_average given the sines of the bounds.

        implicit none

        real(wp),     intent(IN) :: air_temperature_daily_mean   ! [K] Tbar
        real(wp),     intent(IN) :: amplitude                    ! [K] A_T (half-amplitude)
        real(wp),     intent(IN) :: hour_angle_start             ! [rad]
        real(wp),     intent(IN) :: hour_angle_end               ! [rad]
        real(wp_acc), intent(IN) :: sin_start                    ! sin(h_a)
        real(wp_acc), intent(IN) :: sin_end                      ! sin(h_b)
        real(wp) :: t_air                                        ! [K]

        ! Local variables
        real(wp_acc) :: width, dsin

        t_air = air_temperature_daily_mean

        width = real(hour_angle_end,wp_acc) - real(hour_angle_start,wp_acc)

        if (amplitude .le. 0.0_wp)  return
        if (width .le. 0.0_wp_acc)  return

        ! Difference of nearly equal sines for a narrow interval -> wp_acc.
        dsin = sin_end - sin_start

        t_air = real(real(air_temperature_daily_mean,wp_acc) &
                     + real(amplitude,wp_acc)*dsin/width, wp)

        return

    end function temperature_interval_average_sines

    function diurnal_substep_tiling(n_substeps,geom) result(tiling)
        ! The day's substep tiling (diurnal_tiling_class): bounds from
        ! diurnal_substep_bounds, so the end of substep k is bit for bit the
        ! start of k+1, and each sine formed once.

        implicit none

        integer,                      intent(IN) :: n_substeps
        type(diurnal_geometry_class), intent(IN) :: geom
        type(diurnal_tiling_class) :: tiling

        ! Local variables
        integer  :: k
        real(wp) :: hour_angle_start, hour_angle_end

        if (n_substeps .gt. DIURNAL_MAX_SUBSTEPS) then
            write(io_unit_err,*) "diurnal_substep_tiling:: Error: n_substeps exceeds ", &
                                 DIURNAL_MAX_SUBSTEPS
            write(io_unit_err,*) "n_substeps = ", n_substeps
            stop "Program stopped."
        end if

        tiling%n_substeps = n_substeps

        do k = 1, n_substeps
            call diurnal_substep_bounds(k,n_substeps,hour_angle_start,hour_angle_end)
            if (k .eq. 1) tiling%bounds(0) = hour_angle_start
            tiling%bounds(k) = hour_angle_end
        end do

        do k = 0, n_substeps
            tiling%sin_bounds(k) = sin(real(tiling%bounds(k),wp_acc))
        end do

        tiling%sin_h0     = sin(real(geom%h0,wp_acc))
        tiling%sin_neg_h0 = sin(-real(geom%h0,wp_acc))

        return

    end function diurnal_substep_tiling

    pure function diurnal_shortwave_substep_average(shortwave_daily_mean,geom,tiling,k) result(q_sw)
        ! diurnal_shortwave_interval_average over substep k of the tiling.

        implicit none

        real(wp),                     intent(IN) :: shortwave_daily_mean   ! [W m-2]
        type(diurnal_geometry_class), intent(IN) :: geom
        type(diurnal_tiling_class),   intent(IN) :: tiling
        integer,                      intent(IN) :: k
        real(wp) :: q_sw                                                   ! [W m-2]

        q_sw = shortwave_interval_average_sines(shortwave_daily_mean,geom, &
                   tiling%bounds(k-1),tiling%bounds(k), &
                   tiling%sin_bounds(k-1),tiling%sin_bounds(k), &
                   tiling%sin_h0,tiling%sin_neg_h0)

        return

    end function diurnal_shortwave_substep_average

    pure function diurnal_temperature_substep_average(air_temperature_daily_mean,amplitude, &
                                                      tiling,k) result(t_air)
        ! diurnal_temperature_interval_average over substep k of the tiling.

        implicit none

        real(wp),                   intent(IN) :: air_temperature_daily_mean   ! [K]
        real(wp),                   intent(IN) :: amplitude                    ! [K]
        type(diurnal_tiling_class), intent(IN) :: tiling
        integer,                    intent(IN) :: k
        real(wp) :: t_air                                                      ! [K]

        t_air = temperature_interval_average_sines(air_temperature_daily_mean,amplitude, &
                    tiling%bounds(k-1),tiling%bounds(k), &
                    tiling%sin_bounds(k-1),tiling%sin_bounds(k))

        return

    end function diurnal_temperature_substep_average

    pure function diurnal_substep_count_lat(dt_days,shortwave_daily_mean,air_temperature, &
                                            min_air_temperature,latitude_deg,solar_longitude_deg, &
                                            threshold,max_substeps) result(n_substeps)
        ! diurnal_substep_count from latitude and solar longitude.

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

        n_substeps = diurnal_substep_count_geom(dt_days,shortwave_daily_mean,air_temperature, &
                         min_air_temperature,diurnal_geometry(latitude_deg,solar_longitude_deg), &
                         threshold,max_substeps)

        return

    end function diurnal_substep_count_lat

    pure function diurnal_substep_count_geom(dt_days,shortwave_daily_mean,air_temperature, &
                                             min_air_temperature,geom, &
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

        real(wp),                     intent(IN) :: dt_days                ! [d]
        real(wp),                     intent(IN) :: shortwave_daily_mean   ! [W m-2]
        real(wp),                     intent(IN) :: air_temperature        ! [K]
        real(wp),                     intent(IN) :: min_air_temperature    ! [K]
        type(diurnal_geometry_class), intent(IN) :: geom
        real(wp),                     intent(IN) :: threshold              ! [W m-2] peak-minus-mean excess
        integer,                      intent(IN) :: max_substeps           ! [1]
        integer :: n_substeps

        ! Local variables
        real(wp) :: q_peak

        n_substeps = 1

        if (max_substeps .le. 1) return

        if (dt_days .lt. DIURNAL_DT_DAYS_MIN)          return
        if (dt_days .gt. DIURNAL_DT_DAYS_MAX)          return
        if (shortwave_daily_mean .le. 0.0_wp)          return
        if (air_temperature .le. min_air_temperature)  return
        if (.not. ieee_is_finite(air_temperature))     return
        if (.not. geom%defined)                        return

        q_peak = diurnal_shortwave_peak_flux_geom(shortwave_daily_mean,geom)

        if (max(q_peak - shortwave_daily_mean,0.0_wp) .le. threshold) return

        if (geom%h0 .le. 0.0_wp) return

        n_substeps = max_substeps

        return

    end function diurnal_substep_count_geom

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
