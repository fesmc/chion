module snow_albedo_semix
    ! CLIMBER-X SEMIX spectral snow albedo, collapsed to a broadband value.
    !
    ! Port of CLIMBER-X src/smb/smb_surface_par.f90 (surface_albedo /
    ! snow_albedo_ww; Willeit, Calov, Ganopolski). SEMIX carries the snow
    ! albedo in four bands -- {visible, near-IR} x {direct, diffuse} -- and
    ! builds net shortwave by weighting each band with its share of the
    ! incoming surface flux (downscaling.f90:212). chion carries a single
    ! broadband albedo through the energy balance, so this module computes the
    ! four bands and collapses them with the SAME weights the coupled model
    ! uses:
    !
    !   f_vis_dir = (1-cloud)*frac_vu     f_nir_dir = (1-cloud)*(1-frac_vu)
    !   f_vis_dif =    cloud *frac_vu     f_nir_dif =    cloud *(1-frac_vu)
    !   alpha_bb  = sum(f_band * alpha_band)
    !
    ! so (1-alpha_bb)*SW_down reproduces SEMIX's swnet under the assumption that
    ! the broadband SW_down partitions into the bands by those same weights.
    ! When spectral SW-down inputs become available (coupled runs) the bands can
    ! instead be applied to the spectral fluxes directly; at daily steps the
    ! broadband collapse is a faithful approximation. See docs/semix_port_scope.md.
    !
    ! Two band parameterizations are available, selected by
    ! c%semix_snow_albedo: Warren & Wiscombe 1980 (the CLIMBER-2 form) and
    ! Dang et al. 2015, which is what CLIMBER-X itself defaults to. Both take
    ! the diagnosed grain size and the dust concentration; grain aging and
    ! dust-in-snow darkening are shared by the two.

    use chion_defs, only : wp, chion_const_class, SEMIX_SNOW_ALBEDO_DANG

    implicit none

    private

    ! Upper bound on dust concentration in snow (smb_surface_par.f90:191).
    real(wp), parameter, public :: SEMIX_DUST_CON_MAX = 1000.0e-6_wp   ! [kg kg-1]

    public :: semix_surface_albedo
    public :: semix_snow_grain_size
    public :: semix_dust_concentration
    public :: semix_snow_albedo_bands
    public :: semix_bands_ww
    public :: semix_bands_dang
    public :: semix_broadband_albedo
    public :: semix_daily_coszm

contains

    subroutine semix_surface_albedo(c, t_skin, snowfall_rate, coszm, cloud, z_sur_std, &
                                    dust_con, albedo)
        ! Broadband SEMIX snow albedo for one column. Grain size is diagnosed
        ! from the skin temperature and the current snowfall rate; the caller
        ! resolves coszm (measured or the lat/solar-longitude fallback below),
        ! the cloud fraction and the subgrid orography, and supplies the dust
        ! concentration.
        !
        ! Deliberately takes plain scalars rather than chion's forcing type, as
        ! snow_albedo.f90 does: this module is pure physics and stays decoupled
        ! from the forcing contract.

        implicit none

        type(chion_const_class), intent(IN)  :: c
        real(wp),                intent(IN)  :: t_skin         ! [K] surface temperature
        real(wp),                intent(IN)  :: snowfall_rate  ! [kg m-2 s-1]
        real(wp),                intent(IN)  :: coszm          ! [1] daily-mean cos(zenith)
        real(wp),                intent(IN)  :: cloud          ! [1] cloud fraction
        real(wp),                intent(IN)  :: z_sur_std      ! [m] subgrid height std dev
        real(wp),                intent(IN)  :: dust_con       ! [kg kg-1]
        real(wp),                intent(OUT) :: albedo         ! [1] broadband

        real(wp) :: snow_grain
        real(wp) :: av_dir, an_dir, av_dif, an_dif

        snow_grain = semix_snow_grain_size(t_skin, snowfall_rate, c)

        call semix_snow_albedo_bands(snow_grain, dust_con, coszm, z_sur_std, c, &
                                     av_dir, an_dir, av_dif, an_dif)

        albedo = semix_broadband_albedo(av_dir, an_dir, av_dif, an_dif, cloud, c%frac_vu)

        return

    end subroutine semix_surface_albedo

    pure function semix_dust_concentration(dust_dep, snowfall_rate, w_snow, w_snow_max, c) &
                  result(dust_con)
        ! Dust concentration in the surface snow (smb_surface_par.f90:179-218).
        !
        ! The base concentration is the ratio of the dust deposition rate to the
        ! snowfall rate (kg dust per kg snow). Meltwater scavenges only 10-30 %
        ! of the dust (Doherty 2013), so as the pack melts down from its
        ! seasonal peak the remaining dust concentrates: the amplification is
        ! driven by the drawdown (w_snow_max - w_snow), with w_snow_dust the SWE
        ! whose melt doubles the concentration, capped at five-fold.
        !
        ! Only the drawdown enters, never the absolute column mass, so this
        ! carries over unchanged from SEMIX's capped bulk snow reservoir to
        ! chion's much deeper firn column.

        implicit none

        real(wp),                intent(IN) :: dust_dep       ! [kg m-2 s-1]
        real(wp),                intent(IN) :: snowfall_rate  ! [kg m-2 s-1]
        real(wp),                intent(IN) :: w_snow         ! [kg m-2] column SWE
        real(wp),                intent(IN) :: w_snow_max     ! [kg m-2] seasonal peak SWE
        type(chion_const_class), intent(IN) :: c
        real(wp) :: dust_con                                  ! [kg kg-1]

        real(wp) :: melt_fac

        dust_con = dust_dep/max(1.0e-7_wp, snowfall_rate)

        if (w_snow .gt. 1.0_wp) then
            melt_fac = 1.0_wp + (w_snow_max - w_snow)/c%w_snow_dust
            melt_fac = max(melt_fac, 1.0_wp)
            melt_fac = min(melt_fac, 5.0_wp)
        else
            melt_fac = 1.0_wp
        end if

        dust_con = melt_fac*dust_con

        ! Scaling that stands in for different dust imaginary refractive indices.
        dust_con = c%dust_con_scale*dust_con

        if (dust_con .lt. 1.0e-15_wp) dust_con = 0.0_wp
        dust_con = min(dust_con, SEMIX_DUST_CON_MAX)

        return

    end function semix_dust_concentration

    pure function semix_snow_grain_size(t_skin, snowfall_rate, c) result(snow_grain)
        ! Diagnostic snow grain size (aging proxy), CLIMBER-2 form tuned to
        ! MARv3.6/CROCUS over Greenland (smb_surface_par.f90:143-172). Grain
        ! grows from fresh toward old as the surface warms and snowfall thins:
        ! it is a function of the current skin temperature and snowfall rate
        ! only, carrying no memory. snow_0 is stored per day (as in CLIMBER-X's
        ! namelist), so the snowfall rate is converted to per-day to match.

        implicit none

        real(wp),                intent(IN) :: t_skin         ! [K]
        real(wp),                intent(IN) :: snowfall_rate  ! [kg m-2 s-1]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: snow_grain                                ! [um]

        real(wp) :: f_tage1, f_tage2, f_tage, f_p, f_age, snow_per_day

        f_tage1 = exp(c%f_age_t*min(0.0_wp, t_skin - (c%T0 - c%dT_age)))
        f_tage2 = exp(          min(0.0_wp, t_skin - (c%T0 - c%dT_age)))
        f_tage  = f_tage1 + f_tage2

        snow_per_day = snowfall_rate*c%seconds_per_day
        f_p   = f_tage*(c%snow_0/max(1.0e-20_wp, snow_per_day))**c%snow_1
        f_age = 1.0_wp - log(1.0_wp + f_p)/f_p

        snow_grain = c%snow_grain_fresh + (c%snow_grain_old - c%snow_grain_fresh)*f_age

        return

    end function semix_snow_grain_size

    subroutine semix_snow_albedo_bands(snow_grain, dust_con, coszm, z_sur_std, c, &
                                       alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
        ! Four-band snow albedo, dispatching on the configured parameterization:
        ! Warren & Wiscombe 1980 (CLIMBER-2 form) or Dang et al. 2015. CLIMBER-X
        ! defaults to Dang (isnow_albedo = 2).

        implicit none

        real(wp),                intent(IN)  :: snow_grain, dust_con, coszm, z_sur_std
        type(chion_const_class), intent(IN)  :: c
        real(wp),                intent(OUT) :: alb_vis_dir, alb_nir_dir
        real(wp),                intent(OUT) :: alb_vis_dif, alb_nir_dif

        if (c%semix_snow_albedo .eq. SEMIX_SNOW_ALBEDO_DANG) then
            call semix_bands_dang(snow_grain, dust_con, coszm, z_sur_std, c, &
                                        alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
        else
            call semix_bands_ww(snow_grain, dust_con, coszm, c, &
                                      alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
        end if

        return

    end subroutine semix_snow_albedo_bands

    subroutine semix_bands_dang(snow_grain, dust_con, coszm, z_sur_std, c, &
                                      alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
        ! Four-band snow albedo, Dang et al. 2015 (smb_surface_par.f90:288-389).
        ! Albedo is a quadratic in log10(grain radius / 100 um); dust enters as a
        ! black-carbon-equivalent concentration (their eq. 9) and darkens only the
        ! visible bands. The direct bands use a zenith-corrected effective grain
        ! size. A tanh orographic term reduces all four bands over rough terrain
        ! (off by default: k_sigma_orog = 0).

        implicit none

        real(wp),                intent(IN)  :: snow_grain, dust_con, coszm, z_sur_std
        type(chion_const_class), intent(IN)  :: c
        real(wp),                intent(OUT) :: alb_vis_dir, alb_nir_dir
        real(wp),                intent(OUT) :: alb_vis_dif, alb_nir_dif

        real(wp), parameter :: r0       = 100.0_wp    ! [um]
        real(wp), parameter :: c0       = 1.0e-6_wp   ! [kg kg-1]
        real(wp), parameter :: dust_min = 1.0e-8_wp   ! [kg kg-1]

        real(wp) :: r, rn, cc, x, f, h, p, dalb_orog
        real(wp) :: dalpha_vis_dif, dalpha_vis_dir

        x = 0.0_wp
        if (dust_con .gt. dust_min) x = log10(dust_con*1.0e6_wp)

        dalb_orog = c%k_sigma_orog*tanh(z_sur_std/c%sigma_orog_crit)

        ! --- diffuse visible: aging, then black-carbon-equivalent darkening
        r  = snow_grain
        rn = log10(r/r0)
        alb_vis_dif = 0.9856_wp + c%dalb_snow_vis - 0.0202_wp*rn - 0.0125_wp*rn**2
        if (dust_con .gt. dust_min) then
            f  = 152.0_wp + 15.92_wp*x - 0.39_wp*x**2
            cc = dust_con/f
            h  = cc/c0*(r/r0)**0.73_wp
            p  = log10(h)
            dalpha_vis_dif = 10.0_wp**(-0.050_wp*p**2 + 0.514_wp*p - 0.890_wp)
        else
            dalpha_vis_dif = 0.0_wp
        end if
        alb_vis_dif = min(1.0_wp, alb_vis_dif - dalpha_vis_dif)

        ! --- diffuse near-IR: aging only, dust effect negligible
        alb_nir_dif = min(1.0_wp, 0.7493_wp + c%dalb_snow_nir - 0.1820_wp*rn - 0.0388_wp*rn**2)

        ! --- direct visible: zenith-corrected effective grain
        r  = snow_grain*(1.0_wp + 0.781_wp*(coszm - 0.65_wp)**2)
        rn = log10(r/r0)
        alb_vis_dir = 0.9849_wp + c%dalb_snow_vis - 0.0215_wp*rn - 0.0132_wp*rn**2
        if (dust_con .gt. dust_min) then
            f  = 155.0_wp + 17.15_wp*x + 0.27_wp*x**2
            cc = dust_con/f
            h  = cc/c0*(r/r0)**0.73_wp
            p  = log10(h)
            dalpha_vis_dir = 10.0_wp**(-0.049_wp*p**2 + 0.525_wp*p - 0.893_wp)
        else
            dalpha_vis_dir = 0.0_wp
        end if
        alb_vis_dir = min(1.0_wp, alb_vis_dir - dalpha_vis_dir)

        ! --- direct near-IR: its own zenith correction, aging only
        r  = snow_grain*(1.0_wp + 0.791_wp*(coszm - 0.65_wp)**2)
        rn = log10(r/r0)
        alb_nir_dir = min(1.0_wp, 0.6596_wp + c%dalb_snow_nir - 0.1927_wp*rn - 0.0229_wp*rn**2)

        ! --- orographic roughness reduction, all bands
        alb_vis_dif = alb_vis_dif - dalb_orog
        alb_vis_dir = alb_vis_dir - dalb_orog
        alb_nir_dif = alb_nir_dif - dalb_orog
        alb_nir_dir = alb_nir_dir - dalb_orog

        return

    end subroutine semix_bands_dang

    subroutine semix_bands_ww(snow_grain, dust_con, coszm, c, &
                                    alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
        ! Four-band snow albedo, Warren & Wiscombe 1980 form
        ! (smb_surface_par.f90:225-281). Diffuse bands age from the fresh values
        ! toward darker ones; direct bands add a solar-zenith brightening. The
        ! dust terms (c_dust_*) darken the visible band most; with dust_con = 0
        ! they vanish and this reduces to the clean-snow albedo.

        implicit none

        real(wp),                intent(IN)  :: snow_grain   ! [um] grain size
        real(wp),                intent(IN)  :: dust_con     ! [kg kg-1] dust in snow
        real(wp),                intent(IN)  :: coszm        ! [1] daily-mean cos(zenith)
        type(chion_const_class), intent(IN)  :: c
        real(wp),                intent(OUT) :: alb_vis_dir, alb_nir_dir
        real(wp),                intent(OUT) :: alb_vis_dif, alb_nir_dif

        ! Dust reduction lookup (Warren & Wiscombe 1980, Fig. 5), as in SEMIX.
        real(wp), parameter :: tab0(4) = (/1.001_wp,  10.0_wp, 100.0_wp, 1000.0_wp/)
        real(wp), parameter :: tab1(4) = (/0.00_wp,   0.02_wp,  0.10_wp,    0.30_wp/)
        real(wp), parameter :: tab2(4) = (/0.01_wp,   0.05_wp,  0.15_wp,    0.30_wp/)

        integer  :: k
        real(wp) :: f_cosz, f_age, d1, rint
        real(wp) :: c_dust_new_vis, c_dust_age_vis, c_dust_new_nir, c_dust_age_nir
        real(wp) :: c_age_vis, c_age_nir

        ! Solar-zenith factor (BATS-modified, tuned to Gardner & Sharp 2010).
        f_cosz = 0.5_wp*(3.0_wp/(1.0_wp + 2.0_wp*coszm) - 1.0_wp)
        f_cosz = max(0.0_wp, f_cosz)

        ! Dust darkening (visible), log-interpolated over the lookup table.
        d1 = min(dust_con*1.0e6_wp, 999.0_wp)
        if (d1 .gt. 1.0001_wp) then
            if (d1 .lt. 10.0_wp)   k = 1
            if (d1 .ge. 10.0_wp  .and. d1 .lt. 100.0_wp) k = 2
            if (d1 .ge. 100.0_wp) k = 3
            rint = (log(d1) - log(tab0(k)))/(log(tab0(k+1)) - log(tab0(k)))
            c_dust_new_vis = (1.0_wp - rint)*tab1(k) + rint*tab1(k+1)
            c_dust_age_vis = (1.0_wp - rint)*tab2(k) + rint*tab2(k+1)
        else
            c_dust_new_vis = 0.0_wp
            c_dust_age_vis = 0.0_wp
        end if
        ! NIR dust reduction is halved (the scheme overestimates it otherwise).
        c_dust_new_nir = 0.5_wp*c_dust_new_vis
        c_dust_age_nir = 0.5_wp*c_dust_age_vis

        c_age_vis = c%d_alb_age_vis + c_dust_age_vis
        c_age_nir = c%d_alb_age_nir + c_dust_age_nir

        ! Aging factor: 0 at fresh grain, 1 at old grain (log form as in SEMIX).
        f_age = log10(1.0_wp + (snow_grain - c%snow_grain_fresh)/200.0_wp) &
              / log10(1.0_wp + (c%snow_grain_old - c%snow_grain_fresh)/200.0_wp)

        alb_vis_dif = c%alb_snow_vis_new - f_age*c_age_vis - c_dust_new_vis
        alb_nir_dif = c%alb_snow_nir_new - f_age*c_age_nir - c_dust_new_nir
        alb_vis_dir = alb_vis_dif + 0.4_wp*f_cosz*(1.0_wp - alb_vis_dif)
        alb_nir_dir = alb_nir_dif + 0.4_wp*f_cosz*(1.0_wp - alb_nir_dif)

        return

    end subroutine semix_bands_ww

    pure function semix_broadband_albedo(alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
                                         cloud, frac_vu) result(alb)
        ! Collapse the four bands to a broadband albedo using the incoming-SW
        ! spectral weights: cloud fraction sets the direct/diffuse split
        ! (clear -> direct, overcast -> diffuse) and frac_vu the visible share.
        ! Weights sum to one, so alb is a proper broadband albedo in [0,1].

        implicit none

        real(wp), intent(IN) :: alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif
        real(wp), intent(IN) :: cloud     ! [1] cloud fraction, clamped to [0,1]
        real(wp), intent(IN) :: frac_vu   ! [1] visible+UV fraction of incoming solar
        real(wp) :: alb

        real(wp) :: cl, fv

        cl = min(max(cloud,   0.0_wp), 1.0_wp)
        fv = min(max(frac_vu, 0.0_wp), 1.0_wp)

        alb = (1.0_wp - cl)*(fv*alb_vis_dir + (1.0_wp - fv)*alb_nir_dir) &
            +          cl  *(fv*alb_vis_dif + (1.0_wp - fv)*alb_nir_dif)

        return

    end function semix_broadband_albedo

    pure function semix_daily_coszm(latitude_deg, solar_longitude_deg) result(coszm)
        ! Daylight-mean cosine of the solar zenith angle from latitude and solar
        ! longitude (chion already carries both per column). Declination from
        ! the obliquity relation sin(delta) = sin(eps)*sin(lambda_s); the
        ! daily mean over the sunlit part of the day is the standard closed form
        !   coszm = (H0 sinφ sinδ + cosφ cosδ sinH0) / H0 ,
        ! with H0 the half-day sunset hour angle. Used as the fallback when the
        ! host does not supply coszm.

        implicit none

        real(wp), intent(IN) :: latitude_deg
        real(wp), intent(IN) :: solar_longitude_deg
        real(wp) :: coszm

        real(wp), parameter :: pi  = 3.14159265358979_wp
        real(wp), parameter :: d2r = pi/180.0_wp
        real(wp), parameter :: obliquity = 23.44_wp*d2r
        real(wp) :: phi, delta, arg, h0

        phi   = latitude_deg*d2r
        delta = asin(sin(obliquity)*sin(solar_longitude_deg*d2r))

        ! Sunset hour angle; polar day/night clamp the argument to [-1,1].
        arg = -tan(phi)*tan(delta)
        arg = min(max(arg, -1.0_wp), 1.0_wp)
        h0  = acos(arg)

        if (h0 .le. 0.0_wp) then
            coszm = 0.0_wp   ! polar night
        else
            coszm = (h0*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(h0))/h0
            coszm = max(0.0_wp, coszm)
        end if

        return

    end function semix_daily_coszm

end module snow_albedo_semix
