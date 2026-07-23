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
    ! SCOPE (this commit): clean snow only -- fresh grain, no dust. Grain aging
    ! and the dust-in-snow darkening are added in the following commits; the
    ! band formula below already carries the aging (f_age) and dust terms so
    ! those commits only supply the grain-size and dust-concentration state.

    use chion_defs, only : wp, chion_const_class, chion_step_forcing_class

    implicit none

    private

    public :: semix_surface_albedo
    public :: semix_snow_grain_size
    public :: semix_snow_albedo_bands
    public :: semix_broadband_albedo
    public :: semix_daily_coszm

contains

    subroutine semix_surface_albedo(forc, c, t_skin, dust_con, albedo)
        ! Broadband SEMIX snow albedo for one column. Grain size is diagnosed
        ! from the skin temperature and the current snowfall rate; coszm is
        ! taken from the forcing when supplied, else from the lat/solar-longitude
        ! daily mean; cloud defaults to clear-sky (all-direct) when absent. Dust
        ! concentration is passed in: this commit calls it with zero dust (the
        ! dust-in-snow darkening is added in the following commit).

        implicit none

        type(chion_step_forcing_class), intent(IN)  :: forc
        type(chion_const_class),        intent(IN)  :: c
        real(wp),                       intent(IN)  :: t_skin       ! [K] surface temperature
        real(wp),                       intent(IN)  :: dust_con     ! [kg kg-1]
        real(wp),                       intent(OUT) :: albedo       ! [1] broadband

        real(wp) :: coszm, cloud, snow_grain
        real(wp) :: av_dir, an_dir, av_dif, an_dif

        snow_grain = semix_snow_grain_size(t_skin, forc%snowfall_rate, c)

        if (forc%has_coszm) then
            coszm = forc%coszm
        else
            coszm = semix_daily_coszm(forc%latitude_deg, forc%solar_longitude_deg)
        end if

        cloud = 0.0_wp
        if (forc%has_cloud) cloud = forc%cloud

        call semix_snow_albedo_bands(snow_grain, dust_con, coszm, c, &
                                     av_dir, an_dir, av_dif, an_dif)

        albedo = semix_broadband_albedo(av_dir, an_dir, av_dif, an_dif, cloud, c%frac_vu)

        return

    end subroutine semix_surface_albedo

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

    subroutine semix_snow_albedo_bands(snow_grain, dust_con, coszm, c, &
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

    end subroutine semix_snow_albedo_bands

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
