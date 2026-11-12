module chion_forcing_host
    ! Forcing a host can derive when its climate supplies only monthly air
    ! temperature and precipitation (library level; Chion.jl src/forcing.jl has
    ! the same functions). None is called by chion itself.
    !
    ! - forcing_shortwave_cloud: the surface shortwave from the top-of-atmosphere
    !   one, SWdn = TOA*tau_clear(z)*(1 - cloud), with tau_clear the cloud proxy's
    !   clear-sky transmissivity (clear_sky_transmissivity), so the proxy's
    !   cloudiness 1 - SWdn/(TOA*tau_clear) is the given cloud fraction. MAR
    !   v3.14 1981-2010 over Greenland's ice: SWdn/TOA = 0.61 + 0.05 z[km] in
    !   summer, cloud = 0.31 (with chion's TOA); CERES 2001-2020: Greenland 0.33,
    !   Antarctica 0.23.
    ! - forcing_wind_speed: the wind from the surface elevation and slope on a
    !   regular grid, U = u0 + u_z*max(z,0)[km] + u_s*|grad z|[m/km], the slope a
    !   centred difference over a fixed baseline of about L (2n cells,
    !   n = max(1, nint(L/(2 dx)))), so that it does not depend on the grid's
    !   resolution. ERA5 1981-2010 summer 10-m wind over the ice: Greenland
    !   2.5 + 1.27 z[km], Antarctica 5.5 + 0.2 z[km] (u_s = 0).
    ! - forcing_tvar_sigma, forcing_red_noise_monthly: daily air temperature
    !   variability for a smooth (monthly-interpolated) daily series,
    !   T = T_smooth + sigma(T_mon)*xi. sigma falls with the monthly mean
    !   temperature, lowest over melting surfaces (Seguinot and Rogozhina 2014,
    !   J. Glaciol. 60, 603): sigma = clamp(sigma0 + slope*(T_mon - T0), smin,
    !   smax); ERA5 1981-2010 over the Greenland and Antarctic ice:
    !   2.0 - 0.14 (T_mon - T0) in [2, 6] K. xi is AR(1) red noise (e-folding
    !   tau; ERA5: 4-6 d) of zero mean in every month and unit variance, from a
    !   portable Park-Miller generator and Box-Muller, so a seed gives the same
    !   sequence with any compiler and in Chion.jl.

    use chion_defs,          only : wp, chion_const_class
    use snow_surface_fluxes, only : clear_sky_transmissivity

    implicit none

    private

    public :: forcing_shortwave_cloud
    public :: forcing_wind_speed
    public :: forcing_tvar_sigma
    public :: forcing_red_noise_monthly

contains

    elemental function forcing_shortwave_cloud(c,toa,height,cloud) result(sw_down)
        ! SWdn = TOA*tau_clear(z)*(1 - cloud) [W m-2].

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: toa       ! [W m-2] top-of-atmosphere shortwave
        real(wp),                intent(IN) :: height    ! [m] surface height
        real(wp),                intent(IN) :: cloud     ! [1] cloud fraction
        real(wp) :: sw_down

        sw_down = toa*clear_sky_transmissivity(c,height)*(1.0_wp - cloud)

        return

    end function forcing_shortwave_cloud

    function forcing_wind_speed(z_srf,dx,u0,u_z,u_s,length_km) result(u)
        ! U = u0 + u_z*max(z,0)[km] + u_s*|grad z|[m/km], at least 0 [m s-1], on a
        ! regular grid of spacing dx [m]; the gradient over 2n cells,
        ! n = max(1, nint(L/(2 dx))), one-sided over the cells available at the
        ! edges. u_s = 0 needs no dx.

        implicit none

        real(wp), intent(IN) :: z_srf(:,:)      ! [m]
        real(wp), intent(IN) :: dx              ! [m]
        real(wp), intent(IN) :: u0              ! [m s-1]
        real(wp), intent(IN) :: u_z             ! [m s-1 km-1]
        real(wp), intent(IN) :: u_s             ! [m s-1 per m/km]
        real(wp), intent(IN) :: length_km       ! [km] baseline of the slope
        real(wp) :: u(size(z_srf,1),size(z_srf,2))

        integer  :: i, j, n, nx, ny, i0, i1, j0, j1
        real(wp) :: gx, gy

        u = u0 + u_z*max(z_srf,0.0_wp)*1.0e-3_wp

        if (u_s .ne. 0.0_wp) then
            nx = size(z_srf,1)
            ny = size(z_srf,2)
            n  = max(1,nint(length_km*1.0e3_wp/(2.0_wp*dx)))
            do j = 1, ny
                j0 = max(j-n,1)
                j1 = min(j+n,ny)
                do i = 1, nx
                    i0 = max(i-n,1)
                    i1 = min(i+n,nx)
                    gx = 0.0_wp
                    gy = 0.0_wp
                    if (i1 .gt. i0) gx = (z_srf(i1,j) - z_srf(i0,j))/(real(i1-i0,wp)*dx)
                    if (j1 .gt. j0) gy = (z_srf(i,j1) - z_srf(i,j0))/(real(j1-j0,wp)*dx)
                    u(i,j) = u(i,j) + u_s*1.0e3_wp*sqrt(gx*gx + gy*gy)
                end do
            end do
        end if

        u = max(u,0.0_wp)

        return

    end function forcing_wind_speed

    elemental function forcing_tvar_sigma(t_mon,T0,sigma0,slope,smin,smax) result(sigma)
        ! sigma = min(max(sigma0 + slope*(t_mon - T0), smin), smax) [K].

        implicit none

        real(wp), intent(IN) :: t_mon       ! [K] monthly mean air temperature
        real(wp), intent(IN) :: T0          ! [K] melting point
        real(wp), intent(IN) :: sigma0      ! [K] sigma at t_mon = T0
        real(wp), intent(IN) :: slope       ! [K K-1]
        real(wp), intent(IN) :: smin        ! [K]
        real(wp), intent(IN) :: smax        ! [K]
        real(wp) :: sigma

        sigma = min(max(sigma0 + slope*(t_mon - T0),smin),smax)

        return

    end function forcing_tvar_sigma

    function forcing_red_noise_monthly(nday,dt_days,nday_mon,tau,seed) result(xi)
        ! One year of red noise, a value per step of dt_days: AR(1) with the
        ! lag-one correlation exp(-dt_days/tau), driven by standard normals
        ! (Park-Miller minimal standard uniforms x <- 16807 x mod (2^31 - 1),
        ! x0 = mod(|seed|, 2^31 - 2) + 1, two per normal, Box-Muller cosine
        ! branch); then its mean removed in every month of nday_mon days and the
        ! whole scaled to unit variance.

        implicit none

        integer,  intent(IN) :: nday        ! [d] days per year
        integer,  intent(IN) :: dt_days     ! [d] step
        integer,  intent(IN) :: nday_mon    ! [d] days per month
        real(wp), intent(IN) :: tau         ! [d] e-folding time
        integer,  intent(IN) :: seed
        real(wp) :: xi(nday/dt_days)

        integer, parameter :: i8 = selected_int_kind(18)
        integer(i8) :: x
        integer  :: n, k, m, nm, m1
        real(wp) :: r, u1, u2, g

        n = nday/dt_days
        x = int(mod(abs(seed),2147483646),i8) + 1_i8
        r = exp(-real(dt_days,wp)/max(tau,1.0e-6_wp))

        do k = 1, n
            x  = mod(16807_i8*x,2147483647_i8)
            u1 = real(x,wp)/2147483647.0_wp
            x  = mod(16807_i8*x,2147483647_i8)
            u2 = real(x,wp)/2147483647.0_wp
            g  = sqrt(-2.0_wp*log(u1))*cos(8.0_wp*atan(1.0_wp)*u2)
            if (k .eq. 1) then
                xi(k) = g
            else
                xi(k) = r*xi(k-1) + sqrt(1.0_wp - r*r)*g
            end if
        end do

        nm = max(1,nday_mon/dt_days)
        do m = 1, n, nm
            m1 = min(m+nm-1,n)
            xi(m:m1) = xi(m:m1) - sum(xi(m:m1))/real(m1-m+1,wp)
        end do

        xi = xi/sqrt(max(sum(xi*xi)/real(n,wp),tiny(1.0_wp)))

        return

    end function forcing_red_noise_monthly

end module chion_forcing_host
