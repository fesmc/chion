module snow_seb_semix
    ! CLIMBER-X SEMIX surface energy balance: aerodynamic turbulent exchange.
    !
    ! Port of CLIMBER-X src/smb/smb_surface_par.f90 (resistance, :396-438) and
    ! the surface-flux coefficients of src/smb/smb_ebal.f90 (ebal, :100-120);
    ! Willeit, Calov, Ganopolski. The saturation-humidity helpers are ports of
    ! q_sat_i / dqsat_dT_i from src/main/constants.f90.
    !
    ! SEMIX replaces BESSI's single bulk exchange coefficient D_sh by a
    ! Monin-Obukhov-style aerodynamic resistance r_a: a snow-depth-weighted
    ! roughness length, a neutral exchange coefficient from the log-law, and a
    ! bulk-Richardson stability correction. Both schemes end up with the same
    ! algebraic SHAPE -- a coefficient times an air-minus-surface difference --
    ! so nothing downstream of the coefficient changes:
    !
    !     BESSI    Q_sh = D_sh*(T_air - T_s)
    !     SEMIX    Q_sh = f_sh*(T_air - T_s),   f_sh = rho_a*cp_air/r_a
    !
    ! and likewise for the latent flux, where SEMIX works in specific humidity
    ! (f_lh = Ls*rho_a/r_a, times qsat(T_s) - q_air) instead of BESSI's vapour
    ! pressure difference.
    !
    ! The two are NOT interchangeable in practice. At CLIMBER-X defaults, deep
    ! snow, 5 m s-1 wind and 263 K air, f_sh is about 8.0 W m-2 K-1 with the
    ! surface 5 K COLDER than the air (stable, r_a = 168 s m-1) and about
    ! 19.9 W m-2 K-1 with it 5 K warmer (unstable, r_a = 67 s m-1), against
    ! BESSI's fixed D_sh = 10.0.
    !
    ! Which side dominates is an empirical question, and over Greenland it is
    ! the STABLE one: a melting surface is pinned at T0 while the summer air
    ! above it is warmer, so T_s < T_air is the melting-season norm and the
    ! unstable branch is rare. The net effect on GRL-16KM is therefore LESS
    ! sensible heating and less melt, not more (see docs/semix_port_scope.md).
    !
    ! COUPLING (docs/semix_port_scope.md, decision alpha): chion has no massless
    ! skin node -- the top firn layer IS the surface, and the balance is
    ! linearized into row 1 of the conduction matrix. SEMIX's per-flux num/denom
    ! map one-for-one onto chion's q_const/q_lin, so this module returns the
    ! coefficients and the CALLERS assemble them; the ground-flux num_g/denom_g
    ! of ebal has no analog here, because chion's row-1 conduction coupling
    ! already plays that role.
    !
    ! r_a depends on the surface temperature through the Richardson number, so
    ! it is recomputed at each of the three flux sites from whatever surface
    ! temperature that site is evaluating at (linearization point T^n for the
    ! energy solve, T0 for bare ice, T^{n+1} for the post-solve vapour mass),
    ! rather than being frozen once per step as in SEMIX. The stability
    ! correction is weak enough at daily steps that the two agree closely, and
    ! this keeps each site self-contained.
    !
    ! SNOW DEPTH: SEMIX uses h_snow = w_snow/rho_snow with a fixed snow density
    ! (semi.f90:247). chion carries per-layer densities, so callers pass the
    ! true column depth sum(mass/density). Only the roughness blend uses it,
    ! and it saturates quickly (fsnow = h/(h+10*z0m_ice), i.e. half-way at 2 cm).

    use chion_defs, only : wp, chion_const_class, SEMIX_QSAT_BESSI
    use snow_vapor, only : safe_positive, relative_humidity_fraction, &
                           ice_saturation_vapor_pressure, &
                           ice_saturation_vapor_pressure_derivative

    implicit none

    private

    ! Dry-air / water-vapour molar mass ratio, as CLIMBER-X spells it in
    ! q_sat_i and dqsat_dT_i (and as chion spells it in
    ! latent_exchange_coefficient).
    real(wp), parameter, public :: SEMIX_EPS_VAPOR = 0.622_wp

    type semix_flux_lin_class
        ! One flux linearized about a reference temperature, in chion's
        ! convention  Q(T) = constant - linear*T,  so that the pair maps
        ! straight onto q_const / q_lin. This is exactly ebal's num/denom for
        ! that flux (smb_ebal.f90:122-129).
        real(wp) :: constant           ! [W m-2]      ebal num_*
        real(wp) :: linear             ! [W m-2 K-1]  ebal denom_*
    end type semix_flux_lin_class

    type semix_exchange_class
        ! Turbulent exchange coefficients at one surface temperature. f_sh and
        ! f_lh are ebal's (smb_ebal.f90:101-110); qsat/dqsatdT are the
        ! saturation humidity and its derivative at that same temperature, kept
        ! so the caller can build either the exact flux or its linearization
        ! without recomputing them.
        real(wp) :: r_a          ! [s m-1]        aerodynamic resistance
        real(wp) :: f_sh         ! [W m-2 K-1]    sensible heat coefficient
        real(wp) :: f_lh         ! [W m-2]        latent heat coefficient, per unit q
        real(wp) :: qsat         ! [kg kg-1]      saturation humidity at the surface
        real(wp) :: dqsatdT      ! [kg kg-1 K-1]  its temperature derivative
        real(wp) :: q_air        ! [kg kg-1]      air specific humidity
    end type semix_exchange_class

    public :: semix_flux_lin_class
    public :: semix_exchange_class
    public :: semix_snow_depth
    public :: semix_resistance
    public :: semix_air_density
    public :: semix_q_sat
    public :: semix_dqsat_dT
    public :: semix_air_humidity
    public :: semix_turbulent_exchange
    public :: semix_sensible_heat_flux
    public :: semix_latent_heat_flux
    public :: semix_surface_emissivity
    public :: semix_longwave_down
    public :: semix_longwave_linearized
    public :: semix_longwave_flux

contains

    ! =====================================================================
    ! Aerodynamic resistance
    ! =====================================================================

    pure function semix_snow_depth(mass,density,n) result(h_snow)
        ! Physical depth of the snow column, the only column property the
        ! aerodynamic scheme needs. SEMIX gets it as w_snow/rho_snow with a
        ! single fixed density; chion sums the true per-layer thicknesses.
        !
        ! Returns zero for an empty column, which is the bare-ice case: the
        ! roughness blend then collapses to z0m_ice alone.

        implicit none

        real(wp), intent(IN) :: mass(:)      ! (Ntot) [kg m-2]
        real(wp), intent(IN) :: density(:)   ! (Ntot) [kg m-3]
        integer,  intent(IN) :: n            ! active layer count
        real(wp) :: h_snow                   ! [m]

        ! Local variables
        integer :: k

        h_snow = 0.0_wp

        do k = 1, n
            h_snow = h_snow + mass(k)/safe_positive(density(k))
        end do

        return

    end function semix_snow_depth

    pure function semix_resistance(h_snow,t2m,t_skin,wind,c) result(r_a)
        ! CLIMBER-X smb_surface_par.f90:396-438, expression for expression.
        !
        ! The roughness blend weights snow against ice by a saturating function
        ! of snow depth; the neutral exchange coefficient is the product of the
        ! momentum and heat log-law factors; the stability correction is a bulk
        ! Richardson number, applied only on the unstable side.
        !
        ! Returns a very large resistance (1e20) when Ch*wind is not positive,
        ! which is SEMIX's way of switching the turbulent fluxes off in a dead
        ! calm rather than dividing by zero.

        implicit none

        real(wp),                intent(IN) :: h_snow   ! [m] snow depth
        real(wp),                intent(IN) :: t2m      ! [K] air temperature
        real(wp),                intent(IN) :: t_skin   ! [K] surface temperature
        real(wp),                intent(IN) :: wind     ! [m s-1]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: r_a                                 ! [s m-1]

        ! Local variables
        real(wp) :: fsnow, rough_m, rough_h
        real(wp) :: Ch, Ch_neutral, Ri

        fsnow   = h_snow/safe_positive(h_snow + 10.0_wp*c%z0m_ice)
        rough_m = fsnow*c%z0m_snow + (1.0_wp - fsnow)*c%z0m_ice
        rough_h = rough_m*c%zm_to_zh

        ! Neutral exchange coefficient
        Ch_neutral = (c%karman/log(c%z_sfl/rough_m))*(c%karman/log(c%z_sfl/rough_h))

        ! Richardson stability number
        Ri = c%grav*c%z_sfl*(1.0_wp - t_skin/safe_positive(t2m))/safe_positive(wind*wind)

        if (c%l_neutral) then
            Ch = Ch_neutral
        else
            if (Ri .lt. 0.0_wp) then
                ! "Unstable" stratification: enhanced exchange
                Ch = Ch_neutral*(1.0_wp - 2.0_wp*Ri)
            else
                ! "Stable" stratification: SEMIX leaves the neutral value alone
                Ch = Ch_neutral
            end if
        end if

        if (Ch*wind .gt. 0.0_wp) then
            r_a = 1.0_wp/(Ch*wind)
        else
            r_a = 1.0e20_wp
        end if

        return

    end function semix_resistance

    pure function semix_air_density(t2m,air_pressure,c) result(rhoa)
        ! CLIMBER-X constants.f90 rho_a: the ideal gas law for dry air.

        implicit none

        real(wp),                intent(IN) :: t2m           ! [K]
        real(wp),                intent(IN) :: air_pressure  ! [Pa]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: rhoa                                     ! [kg m-3]

        rhoa = air_pressure/(c%R_dry*safe_positive(t2m))

        return

    end function semix_air_density

    ! =====================================================================
    ! Saturation specific humidity over ice
    !
    ! Two parameterizations, selected by c%semix_qsat:
    !
    !   SEMIX_QSAT_SEMIX  CLIMBER-X q_sat_i / dqsat_dT_i (constants.f90), the
    !                     approximate forms actually used by SEMIX.
    !   SEMIX_QSAT_BESSI  chion's own ice saturation vapour pressure
    !                     (energy_flux.jl:62-71), converted to specific
    !                     humidity by the same 0.622/p that BESSI's
    !                     latent_exchange_coefficient folds in.
    !
    ! The two differ only in their fit coefficients (22.587/(T+0.71) against
    ! 22.46/(T-T0+272.62)) and are within about 1% of each other over the
    ! snowpack temperature range; the option exists so the SEMIX SEB can be
    ! run against either without a rebuild.
    ! =====================================================================

    pure function semix_q_sat(temperature,air_pressure,c) result(q)
        ! Saturation specific humidity over ice.

        implicit none

        real(wp),                intent(IN) :: temperature   ! [K]
        real(wp),                intent(IN) :: air_pressure  ! [Pa]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: q                                        ! [kg kg-1]

        ! Local variables
        real(wp) :: es

        if (c%semix_qsat .eq. SEMIX_QSAT_BESSI) then
            es = ice_saturation_vapor_pressure(temperature,c%T0)
            q  = SEMIX_EPS_VAPOR*es/safe_positive(air_pressure)
        else
            ! constants.f90 q_sat_i, with the 380.1726 = 611.2*0.622 prefactor
            ! already carrying the vapour-pressure conversion.
            q = 380.1726_wp*exp(22.587_wp*(temperature - c%T0) &
                                /safe_positive(temperature + 0.71_wp)) &
                /safe_positive(air_pressure)
        end if

        return

    end function semix_q_sat

    pure function semix_dqsat_dT(temperature,air_pressure,c) result(dqdT)
        ! d(q_sat)/dT over ice, matching the parameterization semix_q_sat used.
        !
        ! The SEMIX branch is constants.f90 dqsat_dT_i: Clausius-Clapeyron on
        ! the (unapproximated) saturation vapour pressure e_sat_i, divided by p.
        ! Note CLIMBER-X evaluates e_sat_i with the 273.86-offset form while
        ! q_sat_i itself uses the T+0.71 approximation -- the pair is
        ! deliberately not self-consistent in the source and is reproduced as
        ! written.

        implicit none

        real(wp),                intent(IN) :: temperature   ! [K]
        real(wp),                intent(IN) :: air_pressure  ! [Pa]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: dqdT                                     ! [kg kg-1 K-1]

        ! Local variables. Ls_ice and R_vapor are CLIMBER-X's dqsat_dT_i
        ! locals (its Lv = 2834e3 J kg-1 for ice, and Rv from constants.f90).
        real(wp), parameter :: LS_ICE  = 2834.0e3_wp   ! [J kg-1]
        real(wp), parameter :: R_VAPOR = 461.5_wp      ! [J kg-1 K-1]

        real(wp) :: es, des_dT

        if (c%semix_qsat .eq. SEMIX_QSAT_BESSI) then
            es     = ice_saturation_vapor_pressure(temperature,c%T0)
            des_dT = ice_saturation_vapor_pressure_derivative(temperature,c%T0,es)
            dqdT   = SEMIX_EPS_VAPOR*des_dT/safe_positive(air_pressure)
        else
            es     = semix_e_sat_i(temperature,c%T0)
            des_dT = LS_ICE*es/(R_VAPOR*safe_positive(temperature*temperature))
            dqdT   = SEMIX_EPS_VAPOR*des_dT/safe_positive(air_pressure)
        end if

        return

    end function semix_dqsat_dT

    pure function semix_e_sat_i(temperature,T0) result(es)
        ! CLIMBER-X constants.f90 e_sat_i: saturation vapour pressure over ice.
        ! Private -- only dqsat_dT_i consumes it.

        implicit none

        real(wp), intent(IN) :: temperature   ! [K]
        real(wp), intent(IN) :: T0            ! [K]
        real(wp) :: es                        ! [Pa]

        ! Local variables
        real(wp) :: tc

        tc = temperature - T0
        es = 611.2_wp*exp(22.587_wp*tc/(273.86_wp + tc))

        return

    end function semix_e_sat_i

    pure function semix_air_humidity(t2m,relative_humidity,air_pressure,c) result(q_air)
        ! Air specific humidity from chion's relative-humidity forcing.
        !
        ! semi.f90:201 forms q2m = r2m*q_sat_i(t2m,pressure), i.e. relative
        ! humidity WITH RESPECT TO ICE. The r2m there is a downscaling product
        ! (a blend of the free-atmosphere and skin relative humidities) that
        ! stays host-side in chion's contract, so the forcing value is used
        ! directly. Percent-valued input is handled by chion's existing
        ! relative_humidity_fraction.
        !
        ! NOTE this differs from BESSI, which takes the air vapour pressure
        ! relative to saturation over WATER (energy_flux.jl:57-60). Under the
        ! SEMIX scheme the same forcing number therefore means something
        ! slightly different -- deliberately, since it is SEMIX's own reading.

        implicit none

        real(wp),                intent(IN) :: t2m                ! [K]
        real(wp),                intent(IN) :: relative_humidity  ! [1] or [%]
        real(wp),                intent(IN) :: air_pressure       ! [Pa]
        type(chion_const_class), intent(IN) :: c
        real(wp) :: q_air                                         ! [kg kg-1]

        q_air = relative_humidity_fraction(relative_humidity) &
                *semix_q_sat(t2m,air_pressure,c)

        return

    end function semix_air_humidity

    ! =====================================================================
    ! Turbulent exchange coefficients
    ! =====================================================================

    pure function semix_turbulent_exchange(c,h_snow,t2m,t_skin,wind,air_pressure, &
                                           relative_humidity,has_humidity) result(x)
        ! smb_ebal.f90:100-110. The dew-inhibition switch (l_dew) zeroes the
        ! latent coefficient whenever the air is more humid than the saturated
        ! surface, which is exactly the condition for deposition.
        !
        ! With no humidity forcing the latent coefficient is zero, so the
        ! latent flux vanishes -- the same three-way selection BESSI makes
        ! (energy_flux.jl:379, surface_fluxes.jl:37-43).
        !
        ! SEMIX uses the latent heat of SUBLIMATION for every latent exchange,
        ! whatever the surface temperature. chion has no Ls constant; Lv + Lm
        ! is the same quantity and is what BESSI's own exchange coefficient
        ! already uses.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: h_snow             ! [m]
        real(wp),                intent(IN) :: t2m                ! [K]
        real(wp),                intent(IN) :: t_skin             ! [K]
        real(wp),                intent(IN) :: wind               ! [m s-1]
        real(wp),                intent(IN) :: air_pressure       ! [Pa]
        real(wp),                intent(IN) :: relative_humidity  ! [1] or [%]
        logical,                 intent(IN) :: has_humidity
        type(semix_exchange_class) :: x

        ! Local variables
        real(wp) :: rhoa

        x%r_a  = semix_resistance(h_snow,t2m,t_skin,wind,c)
        rhoa   = semix_air_density(t2m,air_pressure,c)

        x%f_sh = rhoa*c%cp_air/x%r_a

        ! The saturation humidity and its derivative feed the LATENT flux and
        ! nothing else. With no humidity forcing that flux is identically zero
        ! and every consumer multiplies these by f_lh, so evaluating them would
        ! be two wasted exp() -- per call, per flux site, per substep, per
        ! column. Returning early instead is measurable on a whole-domain run.
        !
        ! They are zeroed rather than left undefined: qsat and dqsatdT are only
        ! meaningful when f_lh is non-zero, and a caller that ignores that gets
        ! a harmless zero instead of a stale value.
        if (.not. has_humidity) then
            x%qsat    = 0.0_wp
            x%dqsatdT = 0.0_wp
            x%q_air   = 0.0_wp
            x%f_lh    = 0.0_wp
            return
        end if

        x%qsat    = semix_q_sat(t_skin,air_pressure,c)
        x%dqsatdT = semix_dqsat_dT(t_skin,air_pressure,c)
        x%q_air   = semix_air_humidity(t2m,relative_humidity,air_pressure,c)

        if (.not. c%l_dew .and. x%q_air .gt. x%qsat) then
            ! Inhibit dew/frost deposition
            x%f_lh = 0.0_wp
        else
            x%f_lh = (c%Lv + c%Lm)/x%r_a*rhoa
        end if

        return

    end function semix_turbulent_exchange

    ! =====================================================================
    ! Longwave
    !
    ! ebal num_lw/denom_lw (smb_ebal.f90:128-129). ONE thing distinguishes
    ! these from the BESSI expressions they replace: the DOWNWELLING flux is
    ! absorbed with the surface emissivity, emiss*lwdown, rather than in full.
    ! BESSI takes lwdown at face value and applies eps_snow only to the
    ! emitted term, which is an absorptivity of 1 against an emissivity of
    ! 0.98 -- not Kirchhoff-consistent, but it is what Chion.jl does and the
    ! bessi scheme keeps it.
    !
    ! At eps = 0.98 and a typical 250 W m-2 of downwelling longwave the
    ! difference is about 5 W m-2 LESS energy into the surface, in one
    ! direction, all year.
    !
    ! SEMIX also carries separate snow and ice emissivities (both 0.99 in
    ! CLIMBER-X). chion reuses its existing eps_snow and adds eps_ice
    ! alongside it, both defaulting to chion's 0.98; setting the pair to 0.99
    ! reproduces SEMIX exactly.
    ! =====================================================================

    pure function semix_surface_emissivity(c,has_snow) result(emiss)
        ! ebal's mask_snow branch (smb_ebal.f90:88-95).

        implicit none

        type(chion_const_class), intent(IN) :: c
        logical,                 intent(IN) :: has_snow
        real(wp) :: emiss                                     ! [1]

        if (has_snow) then
            emiss = c%eps_snow
        else
            emiss = c%eps_ice
        end if

        return

    end function semix_surface_emissivity

    pure function semix_longwave_down(c,forc_q_lw_down,has_q_lw_down,t2m) result(lw_down)
        ! Resolve the downwelling longwave the same way both schemes do:
        ! prescribed if the host supplies it, otherwise chion's grey-air
        ! parameterization eps_air*sigma*T_air^4 (energy_flux.jl:374).
        !
        ! SEMIX always receives lwdown from its atmosphere and has no fallback;
        ! the fallback is chion's, kept so the SEMIX scheme still runs in
        ! minimal-input mode.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: forc_q_lw_down   ! [W m-2]
        logical,                 intent(IN) :: has_q_lw_down
        real(wp),                intent(IN) :: t2m              ! [K]
        real(wp) :: lw_down                                     ! [W m-2]

        if (has_q_lw_down) then
            lw_down = forc_q_lw_down
        else
            lw_down = c%sigma_sb*c%eps_air*t2m**4
        end if

        return

    end function semix_longwave_down

    pure function semix_longwave_linearized(c,emiss,lw_down,t_ref) result(coef)
        ! num_lw / denom_lw, linearized about t_ref, in chion's convention
        !     Q_lw(T) = constant - linear*T
        ! so the two map straight onto q_const / q_lin.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: emiss     ! [1]
        real(wp),                intent(IN) :: lw_down   ! [W m-2]
        real(wp),                intent(IN) :: t_ref     ! [K] linearization point
        type(semix_flux_lin_class) :: coef

        coef%constant = emiss*lw_down + 3.0_wp*emiss*c%sigma_sb*t_ref**4
        coef%linear   = 4.0_wp*emiss*c%sigma_sb*t_ref**3

        return

    end function semix_longwave_linearized

    pure function semix_longwave_flux(c,emiss,lw_down,t_surface) result(q_lw)
        ! EXACT net longwave at a known surface temperature, positive into the
        ! surface. Reproduces semix_longwave_linearized at t_ref = t_surface.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: emiss        ! [1]
        real(wp),                intent(IN) :: lw_down      ! [W m-2]
        real(wp),                intent(IN) :: t_surface    ! [K]
        real(wp) :: q_lw                                    ! [W m-2]

        q_lw = emiss*lw_down - emiss*c%sigma_sb*t_surface**4

        return

    end function semix_longwave_flux

    ! =====================================================================
    ! Fluxes, in chion's sign convention
    !
    ! chion counts a flux POSITIVE INTO THE SURFACE; SEMIX's diagnostics
    ! (smb_ebal.f90:139-142) count sensible and latent positive AWAY from it.
    ! The sign is flipped here, once, so callers never have to think about it.
    ! =====================================================================

    pure function semix_sensible_heat_flux(x,t2m,t_surface) result(q_sh)
        ! Positive into the surface: f_sh*(T_air - T_s).

        implicit none

        type(semix_exchange_class), intent(IN) :: x
        real(wp),                   intent(IN) :: t2m         ! [K]
        real(wp),                   intent(IN) :: t_surface   ! [K]
        real(wp) :: q_sh                                      ! [W m-2]

        q_sh = x%f_sh*(t2m - t_surface)

        return

    end function semix_sensible_heat_flux

    pure function semix_latent_heat_flux(x) result(q_lh)
        ! Positive into the surface (deposition): -f_lh*(qsat(T_s) - q_air).
        !
        ! EXACT at the temperature the exchange coefficients were built at:
        ! x%qsat is already qsat(T_s), so no linearization is involved. The
        ! linearized form the energy solve needs is assembled by the caller
        ! from f_lh, qsat and dqsatdT, mirroring ebal's num_lh/denom_lh.

        implicit none

        type(semix_exchange_class), intent(IN) :: x
        real(wp) :: q_lh                                      ! [W m-2]

        q_lh = -x%f_lh*(x%qsat - x%q_air)

        return

    end function semix_latent_heat_flux

end module snow_seb_semix
