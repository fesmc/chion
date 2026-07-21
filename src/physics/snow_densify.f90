module snow_densify
    ! Firn densification: BESSI / HTESSEL low-density branches plus the shared
    ! mid- and high-density parameterizations, and the HTESSEL liquid-water
    ! compaction correction.
    !
    ! Port of Chion.jl/src/processes/densification.jl.
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n, per docs/porting_notes.md D8.
    !
    ! MAGIC CONSTANTS CARRIED OVER VERBATIM (docs/PLAN.md section 5, item 3).
    ! Two of these disagree with the constants chion uses everywhere else, and
    ! that disagreement is NOT resolved here:
    !
    !   DENSIFY_GRAVITY   = 9.81   while DEF_GRAVITY               = 9.80665
    !   DENSIFY_R_GAS     = 8.13   while DEF_UNIVERSAL_GAS_CONSTANT = 8.31446...
    !
    ! docs/PLAN.md section 4.1 permits reconciling both, BUT ONLY after
    ! measuring the impact on a decade-long column run and recording it. That
    ! measurement has not been done (it needs WP8), so both are ported exactly
    ! as Chion.jl has them and are named, not inlined, so the experiment is a
    ! one-line change when WP8 exists.
    !
    !   * 9.81 vs 9.80665 is a 3.6e-4 relative change in overburden, entering
    !     the mid/high branches cubed -> ~1.1e-3 relative on those tendencies.
    !   * 8.13 vs 8.314 is the load-bearing one: it sits inside an Arrhenius
    !     exponent, so exp(-10160/(8.13*T)) / exp(-10160/(8.314*T)) at T=260 K
    !     is a factor of ~2.6. Substituting the gas constant would change the
    !     low-density BESSI rate by more than a factor of two. This is a
    !     modelling question, not a typo to fix silently.
    !
    ! PRESERVED QUIRKS:
    !   * The overburden is sigma = g*(M_above + m/2), with M_above accumulated
    !     AFTER the current layer is processed (so layer 1 sees only half its
    !     own mass). Skipped layers still contribute to M_above.
    !   * The scheme flag selects ONLY the rho < 550 branch. The mid and high
    !     branches are shared.
    !   * Trap 8: Chion.jl dispatches with "if htessel ... else bessi", so any
    !     unrecognized flag silently selects BESSI. chion_defs validates the
    !     flag at parse time (docs/porting_notes.md D5), so this module uses a
    !     select case with an error default instead.
    !   * The update is monotone by construction: max(rho, rho + drho*dt), then
    !     capped at rho_i. A negative tendency therefore cannot thin a layer.

    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

    use chion_defs, only : wp, wp_acc, TOL_TINY, chion_const_class, &
                           io_unit_err, CHION_DENSIFY_BESSI, CHION_DENSIFY_HTESSEL

    implicit none

    private

    ! --- Constants that differ from chion's standard ones. See header. -----
    real(wp_acc), parameter, public :: DENSIFY_GRAVITY = 9.81_wp_acc      ! [m s-2] NOT DEF_GRAVITY
    real(wp_acc), parameter, public :: DENSIFY_R_GAS   = 8.13_wp_acc      ! [-]     NOT the gas constant

    ! --- Density regime thresholds -----------------------------------------
    real(wp), parameter, public :: DENSIFY_RHO_LOW = 550.0_wp   ! [kg m-3] low  | mid boundary
    real(wp), parameter, public :: DENSIFY_RHO_MID = 800.0_wp   ! [kg m-3] mid  | high boundary

    ! --- Hard-coded firn parameters (Chion.jl densification.jl:202-203) -----
    real(wp_acc), parameter, public :: DENSIFY_RHO_E = 815.0_wp_acc     ! [kg m-3] close-off density
    real(wp_acc), parameter, public :: DENSIFY_P_ATM = 101325.0_wp_acc  ! [Pa] atmospheric pressure

    public :: bessi_low_density_rate
    public :: htessel_low_density_rate
    public :: bubble_pressure_mpa
    public :: mid_density_tendency
    public :: high_density_tendency
    public :: densify_column
    public :: apply_htessel_liquid_water_compaction

contains

    pure function bessi_low_density_rate(density,temperature,rho_i,accumulation_rate) result(drho)
        ! Chion.jl/src/processes/densification.jl:5-11:
        !     0.011*exp(-10160/(8.13*T))*(rho_i - rho)*max(A_t, 0)
        ! NOTE 8.13, not the gas constant. See the module header.

        implicit none

        real(wp), intent(IN) :: density            ! [kg m-3]
        real(wp), intent(IN) :: temperature        ! [K]
        real(wp), intent(IN) :: rho_i              ! [kg m-3]
        real(wp), intent(IN) :: accumulation_rate  ! [kg m-2 s-1] accumulation proxy A_t
        real(wp) :: drho                           ! [kg m-3 s-1]

        drho = real(0.011_wp_acc &
                    *exp(-10160.0_wp_acc/(DENSIFY_R_GAS*real(temperature,wp_acc))) &
                    *real(rho_i - density,wp_acc) &
                    *max(real(accumulation_rate,wp_acc),0.0_wp_acc), wp)

        return

    end function bessi_low_density_rate

    pure function htessel_low_density_rate(c,overburden_pressure,temperature,density) result(drho)
        ! Chion.jl/src/processes/densification.jl:13-41:
        !     eta  = 3.7e7*exp(8.1e-2*(T0-T) + 1.8e-2*rho)      snow viscosity
        !     xi   = 2.8e-6*exp(-4.2e-2*(T0-T) - 460*max(0,rho-150))   thermal metamorphism
        !     drho = rho*(sigma/eta + xi)
        !
        ! NOTE the 460 coefficient multiplies a density excess in kg m-3, so xi
        ! underflows to exactly zero for any rho above ~150.15 kg m-3. That is
        ! upstream behaviour, not a transcription error; it means xi only ever
        ! acts on the very freshest snow. Reported upstream.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: overburden_pressure  ! [Pa]
        real(wp),                intent(IN) :: temperature          ! [K]
        real(wp),                intent(IN) :: density              ! [kg m-3]
        real(wp) :: drho                                            ! [kg m-3 s-1]

        ! Local variables
        real(wp_acc) :: dT, eta, xi

        dT = real(c%T0,wp_acc) - real(temperature,wp_acc)

        eta = 3.7e7_wp_acc*exp( 8.1e-2_wp_acc*dT + 1.8e-2_wp_acc*real(density,wp_acc))
        xi  = 2.8e-6_wp_acc*exp(-4.2e-2_wp_acc*dT &
                                -460.0_wp_acc*max(0.0_wp_acc,real(density,wp_acc)-150.0_wp_acc))

        drho = real(real(density,wp_acc)*(real(overburden_pressure,wp_acc)/eta + xi),wp)

        return

    end function htessel_low_density_rate

    pure function bubble_pressure_mpa(density,rho_i) result(p_bubble)
        ! Chion.jl/src/processes/densification.jl:90-98:
        !     0                                       for rho <= rho_e
        !     P_atm*((1/rho_e - 1/rho_i)/(1/rho - 1/rho_i) - 1)/1e6  otherwise
        !
        ! Computed in wp_acc: both the numerator and the denominator are
        ! differences of nearly equal reciprocals, and the denominator tends to
        ! zero as rho -> rho_i (docs/PLAN.md section 3.1).

        implicit none

        real(wp), intent(IN) :: density     ! [kg m-3]
        real(wp), intent(IN) :: rho_i       ! [kg m-3]
        real(wp_acc) :: p_bubble            ! [MPa]

        ! Local variables
        real(wp_acc) :: rho, rhoi, num, den

        rho  = real(density,wp_acc)
        rhoi = real(rho_i,wp_acc)

        p_bubble = 0.0_wp_acc
        if (rho .le. DENSIFY_RHO_E) return

        num = 1.0_wp_acc/DENSIFY_RHO_E - 1.0_wp_acc/rhoi
        den = 1.0_wp_acc/rho           - 1.0_wp_acc/rhoi

        p_bubble = DENSIFY_P_ATM*(num/den - 1.0_wp_acc)/1.0e6_wp_acc

        return

    end function bubble_pressure_mpa

    pure function mid_density_tendency(density,temperature,dP,rho_i) result(drho)
        ! Chion.jl/src/processes/densification.jl:138-153.
        !     r    = rho/rho_i
        !     f    = 10**(-29.166*r**3 + 84.422*r**2 - 87.425*r + 30.673)
        !     drho = 25400*exp(-60000/(8.13*T))*rho*f*dP**3
        ! NOTE 8.13 again; see the module header.

        implicit none

        real(wp),     intent(IN) :: density      ! [kg m-3]
        real(wp),     intent(IN) :: temperature  ! [K]
        real(wp_acc), intent(IN) :: dP           ! [MPa] pressure excess
        real(wp),     intent(IN) :: rho_i        ! [kg m-3]
        real(wp) :: drho                         ! [kg m-3 s-1]

        ! Local variables
        real(wp_acc) :: r, f

        r = real(density,wp_acc)/real(rho_i,wp_acc)

        f = 10.0_wp_acc**( -29.166_wp_acc*r**3 &
                          + 84.422_wp_acc*r**2 &
                          - 87.425_wp_acc*r    &
                          + 30.673_wp_acc)

        drho = real(25400.0_wp_acc &
                    *exp(-60000.0_wp_acc/(DENSIFY_R_GAS*real(temperature,wp_acc))) &
                    *real(density,wp_acc)*f*dP**3, wp)

        return

    end function mid_density_tendency

    pure function high_density_tendency(density,temperature,dP,rho_i) result(drho)
        ! Chion.jl/src/processes/densification.jl:160-176.
        !     phi  = clamp(1 - rho/rho_i, 0, 1)
        !     den  = 1 - phi**(1/3)
        !     if abs(den) <= TOL_TINY -> 0
        !     f    = (3/16)*phi/den**3
        !     drho = 25400*exp(-60000/(8.13*T))*rho*f*dP**3
        !
        ! The porosity and the 1 - phi**(1/3) difference are computed in wp_acc,
        ! so that the abs(den) <= TOL_TINY guard can actually fire (see
        ! docs/porting_notes.md D1b: a dp tolerance is only meaningful against a
        ! dp-computed quantity).

        implicit none

        real(wp),     intent(IN) :: density      ! [kg m-3]
        real(wp),     intent(IN) :: temperature  ! [K]
        real(wp_acc), intent(IN) :: dP           ! [MPa] pressure excess
        real(wp),     intent(IN) :: rho_i        ! [kg m-3]
        real(wp) :: drho                         ! [kg m-3 s-1]

        ! Local variables
        real(wp_acc) :: phi, den, f

        phi = min(max(1.0_wp_acc - real(density,wp_acc)/real(rho_i,wp_acc), &
                      0.0_wp_acc),1.0_wp_acc)

        den = 1.0_wp_acc - phi**(1.0_wp_acc/3.0_wp_acc)

        drho = 0.0_wp
        if (abs(den) .le. TOL_TINY) return

        f = (3.0_wp_acc/16.0_wp_acc)*phi/den**3

        drho = real(25400.0_wp_acc &
                    *exp(-60000.0_wp_acc/(DENSIFY_R_GAS*real(temperature,wp_acc))) &
                    *real(density,wp_acc)*f*dP**3, wp)

        return

    end function high_density_tendency

    subroutine densify_column(mass,density,temperature,n,c,accumulation_rate,dt_seconds)
        ! Chion.jl/src/processes/densification.jl:185-296. Advance all active
        ! layer densities by one step. Mutates density only.
        !
        ! Per-layer flow:
        !   non-finite rho or T          -> skip (but still count toward M_above)
        !   rho >= rho_i - TOL_TINY      -> snap to rho_i, skip
        !   mass <= 0                    -> tendency 0 (still monotone-updated)
        !   rho <  550                   -> scheme-dependent low-density branch
        !   550 <= rho < 800             -> mid branch
        !   rho >= 800                   -> high branch
        ! then rho <- min(rho_i, max(rho, rho + drho*dt)).

        implicit none

        real(wp),                intent(IN)    :: mass(:)         ! (Ntot) [kg m-2] solid
        real(wp),                intent(INOUT) :: density(:)      ! (Ntot) [kg m-3]
        real(wp),                intent(IN)    :: temperature(:)  ! (Ntot) [K]
        integer,                 intent(IN)    :: n               ! active layers
        type(chion_const_class), intent(IN)    :: c
        real(wp),                intent(IN)    :: accumulation_rate ! [kg m-2 s-1] A_t
        real(wp),                intent(IN)    :: dt_seconds        ! [s]

        ! Local variables
        integer      :: k
        real(wp)     :: rho, temp, m_s, drho, sigma, rho_upd
        real(wp_acc) :: mass_above, p_ice_mpa, p_bubble_mpa, dP

        if (n .le. 0) return

        ! Accumulated overburden mass. wp_acc: a running sum over layers that
        ! feeds a pressure entering the mid/high branches cubed.
        mass_above = 0.0_wp_acc

        do k = 1, n

            rho  = density(k)
            temp = temperature(k)
            m_s  = mass(k)

            if (.not. ieee_is_finite(rho) .or. .not. ieee_is_finite(temp)) then
                mass_above = mass_above + max(real(m_s,wp_acc),0.0_wp_acc)
                cycle
            end if

            if (real(rho,wp_acc) .ge. real(c%rho_i,wp_acc) - TOL_TINY) then
                density(k) = c%rho_i
                mass_above = mass_above + max(real(m_s,wp_acc),0.0_wp_acc)
                cycle
            end if

            drho = 0.0_wp

            if (m_s .gt. 0.0_wp) then

                ! Overburden: g*(M_above + m/2). NOTE g = 9.81, see header.
                sigma = real(DENSIFY_GRAVITY*(mass_above + real(m_s,wp_acc)/2.0_wp_acc),wp)

                if (rho .lt. DENSIFY_RHO_LOW) then

                    ! The scheme flag selects ONLY this branch.
                    select case(c%low_density_densification)
                        case(CHION_DENSIFY_BESSI)
                            drho = bessi_low_density_rate(rho,temp,c%rho_i,accumulation_rate)
                        case(CHION_DENSIFY_HTESSEL)
                            drho = htessel_low_density_rate(c,sigma,temp,rho)
                        case DEFAULT
                            ! Chion.jl would silently fall back to BESSI here
                            ! (trap 8). chion refuses instead; the flag is
                            ! validated at parse time, so this is unreachable
                            ! unless the constants object was hand-built.
                            write(io_unit_err,*) "densify_column:: Error: &
                                                 &low_density_densification not recognized."
                            write(io_unit_err,*) "low_density_densification should be one of: &
                                                 &[CHION_DENSIFY_BESSI, CHION_DENSIFY_HTESSEL]"
                            write(io_unit_err,*) "low_density_densification = ", &
                                                 c%low_density_densification
                            stop "Program stopped."
                    end select

                else

                    p_ice_mpa    = real(sigma,wp_acc)/1.0e6_wp_acc
                    p_bubble_mpa = bubble_pressure_mpa(rho,c%rho_i)
                    dP           = p_ice_mpa - p_bubble_mpa

                    if (rho .lt. DENSIFY_RHO_MID) then
                        drho = mid_density_tendency(rho,temp,dP,c%rho_i)
                    else
                        drho = high_density_tendency(rho,temp,dP,c%rho_i)
                    end if

                end if

            end if

            rho_upd    = max(rho, rho + drho*dt_seconds)
            density(k) = min(rho_upd,c%rho_i)

            mass_above = mass_above + max(real(m_s,wp_acc),0.0_wp_acc)

        end do

        return

    end subroutine densify_column

    subroutine apply_htessel_liquid_water_compaction(mass,mass_w,density,n,mass_w_prev,c)
        ! Chion.jl/src/processes/densification.jl:52-82. Called only when the
        ! HTESSEL scheme is active, AFTER the energy solve, melt and
        ! percolation, using a snapshot of mass_w taken BEFORE the energy solve.
        !
        !     dmw = max(mass_w - mass_w_prev, 0)
        !     rho <- min(max(rho, rho*(1 + dmw/m)), rho_i)   for rho < 550, m > TOL_TINY
        !
        ! Layers at or above 550 kg m-3 are untouched, as is any layer that lost
        ! liquid water. The caller owns the snapshot; this routine never takes it.

        implicit none

        real(wp),                intent(IN)    :: mass(:)         ! (Ntot) [kg m-2] solid
        real(wp),                intent(IN)    :: mass_w(:)       ! (Ntot) [kg m-2] liquid, after
        real(wp),                intent(INOUT) :: density(:)      ! (Ntot) [kg m-3]
        integer,                 intent(IN)    :: n
        real(wp),                intent(IN)    :: mass_w_prev(:)  ! (Ntot) [kg m-2] liquid, before
        type(chion_const_class), intent(IN)    :: c

        ! Local variables
        integer  :: k
        real(wp) :: rho, m_s, dmw

        if (n .le. 0) return

        do k = 1, n

            rho = density(k)
            m_s = mass(k)

            if (rho .lt. DENSIFY_RHO_LOW .and. real(m_s,wp_acc) .gt. TOL_TINY) then

                dmw = max(mass_w(k) - mass_w_prev(k),0.0_wp)

                if (dmw .gt. 0.0_wp) then
                    density(k) = min(max(rho, rho + rho*dmw/m_s),c%rho_i)
                end if

            end if

        end do

        return

    end subroutine apply_htessel_liquid_water_compaction

end module snow_densify
