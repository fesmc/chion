program test_energy
    ! WP5 acceptance test: snow_energy (the implicit conduction solve).
    !
    ! Coverage:
    !   1. Thomas solver vs a dense Gaussian-elimination solve, random
    !      strictly-diagonally-dominant systems at several sizes.
    !   2. Pure diffusion: a uniform profile is preserved, and total sensible
    !      energy is conserved  (Chion.jl/docs/src/tests/test_energy_flux_analytical.md).
    !   3. Steady-state conduction under a constant surface flux: the converged
    !      profile reproduces the analytic conductive flux profile, and the
    !      Robin interface temperature is the analytic T1 + F dz1/(2 K1).
    !   4. Energy conservation in the non-melting case.
    !   5. The melting case (Robin, Chion.jl 03bb445): interface pinned exactly
    !      at T0 over a subfreezing top cell, melt energy non-negative, and
    !      sum m ci dT + melt energy == Q(T0) dt (subsurface conduction kept).
    !   6. One layer goes through the same matrix: it agrees with a 2-layer
    !      column whose second layer is thermally negligible, and an insulated
    !      equilibrium is preserved down to a 1e-8 kg m-2 top cell.
    !   7. The n <= 0 and mass(1) <= 0 early exits write nothing.
    !  11. Two-layer diffusion step against its closed form (Chion.jl 03bb445).
    !   8. CLIMBER-X SEMIX (seb_scheme = "semix", turbulent_flux_scheme =
    !      "climberx"): the surface row picks up SEMIX's f_sh and the
    !      ebal num_lh/denom_lh decomposition, and nothing below row 1 moves.
    !   9. Harmonic interface conductance (Chion.jl 81034fa/d0146e1): equals
    !      the former arithmetic form on a uniform column, coefficient and full
    !      step; two layers of contrasting conductivity carry the analytic
    !      series-resistance steady flux.
    !  10. Calonne et al. (2019) conductivity against hand-computed values.
    !  12. Thermal ice substrate (Chion.jl 03bb445): a zero-size substrate is
    !      bit-identical to none; an isothermal column over the substrate is
    !      preserved; snow + substrate conserve energy with the insulated base
    !      at the bottom of the substrate; the steady flux profile runs through
    !      the substrate; bare ice cools below T0 under a negative balance.

    use chion_defs,   only : wp, wp_acc, chion_const_class, chion_step_forcing_class, &
                             chion_const_init, CHION_SEB_BESSI, CHION_SEB_SEMIX, &
                             CHION_TURB_BESSI, CHION_TURB_CLIMBERX
    use snow_energy
    use snow_seb_semix, only : semix_exchange_class, semix_flux_lin_class, &
                               semix_snow_depth, semix_turbulent_exchange, &
                               semix_surface_emissivity, semix_longwave_down, &
                               semix_longwave_linearized

    implicit none

    integer, parameter :: Ntot = 15

    integer  :: nfail

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP5 acceptance test: snow_energy"
    write(*,"(a)") "=========================================================="
    write(*,*)

    call test_thomas(nfail)
    call test_pure_diffusion(nfail)
    call test_steady_state(nfail)
    call test_energy_conservation(nfail)
    call test_melting(nfail)
    call test_single_layer_shortcut(nfail)
    call test_early_exit(nfail)
    call test_semix_surface_row(nfail)
    call test_interface_conductance(nfail)
    call test_conductivity_calonne(nfail)
    call test_two_layer_diffusion(nfail)
    call test_ice_substrate(nfail)

    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP5 (energy): ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP5 (energy): ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    ! =====================================================================
    ! 1. Thomas algorithm vs dense Gaussian elimination
    ! =====================================================================

    subroutine test_thomas(nfail)

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: nsize = 4
        integer, parameter :: sizes(nsize) = [2, 3, 5, 15]

        integer  :: is, itrial, i, j, n
        integer  :: seed
        real(wp) :: lower(Ntot), diag(Ntot), upper(Ntot), rhs(Ntot)
        real(wp) :: diag_work(Ntot), rhs_work(Ntot)
        real(wp) :: A(Ntot,Ntot), b(Ntot), x_dense(Ntot)
        real(wp) :: err, denom, worst

        write(*,"(a)") "--- Thomas solver vs dense Gaussian elimination ---"

        seed  = 20250721
        worst = 0.0_wp

        do is = 1, nsize

            n = sizes(is)

            do itrial = 1, 20

                ! Random tridiagonal system, made strictly diagonally dominant
                ! so that the pivot-free Thomas sweep is well posed -- which is
                ! exactly the situation the assembled conduction matrix is in.
                lower = 0.0_wp ; diag = 0.0_wp ; upper = 0.0_wp ; rhs = 0.0_wp

                do i = 1, n
                    if (i .lt. n) then
                        lower(i) = random_uniform(seed) - 0.5_wp
                        upper(i) = random_uniform(seed) - 0.5_wp
                    end if
                    rhs(i) = random_uniform(seed)*10.0_wp - 5.0_wp
                end do

                do i = 1, n
                    diag(i) = 1.0_wp + random_uniform(seed)
                    if (i .gt. 1) diag(i) = diag(i) + abs(lower(i-1))
                    if (i .lt. n) diag(i) = diag(i) + abs(upper(i))
                end do

                ! Dense image, using the module's indexing convention:
                ! lower(k) is the sub-diagonal entry of row k+1,
                ! upper(k) is the super-diagonal entry of row k.
                A = 0.0_wp
                do i = 1, n
                    A(i,i) = diag(i)
                    if (i .lt. n) then
                        A(i,i+1)   = upper(i)
                        A(i+1,i)   = lower(i)
                    end if
                    b(i) = rhs(i)
                end do

                call dense_solve(A,b,x_dense,n)

                diag_work(1:n) = diag(1:n)
                rhs_work(1:n)  = rhs(1:n)
                call solve_tridiagonal_thomas(lower,diag_work,upper,rhs_work,n)

                do j = 1, n
                    denom = max(abs(x_dense(j)),1.0_wp)
                    err   = abs(rhs_work(j) - x_dense(j))/denom
                    worst = max(worst,err)
                end do

            end do

        end do

        call check("max relative difference vs dense solve < 1e-5", &
                   worst .lt. 1.0e-5_wp, nfail)
        write(*,"(a,g14.6)") "         worst relative difference = ", worst

        return

    end subroutine test_thomas

    ! =====================================================================
    ! 2. Pure diffusion: uniform profile preserved, sensible energy conserved
    ! =====================================================================

    subroutine test_pure_diffusion(nfail)

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: n = 6

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp)     :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp)     :: t_srf, dt
        real(wp_acc) :: e0, e1
        integer      :: k, step
        real(wp)     :: dev

        write(*,*)
        write(*,"(a)") "--- Pure diffusion (no surface flux) ---"

        call chion_const_init(c)
        call quiet_forcing(forc)

        ! All surface fluxes prescribed and zero, and both emissivities zeroed,
        ! so the top row carries no forcing at all and the step is pure
        ! conduction. (A prescribed q_lw_down of zero is NOT enough: the
        ! outgoing eps_snow*sigma*T^4 term is always linearized in.)
        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp
        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        do k = 1, n
            mass(k)        = 100.0_wp
            density(k)     = 350.0_wp
            temperature(k) = 250.0_wp
        end do

        t_srf = 250.0_wp
        dt    = 3600.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.8_wp, &
                              0.0_wp,0.0_wp,dt,res)

        dev = 0.0_wp
        do k = 1, n
            dev = max(dev,abs(temperature(k) - 250.0_wp))
        end do

        call check("uniform profile is preserved to roundoff", &
                   dev .le. 1.0e-3_wp, nfail)
        call check("no melt flagged", .not. res%needs_melt, nfail)

        ! Non-uniform start: conduction redistributes but conserves
        ! E = sum_k m_k*ci*T_k, because the bottom is zero-flux and the
        ! surface flux is identically zero.
        do k = 1, n
            temperature(k) = 240.0_wp + 4.0_wp*real(k,wp)
        end do

        e0 = sensible_energy(mass,temperature,n,c)

        dt = 86400.0_wp

        do step = 1, 500
            call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.8_wp, &
                                  0.0_wp,0.0_wp,dt,res)
        end do

        e1 = sensible_energy(mass,temperature,n,c)

        call check("sensible energy conserved under pure diffusion (1e-5 rel)", &
                   abs(e1-e0)/abs(e0) .lt. 1.0e-5_wp_acc, nfail)
        write(*,"(a,g16.8)") "         relative energy drift = ", abs(e1-e0)/abs(e0)

        dev = 0.0_wp
        do k = 1, n
            dev = max(dev,abs(temperature(k) - temperature(1)))
        end do
        call check("profile has relaxed towards uniform", dev .lt. 0.5_wp, nfail)

        return

    end subroutine test_pure_diffusion

    ! =====================================================================
    ! 3. Steady-state conduction under a constant surface flux
    ! =====================================================================

    subroutine test_steady_state(nfail)
        ! Analytic reference. With a constant surface flux F, a zero-flux
        ! bottom and a uniform column, the column warms uniformly at
        !     r = F/(ci*M_tot)
        ! and the profile converges to the steady state in which the conductive
        ! flux crossing interface k carries exactly the heating of everything
        ! below it:
        !     q_k = F*(1 - M_k/M_tot),      M_k = sum_{j<=k} m_j.
        ! The assembled stencil makes the interface flux G_k*(T_k - T_k+1),
        ! so the analytic temperature drops are
        !     T_k - T_(k+1) = q_k/G_k = F*(1 - k/n)*dz/K           (uniform column)
        ! i.e. the conductive flux is exactly linear in depth and the drops
        ! decrease linearly with depth. This is the strongest analytic
        ! statement available for a zero-flux-bottom column: a genuinely
        ! constant gradient would require a Dirichlet base, which the scheme
        ! deliberately does not have.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: n = 5

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_srf, dt, F, m_layer, rho, dz, G_face
        real(wp) :: drop, drop_expect, worst
        integer  :: k, step

        write(*,*)
        write(*,"(a)") "--- Steady-state conduction under constant surface flux ---"

        call chion_const_init(c)
        call quiet_forcing(forc)

        ! Longwave off entirely: only then is the surface forcing exactly F.
        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp

        m_layer = 100.0_wp
        rho     = 300.0_wp
        F       = 0.3_wp                     ! [W m-2] prescribed, purely constant

        ! q_sh prescribed -> a purely constant surface flux, no linear
        ! feedback, so Q_const = F and Q_lin = 0.
        forc%has_q_sh = .TRUE.
        forc%q_sh     = F

        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        do k = 1, n
            mass(k)        = m_layer
            density(k)     = rho
            temperature(k) = 230.0_wp
        end do
        t_srf = 230.0_wp

        dt = 86400.0_wp

        do step = 1, 1200
            call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.0_wp, &
                                  0.0_wp,0.0_wp,dt,res)
        end do

        call check("no melt reached (stays well below T0)", .not. res%needs_melt, nfail)

        dz = m_layer/rho

        worst = 0.0_wp
        do k = 1, n-1
            G_face      = interface_conductance(conductivity(rho,temperature(k),c),dz, &
                                               conductivity(rho,temperature(k+1),c),dz)
            drop        = temperature(k) - temperature(k+1)
            drop_expect = F*(1.0_wp - real(k,wp)/real(n,wp))/G_face
            worst       = max(worst,abs(drop-drop_expect)/abs(drop_expect))
            write(*,"(a,i2,a,g14.6,a,g14.6)") "         drop across interface ", k, &
                                              " = ", drop, "  analytic ", drop_expect
        end do

        call check("converged drops match the analytic flux profile (1e-2 rel)", &
                   worst .lt. 1.0e-2_wp, nfail)
        write(*,"(a,g14.6)") "         worst relative deviation = ", worst

        ! Robin interface: with Q_lin = 0 the whole flux F crosses the top
        ! half cell, so Ts - T1 = F dz1/(2 K1) (K1 at the top cell's T).
        drop        = t_srf - temperature(1)
        drop_expect = F*dz/(2.0_wp*conductivity(rho,temperature(1),c))
        call check("interface temperature is T1 + F dz1/(2 K1) (1e-3 rel)", &
                   abs(drop - drop_expect)/drop_expect .lt. 1.0e-3_wp, nfail)
        write(*,"(a,g14.6,a,g14.6)") "         Ts - T1 = ", drop, "  analytic ", drop_expect

        return

    end subroutine test_steady_state

    ! =====================================================================
    ! 4. Energy conservation, non-melting
    ! =====================================================================

    subroutine test_energy_conservation(nfail)
        ! Summing the assembled rows telescopes the internal conduction terms
        ! and leaves
        !     sum_k m_k*ci*(T_k^{n+1} - T_k^n) = dt*(Q_const - Q_lin*T_1^{n+1})
        ! which is the net surface energy over the step. Checked both ways:
        ! against the reported `heating`, and against the reported surface flux
        ! coefficients.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: n = 5

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp)     :: mass(Ntot), density(Ntot), temperature(Ntot), t_old(Ntot)
        real(wp)     :: t_srf, dt
        real(wp_acc) :: de, net_surface
        integer      :: k

        write(*,*)
        write(*,"(a)") "--- Energy conservation, non-melting ---"

        call chion_const_init(c)
        call quiet_forcing(forc)

        ! A realistic, fully internal surface energy balance.
        forc%has_q_sh          = .FALSE.
        forc%has_q_lw_down     = .FALSE.
        forc%has_q_sw_net      = .FALSE.
        forc%air_temperature   = 263.15_wp
        forc%shortwave_down    = 150.0_wp

        c%eps_air  = 0.80_wp
        c%eps_snow = 0.98_wp

        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        do k = 1, n
            mass(k)        = 120.0_wp
            density(k)     = 320.0_wp + 20.0_wp*real(k,wp)
            temperature(k) = 258.0_wp - 1.0_wp*real(k,wp)
        end do

        t_old(1:n) = temperature(1:n)
        t_srf      = temperature(1)
        dt         = 3600.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res)

        call check("no melt in this configuration", .not. res%needs_melt, nfail)

        de = 0.0_wp_acc
        do k = 1, n
            de = de + real(mass(k),wp_acc)*real(c%ci,wp_acc) &
                      *(real(temperature(k),wp_acc) - real(t_old(k),wp_acc))
        end do

        net_surface = real(dt,wp_acc) &
                      *(real(res%surface_flux_constant,wp_acc) &
                        - real(res%surface_flux_linear,wp_acc)*real(t_srf,wp_acc))

        call check_rel("column sensible-energy gain == heating", &
                       de, res%heating, 1.0e-4_wp_acc, nfail)
        call check_rel("heating == net surface energy over the step", &
                       res%heating, net_surface, 1.0e-6_wp_acc, nfail)

        return

    end subroutine test_energy_conservation

    ! =====================================================================
    ! 5. Melting case
    ! =====================================================================

    subroutine test_melting(nfail)

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: n = 4

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp)     :: mass(Ntot), density(Ntot), temperature(Ntot), t_old(Ntot)
        real(wp)     :: t_srf, dt
        real(wp_acc) :: de
        integer      :: k
        logical      :: all_clamped

        write(*,*)
        write(*,"(a)") "--- Melting case (Robin interface held at T0) ---"

        call chion_const_init(c)
        call quiet_forcing(forc)

        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp

        ! A large prescribed sensible flux drives the surface past T0.
        forc%has_q_sh = .TRUE.
        forc%q_sh     = 300.0_wp

        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        do k = 1, n
            mass(k)        = 100.0_wp
            density(k)     = 350.0_wp
            temperature(k) = 272.0_wp
        end do

        t_srf = 272.0_wp
        dt    = 3600.0_wp
        t_old(1:n) = temperature(1:n)

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.7_wp, &
                              0.0_wp,0.0_wp,dt,res)

        call check("needs_melt is flagged", res%needs_melt, nfail)
        call check("interface is pinned EXACTLY at T0", t_srf .eq. c%T0, nfail)
        call check("top cell stays subfreezing (T1 < T0 = Ts)", temperature(1) .lt. c%T0, nfail)
        call check("melt_energy_available > 0 under strong forcing", &
                   res%melt_energy_available .gt. 0.0_wp_acc, nfail)
        call check_rel("heating == Q(T0) dt in the melting branch", res%heating, &
                       real(dt,wp_acc)*(real(res%surface_flux_constant,wp_acc) &
                       - real(res%surface_flux_linear,wp_acc)*real(c%T0,wp_acc)), &
                       1.0e-12_wp_acc, nfail)

        all_clamped = .TRUE.
        do k = 1, n
            if (temperature(k) .gt. c%T0) all_clamped = .FALSE.
        end do
        call check("all layers clamped to <= T0", all_clamped, nfail)

        ! Conservation: what the interface does not conduct into the column
        ! is exactly the melt energy.
        de = 0.0_wp_acc
        do k = 1, n
            de = de + real(mass(k),wp_acc)*real(c%ci,wp_acc) &
                      *(real(temperature(k),wp_acc) - real(t_old(k),wp_acc))
        end do
        call check_rel("sum m ci dT + melt energy == Q(T0) dt", &
                       de + res%melt_energy_available, res%heating, 1.0e-4_wp_acc, nfail)

        write(*,"(a,g16.8)") "         conducted into column = ", de
        write(*,"(a,g16.8)") "         melt_energy_available = ", res%melt_energy_available

        ! One layer: same expectations (it goes through the same matrix).
        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp
        mass(1)        = 100.0_wp
        density(1)     = 350.0_wp
        temperature(1) = 272.0_wp
        t_srf          = 272.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,1,c,forc,0.7_wp, &
                              0.0_wp,0.0_wp,dt,res)

        de = real(mass(1),wp_acc)*real(c%ci,wp_acc) &
             *(real(temperature(1),wp_acc) - 272.0_wp_acc)

        call check("single layer: needs_melt flagged", res%needs_melt, nfail)
        call check("single layer: interface pinned exactly at T0", t_srf .eq. c%T0, nfail)
        call check("single layer: melt_energy_available > 0", &
                   res%melt_energy_available .gt. 0.0_wp_acc, nfail)
        call check_rel("single layer: m ci dT + melt energy == Q(T0) dt", &
                       de + res%melt_energy_available, res%heating, 1.0e-4_wp_acc, nfail)

        call test_robin_melt_julia(nfail)

        return

    end subroutine test_melting

    ! =====================================================================
    ! 6. One layer through the matrix vs a 2-layer column, thin-layer limit
    ! =====================================================================

    subroutine test_single_layer_shortcut(nfail)
        ! Since Chion.jl 03bb445 there is no closed-form n == 1 branch: one
        ! layer goes through the same matrix. It must agree with a 2-layer
        ! column whose second layer carries no heat capacity and starts at the
        ! same temperature. And with no surface forcing at all an isothermal
        ! layer stays exactly put, also when it is vanishingly thin (Chion.jl
        ! test_longwave_consistency.jl: the row-1 diagonal is written so that
        ! it stays exactly 1 as m -> 0).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res1, res2

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_srf1, t_srf2, dt
        integer  :: k
        character(len=96) :: label
        real(wp), parameter :: THIN_MASSES(2) = [100.0_wp, 1.0e-8_wp]

        write(*,*)
        write(*,"(a)") "--- One layer vs 2-layer limit; thin-layer insulated equilibrium ---"

        call chion_const_init(c)
        call quiet_forcing(forc)

        forc%has_q_sw_net    = .FALSE.
        forc%has_q_lw_down   = .FALSE.
        forc%has_q_sh        = .FALSE.
        forc%air_temperature = 265.0_wp
        forc%shortwave_down  = 100.0_wp

        dt = 3600.0_wp

        ! Single layer
        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp
        mass(1)        = 100.0_wp
        density(1)     = 330.0_wp
        temperature(1) = 260.0_wp
        t_srf1         = 260.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf1,1,c,forc,0.8_wp, &
                              0.0_wp,0.0_wp,dt,res1)

        ! Two layers, the second one thermally negligible and isothermal with
        ! the first, so it neither stores nor conducts any appreciable energy.
        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp
        mass(1)        = 100.0_wp
        density(1)     = 330.0_wp
        temperature(1) = 260.0_wp
        mass(2)        = 1.0e-4_wp
        density(2)     = 330.0_wp
        temperature(2) = 260.0_wp
        t_srf2         = 260.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf2,2,c,forc,0.8_wp, &
                              0.0_wp,0.0_wp,dt,res2)

        call check_val("surface temperature agrees", t_srf2, t_srf1, 1.0e-4_wp, nfail)
        call check_rel("heating agrees", res2%heating, res1%heating, 1.0e-3_wp_acc, nfail)
        call check_val("surface_flux_constant identical", &
                       res2%surface_flux_constant, res1%surface_flux_constant, &
                       0.0_wp, nfail)
        call check_val("surface_flux_linear identical", &
                       res2%surface_flux_linear, res1%surface_flux_linear, &
                       0.0_wp, nfail)

        ! Insulated equilibrium: every flux prescribed zero, no emission.
        call chion_const_init(c)
        call quiet_forcing(forc)
        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp
        c%D_sh     = 0.0_wp

        do k = 1, 2
            mass        = 0.0_wp
            density     = 0.0_wp
            temperature = 0.0_wp
            mass(1)        = THIN_MASSES(k)
            density(1)     = 300.0_wp
            temperature(1) = 260.0_wp
            t_srf1         = 260.0_wp

            call snow_energy_flux(mass,density,temperature,t_srf1,1,c,forc,0.8_wp, &
                                  0.0_wp,0.0_wp,10800.0_wp,res1)

            write(label,"(a,es8.1,a)") "insulated equilibrium kept, m1 = ", THIN_MASSES(k), &
                                        " (T1, Ts exact; no melt)"
            call check(trim(label), temperature(1) .eq. 260.0_wp .and. t_srf1 .eq. 260.0_wp &
                       .and. res1%melt_energy_available .eq. 0.0_wp_acc, nfail)
        end do

        return

    end subroutine test_single_layer_shortcut

    ! =====================================================================
    ! 7. Early exits
    ! =====================================================================

    subroutine test_early_exit(nfail)
        ! energy_flux.jl:343-355. The guard is mass(1) <= 0, NOT
        ! TOL_EMPTY_LAYER: a surface layer with mass 1e-11 does NOT take the
        ! early exit here, even though surface_has_snow calls it empty.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_srf
        integer  :: k
        logical  :: untouched

        write(*,*)
        write(*,"(a)") "--- Early exits (n <= 0, mass(1) <= 0) ---"

        call chion_const_init(c)
        call quiet_forcing(forc)
        forc%has_q_sh = .TRUE.
        forc%q_sh     = 500.0_wp

        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        do k = 1, 3
            mass(k)        = 100.0_wp
            density(k)     = 350.0_wp
            temperature(k) = 251.0_wp + real(k,wp)
        end do
        t_srf = 999.0_wp

        ! n = 0
        call snow_energy_flux(mass,density,temperature,t_srf,0,c,forc,0.8_wp, &
                              1.5_wp,2.5_wp,3600.0_wp,res)

        untouched = (t_srf .eq. 999.0_wp)
        do k = 1, 3
            if (temperature(k) .ne. 251.0_wp + real(k,wp)) untouched = .FALSE.
        end do

        call check("n=0: temperature and t_srf untouched", untouched, nfail)
        call check("n=0: needs_melt false", .not. res%needs_melt, nfail)
        call check("n=0: all energies zero", &
                   res%melt_energy_available .eq. 0.0_wp_acc .and. &
                   res%heating .eq. 0.0_wp_acc .and. &
                   res%surface_flux_constant .eq. 0.0_wp .and. &
                   res%surface_flux_linear .eq. 0.0_wp, nfail)
        call check("n=0: latent coefficients echoed through", &
                   res%latent_heat_linear_coefficient .eq. 1.5_wp .and. &
                   res%latent_heat_constant_term .eq. 2.5_wp, nfail)

        ! mass(1) = 0 with n > 0
        mass(1) = 0.0_wp
        t_srf   = 999.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,3,c,forc,0.8_wp, &
                              1.5_wp,2.5_wp,3600.0_wp,res)

        untouched = (t_srf .eq. 999.0_wp)
        do k = 1, 3
            if (temperature(k) .ne. 251.0_wp + real(k,wp)) untouched = .FALSE.
        end do

        call check("mass(1)=0: temperature and t_srf untouched", untouched, nfail)
        call check("mass(1)=0: needs_melt false and energies zero", &
                   (.not. res%needs_melt) .and. res%heating .eq. 0.0_wp_acc, nfail)
        call check("mass(1)=0: latent coefficients echoed through", &
                   res%latent_heat_linear_coefficient .eq. 1.5_wp .and. &
                   res%latent_heat_constant_term .eq. 2.5_wp, nfail)

        ! The distinguishing case: mass(1) between 0 and TOL_EMPTY_LAYER does
        ! NOT take the early exit. If anyone "unifies" the thresholds, this
        ! check fails.
        mass(1) = 1.0e-11_wp
        t_srf   = 999.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,3,c,forc,0.8_wp, &
                              0.0_wp,0.0_wp,3600.0_wp,res)

        call check("mass(1)=1e-11 does NOT take the early exit (guard is <= 0)", &
                   t_srf .ne. 999.0_wp, nfail)

        return

    end subroutine test_early_exit

    ! =====================================================================
    ! 9. Harmonic interface conductance
    ! =====================================================================

    subroutine test_interface_conductance(nfail)
        ! (a) On a uniform column (K, dz) the harmonic conductance with
        !     beta = -dt/ci is the former arithmetic form with beta = -2 dt/ci:
        !     G = K/dz = 2*(K dz + K dz)/(2 dz)^2. Checked on the coefficient
        !     and on a full forced step against the old operator assembled
        !     here and solved densely.
        ! (b) Two layers of very different conductivity under a constant
        !     surface flux F and a zero-flux bottom: in the quasi-steady state
        !     the interface carries q = F*m2/(m1+m2), and the two half-layer
        !     resistances in series give
        !         T1 - T2 = q*(dz1/(2 K1) + dz2/(2 K2)).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer, parameter :: n = 4

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: A(n,n), b(n), t_ref(n)
        real(wp) :: K_u, dz, G_old, G_new, beta, lambda, F, t_srf, dt, worst
        real(wp) :: K1, K2, dz1, dz2, q, drop, drop_expect
        integer  :: i, j, k, step

        real(wp), parameter :: Ks(3)  = [0.05_wp, 0.3_wp, 2.1_wp]
        real(wp), parameter :: dzs(3) = [0.02_wp, 0.4_wp, 3.0_wp]

        write(*,*)
        write(*,"(a)") "--- Harmonic interface conductance ---"

        ! (a) coefficient
        worst = 0.0_wp
        do i = 1, size(Ks)
            do j = 1, size(dzs)
                G_old = 2.0_wp*(Ks(i)*dzs(j) + Ks(i)*dzs(j))/((dzs(j) + dzs(j))**2)
                G_new = interface_conductance(Ks(i),dzs(j),Ks(i),dzs(j))
                worst = max(worst,abs(G_new - G_old)/G_old)
            end do
        end do
        call check("uniform: harmonic G == 2 x arithmetic G to round-off (4 eps)", &
                   worst .le. 4.0_wp*epsilon(1.0_wp), nfail)
        write(*,"(a,g14.6)") "         worst relative difference = ", worst

        ! (a) one forced step on a uniform column, against the old operator.
        ! With Q_lin = 0 the Robin boundary delivers exactly F to row 1, the
        ! old surface term, as long as the interface Ts = T1 + F dz/(2K) stays
        ! below T0: F = 5 W m-2 keeps it ~5 K above the cell centre.
        call chion_const_init(c)
        call quiet_forcing(forc)
        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp
        F          = 5.0_wp
        forc%q_sh  = F

        mass = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        mass(1:n)        = 100.0_wp
        density(1:n)     = 350.0_wp
        temperature(1:n) = 250.0_wp
        t_srf = 250.0_wp
        dt    = 86400.0_wp

        K_u    = conductivity(350.0_wp,250.0_wp,c)
        dz     = 100.0_wp/350.0_wp
        G_old  = (K_u*dz + K_u*dz)/((dz + dz)**2)
        beta   = -2.0_wp*dt/c%ci/100.0_wp
        lambda = dt/c%ci/100.0_wp

        A = 0.0_wp
        do k = 1, n-1
            A(k,k+1) = beta*G_old
            A(k+1,k) = beta*G_old
        end do
        do k = 1, n
            A(k,k) = 1.0_wp - sum(A(k,:))
        end do
        b    = 250.0_wp
        b(1) = b(1) + lambda*F
        call dense_solve(A,b,t_ref,n)

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.0_wp, &
                              0.0_wp,0.0_wp,dt,res)

        call check("uniform column: interface stays below T0 (no melt)", &
                   .not. res%needs_melt, nfail)
        worst = maxval(abs(temperature(1:n) - t_ref))
        call check("uniform column: one step equals the arithmetic-form operator", &
                   worst .le. 64.0_wp*epsilon(1.0_wp)*250.0_wp, nfail)
        write(*,"(a,g14.6,a,g14.6)") "         max |dT| = ", worst, &
                                     "   surface warming = ", t_ref(1) - 250.0_wp

        ! (b) two layers, K contrast ~8
        call chion_const_init(c)
        call quiet_forcing(forc)
        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp
        F          = 0.3_wp
        forc%q_sh  = F

        mass = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        mass(1:2)        = 250.0_wp
        density(1)       = 200.0_wp
        density(2)       = 600.0_wp
        temperature(1:2) = 230.0_wp
        t_srf = 230.0_wp

        do step = 1, 1200
            call snow_energy_flux(mass,density,temperature,t_srf,2,c,forc,0.0_wp, &
                                  0.0_wp,0.0_wp,dt,res)
        end do
        call check("two layers: no melt reached", .not. res%needs_melt, nfail)

        K1  = conductivity(density(1),temperature(1),c)
        K2  = conductivity(density(2),temperature(2),c)
        dz1 = mass(1)/density(1)
        dz2 = mass(2)/density(2)

        G_new = interface_conductance(K1,dz1,K2,dz2)
        call check("two layers: G == 1/(dz1/(2 K1) + dz2/(2 K2)) to round-off (4 eps)", &
                   abs(G_new*(dz1/(2.0_wp*K1) + dz2/(2.0_wp*K2)) - 1.0_wp) &
                   .le. 4.0_wp*epsilon(1.0_wp), nfail)

        q           = F*mass(2)/(mass(1) + mass(2))
        drop        = temperature(1) - temperature(2)
        drop_expect = q*(dz1/(2.0_wp*K1) + dz2/(2.0_wp*K2))
        call check("two layers: steady drop is the series-resistance flux (1e-2 rel)", &
                   abs(drop - drop_expect)/drop_expect .lt. 1.0e-2_wp, nfail)
        write(*,"(a,g14.6,a,g14.6,a,g14.6)") "         drop = ", drop, "  analytic ", &
                                             drop_expect, "  K2/K1 = ", K2/K1

        return

    end subroutine test_interface_conductance

    ! =====================================================================
    ! 10. Calonne et al. (2019) conductivity
    ! =====================================================================

    subroutine test_conductivity_calonne(nfail)
        ! Reference values computed by hand (double precision) from Eq. (5) of
        ! Calonne et al. (2019) with rho_i = 917 kg m-3:
        !   theta  = 1/(1 + exp(-0.04 (rho - 450)))
        !   k_ice  = 9.828 exp(-5.7e-3 T),  k_air = 2.334e-3 T^1.5/(164.54 + T)
        !   k_snow = 0.024 - 1.23e-4 rho + 2.5e-6 rho^2
        !   k_firn = 2.107 + 3.618e-3 (rho - rho_i)
        !   K = (1-theta) k_ice k_air/(2.107*0.024) k_snow + theta k_ice/2.107 k_firn
        ! (300,253): theta 0.00247262, k_snow 0.2121           -> 0.2183550623
        ! (500,253): theta 0.880797,   k_firn 0.598294         -> 0.6535479534
        ! (917,263): theta ~1, k_firn = 2.107, k_ice 2.1949    -> 2.194897732
        ! The former Ki*(rho/1000)^1.88 gives 0.2184, 0.5705 and 1.784.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        real(wp), parameter :: rho(3)  = [300.0_wp, 500.0_wp, 917.0_wp]
        real(wp), parameter :: T(3)    = [253.0_wp, 253.0_wp, 263.0_wp]
        real(wp), parameter :: Kref(3) = [0.2183550623_wp, 0.6535479534_wp, 2.194897732_wp]
        real(wp) :: K
        integer  :: i
        character(len=64) :: label

        write(*,*)
        write(*,"(a)") "--- Calonne et al. (2019) conductivity ---"

        do i = 1, size(rho)
            K = snow_thermal_conductivity(rho(i),T(i),917.0_wp)
            write(label,"(a,f5.0,a,f5.0,a)") "K(rho=", rho(i), ", T=", T(i), ") (1e-5 rel)"
            call check_val(trim(label),K,Kref(i),1.0e-5_wp*Kref(i),nfail)
        end do

        call check("K rises with density at fixed T", &
                   snow_thermal_conductivity(400.0_wp,253.0_wp,917.0_wp) .lt. &
                   snow_thermal_conductivity(600.0_wp,253.0_wp,917.0_wp), nfail)
        call check("K of dense firn falls with temperature (ice phase)", &
                   snow_thermal_conductivity(800.0_wp,263.0_wp,917.0_wp) .lt. &
                   snow_thermal_conductivity(800.0_wp,233.0_wp,917.0_wp), nfail)

        return

    end subroutine test_conductivity_calonne

    ! =====================================================================
    ! 5b. Chion.jl 03bb445 Robin melt tests, ported
    ! =====================================================================

    subroutine test_robin_melt_julia(nfail)
        ! test_case_api.jl "Robin surface boundary melts over a subfreezing
        ! top-cell centre" and "Robin melting boundary retains subsurface
        ! conduction": no emission, no sensible exchange, LW down, q_sh and
        ! q_lh prescribed zero.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp)     :: mass(Ntot), density(Ntot), temperature(Ntot), t_old(Ntot)
        real(wp)     :: t_srf, dt, input_flux
        real(wp_acc) :: de

        call chion_const_init(c)
        call quiet_forcing(forc)
        c%eps_snow = 0.0_wp
        c%D_sh     = 0.0_wp

        ! (a) one layer, 400 W m-2 of shortwave at albedo 0.81, one day.
        forc%air_temperature = 280.0_wp
        forc%has_q_sw_net    = .FALSE.
        forc%shortwave_down  = 400.0_wp

        mass = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        mass(1) = 100.0_wp; density(1) = 300.0_wp; temperature(1) = 270.0_wp
        t_srf   = 270.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,1,c,forc,0.81_wp, &
                              0.0_wp,0.0_wp,86400.0_wp,res)

        call check("Julia (a): Tsrf == T0", t_srf .eq. c%T0, nfail)
        call check("Julia (a): top-cell centre stays below T0", temperature(1) .lt. c%T0, nfail)
        call check("Julia (a): melt energy > 0", res%melt_energy_available .gt. 0.0_wp_acc, nfail)

        ! (b) two layers, 100 W m-2 net shortwave prescribed, one day.
        forc%has_q_sw_net = .TRUE.
        input_flux        = 100.0_wp
        forc%q_sw_net     = input_flux
        dt                = 86400.0_wp

        mass = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        mass(1:2)        = [100.0_wp, 300.0_wp]
        density(1:2)     = [300.0_wp, 500.0_wp]
        temperature(1:2) = [270.0_wp, 260.0_wp]
        t_srf            = 270.0_wp
        t_old(1:2)       = temperature(1:2)

        call snow_energy_flux(mass,density,temperature,t_srf,2,c,forc,0.81_wp, &
                              0.0_wp,0.0_wp,dt,res)

        de = sum(real(mass(1:2),wp_acc)*real(c%ci,wp_acc) &
                 *(real(temperature(1:2),wp_acc) - real(t_old(1:2),wp_acc)))

        call check("Julia (b): Tsrf == T0", t_srf .eq. c%T0, nfail)
        call check("Julia (b): top-cell centre stays below T0", temperature(1) .lt. c%T0, nfail)
        call check_rel("Julia (b): sum m ci dT + melt energy == input flux * dt", &
                       de + res%melt_energy_available, &
                       real(input_flux,wp_acc)*real(dt,wp_acc), 1.0e-4_wp_acc, nfail)

        return

    end subroutine test_robin_melt_julia

    ! =====================================================================
    ! 11. Two-layer diffusion step, closed form
    ! =====================================================================

    subroutine test_two_layer_diffusion(nfail)
        ! test_case_api.jl "Energy-flux two-layer diffusion uses physical
        ! interface conductance": with no surface flux at all (Q = 0, Q_lin =
        ! 0) the Robin boundary carries nothing, and two equal layers (m = 300,
        ! rho = 300, dz = 1 m) relax as
        !     T1' = T1 - a dT/(1 + 2a),  T2' = T2 + a dT/(1 + 2a),
        !     a = dt G/(ci m),  G = interface_conductance(K1,1,K2,1).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_srf, dt, G, a, e1, e2

        write(*,*)
        write(*,"(a)") "--- 11. Two-layer diffusion step, closed form ---"

        call chion_const_init(c)
        call quiet_forcing(forc)
        c%eps_snow = 0.0_wp
        c%D_sh     = 0.0_wp
        forc%air_temperature = 270.0_wp

        mass = 0.0_wp; density = 0.0_wp; temperature = 0.0_wp
        mass(1:2)        = 300.0_wp
        density(1:2)     = 300.0_wp
        temperature(1:2) = [270.0_wp, 260.0_wp]
        t_srf            = 270.0_wp
        dt               = 3600.0_wp

        G = interface_conductance(conductivity(300.0_wp,270.0_wp,c),1.0_wp, &
                                  conductivity(300.0_wp,260.0_wp,c),1.0_wp)
        a = dt*G/(c%ci*300.0_wp)
        e1 = 270.0_wp - a*10.0_wp/(1.0_wp + 2.0_wp*a)
        e2 = 260.0_wp + a*10.0_wp/(1.0_wp + 2.0_wp*a)

        call snow_energy_flux(mass,density,temperature,t_srf,2,c,forc,0.8_wp, &
                              0.0_wp,0.0_wp,dt,res)

        call check_val("layer 1 matches the closed form", temperature(1), e1, &
                       64.0_wp*epsilon(1.0_wp)*270.0_wp, nfail)
        call check_val("layer 2 matches the closed form", temperature(2), e2, &
                       64.0_wp*epsilon(1.0_wp)*270.0_wp, nfail)

        return

    end subroutine test_two_layer_diffusion

    ! =====================================================================
    ! 12. Thermal ice substrate
    ! =====================================================================

    subroutine test_ice_substrate(nfail)

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer,  parameter :: n = 3, n_ice = 5
        real(wp), parameter :: h0 = 0.05_wp

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res, res_ref

        real(wp)     :: mass(Ntot), density(Ntot), temperature(Ntot), t_ref(Ntot), t_old(Ntot)
        real(wp)     :: ice_t(n_ice), ice_old(n_ice), no_ice(0)
        real(wp)     :: row_m(n+n_ice), row_rho(n+n_ice), row_t(n+n_ice)
        real(wp)     :: t_srf, t_srf_ref, dt, F, G_face, drop, drop_expect, worst
        real(wp_acc) :: de, de_ice, m_tot, m_above
        integer      :: k, step

        write(*,*)
        write(*,"(a)") "--- Thermal ice substrate ---"

        call chion_const_init(c)

        ! --- (a) zero-size substrate == no substrate, bit for bit ---------
        call quiet_forcing(forc)
        forc%has_q_lw_down = .FALSE.
        forc%has_q_sw_net  = .FALSE.
        forc%has_q_sh      = .FALSE.
        forc%air_temperature = 263.15_wp
        forc%shortwave_down  = 150.0_wp

        call snow_column(mass,density,temperature)
        t_ref     = temperature
        t_srf     = 257.0_wp
        t_srf_ref = t_srf
        dt        = 3600.0_wp

        call snow_energy_flux(mass,density,t_ref,t_srf_ref,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_ref)
        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res,no_ice,h0)

        call check("zero-size substrate: temperatures bit-identical to none", &
                   all(temperature(1:n) .eq. t_ref(1:n)), nfail)
        call check("zero-size substrate: Ts and heating bit-identical to none", &
                   t_srf .eq. t_srf_ref .and. res%heating .eq. res_ref%heating, nfail)

        ! --- (b) isothermal column over the substrate is preserved --------
        ! No surface flux at all (quiet forcing, longwave off): an isothermal
        ! snow + ice column is a steady state of the insulated system.
        call quiet_forcing(forc)
        c%eps_air  = 0.0_wp
        c%eps_snow = 0.0_wp

        call snow_column(mass,density,temperature)
        temperature(1:n) = 255.0_wp
        ice_t = 255.0_wp
        t_srf = 255.0_wp

        do step = 1, 30
            call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.0_wp, &
                                  0.0_wp,0.0_wp,86400.0_wp,res,ice_t,h0)
        end do

        call check("isothermal snow + substrate preserved (1e-4 K)", &
                   maxval(abs(temperature(1:n) - 255.0_wp)) .lt. 1.0e-4_wp .and. &
                   maxval(abs(ice_t - 255.0_wp)) .lt. 1.0e-4_wp, nfail)

        ! --- (c) energy conservation, insulated base below the substrate --
        ! A constant surface flux F into a cold, non-uniform column over a
        ! day: every joule stays in the snow + substrate, sum m ci dT == F dt.
        ! F is small because the top snow cell is a poor conductor (Ts - T1 =
        ! F/Gs, Gs ~ 1 W m-2 K-1), and the run must not melt. The tolerance
        ! is set by sp storage of ~250 K in 1.4 t m-2 of ice (ulp*m*ci ~ 25 J
        ! per substrate layer against F dt = 432 kJ).
        F = 5.0_wp
        forc%q_sh = F

        call snow_column(mass,density,temperature)
        do k = 1, n_ice
            ice_t(k) = 250.0_wp - 2.0_wp*real(k,wp)
        end do
        t_old   = temperature
        ice_old = ice_t
        t_srf   = temperature(1)
        dt      = 86400.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.0_wp, &
                              0.0_wp,0.0_wp,dt,res,ice_t,h0)

        de_ice = 0.0_wp_acc
        do k = 1, n_ice
            de_ice = de_ice + real(c%rho_i*h0*2.0_wp**(k-1),wp_acc)*real(c%ci,wp_acc) &
                              *(real(ice_t(k),wp_acc) - real(ice_old(k),wp_acc))
        end do
        de = de_ice
        do k = 1, n
            de = de + real(mass(k),wp_acc)*real(c%ci,wp_acc) &
                      *(real(temperature(k),wp_acc) - real(t_old(k),wp_acc))
        end do

        call check("substrate: no melt", .not. res%needs_melt, nfail)
        call check("substrate gained energy (the snow base conducts into it)", &
                   de_ice .gt. 0.0_wp_acc, nfail)
        call check_rel("snow + substrate: sum m ci dT == F dt (insulated base)", &
                       de, real(F,wp_acc)*real(dt,wp_acc), 5.0e-4_wp_acc, nfail)

        ! --- (d) steady flux profile through the substrate ----------------
        ! As test 3: under a constant F the column warms uniformly and the
        ! flux crossing row interface k carries the heating of all rows
        ! below it, q_k = F*(1 - M_k/M_tot), down to zero at the insulated
        ! base of the substrate.
        F = 0.3_wp
        forc%q_sh = F

        call snow_column(mass,density,temperature)
        temperature(1:n) = 230.0_wp
        ice_t = 230.0_wp
        t_srf = 230.0_wp

        do step = 1, 2000
            call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.0_wp, &
                                  0.0_wp,0.0_wp,86400.0_wp,res,ice_t,h0)
        end do

        do k = 1, n
            row_m(k) = mass(k) ; row_rho(k) = density(k) ; row_t(k) = temperature(k)
        end do
        do k = 1, n_ice
            row_m(n+k) = c%rho_i*h0*2.0_wp**(k-1) ; row_rho(n+k) = c%rho_i ; row_t(n+k) = ice_t(k)
        end do

        m_tot = 0.0_wp_acc
        do k = 1, n+n_ice
            m_tot = m_tot + real(row_m(k),wp_acc)
        end do

        worst   = 0.0_wp
        m_above = 0.0_wp_acc
        do k = 1, n+n_ice-1
            m_above     = m_above + real(row_m(k),wp_acc)
            G_face      = interface_conductance(conductivity(row_rho(k),row_t(k),c), &
                                                row_m(k)/row_rho(k), &
                                                conductivity(row_rho(k+1),row_t(k+1),c), &
                                                row_m(k+1)/row_rho(k+1))
            drop        = row_t(k) - row_t(k+1)
            drop_expect = F*real(1.0_wp_acc - m_above/m_tot,wp)/G_face
            worst       = max(worst,abs(drop - drop_expect)/abs(drop_expect))
        end do

        call check("steady drops through snow and substrate match the flux profile (1e-2 rel)", &
                   worst .lt. 1.0e-2_wp, nfail)
        write(*,"(a,g14.6)") "         worst relative deviation = ", worst

        ! --- (e) bare ice cools below T0 under a negative balance ---------
        ! Chion.jl test_case_api "Public energy solve couples the ice
        ! substrate on bare ice": no snow, ice at T0, only longwave emission.
        call chion_const_init(c)
        call quiet_forcing(forc)
        forc%air_temperature = 250.0_wp

        ice_t = c%T0
        t_srf = c%T0

        call snow_energy_flux(mass,density,temperature,t_srf,0,c,forc,0.0_wp, &
                              0.0_wp,0.0_wp,86400.0_wp,res,ice_t,h0)

        call check("bare substrate: no melt under emission only", .not. res%needs_melt, nfail)
        call check("bare substrate: Ts < T0", t_srf .lt. c%T0, nfail)
        call check("bare substrate: top ice layer cooled below T0", ice_t(1) .lt. c%T0, nfail)

        return

    end subroutine test_ice_substrate

    subroutine snow_column(mass,density,temperature)
        ! Three snow layers of contrasting mass and density, cold and
        ! stratified, for the substrate tests.

        implicit none

        real(wp), intent(OUT) :: mass(:), density(:), temperature(:)

        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        mass(1:3)        = [60.0_wp, 150.0_wp, 300.0_wp]
        density(1:3)     = [300.0_wp, 420.0_wp, 600.0_wp]
        temperature(1:3) = [256.0_wp, 254.0_wp, 252.0_wp]

        return

    end subroutine snow_column

    ! =====================================================================
    ! Helpers
    ! =====================================================================

    function conductivity(rho,T,c) result(K)
        ! The layer conductivity snow_energy_flux uses for a layer of density
        ! rho at temperature T.

        implicit none

        real(wp),                intent(IN) :: rho, T
        type(chion_const_class), intent(IN) :: c
        real(wp) :: K

        K = snow_thermal_conductivity(rho,T,c%rho_i)

        return

    end function conductivity

    subroutine quiet_forcing(forc)
        ! A forcing with every optional flux prescribed and zero, and no
        ! humidity: the surface energy balance then contributes nothing.

        implicit none

        type(chion_step_forcing_class), intent(OUT) :: forc

        forc%air_temperature  = 260.0_wp
        forc%dt_days          = 1.0_wp
        forc%snowfall_rate    = 0.0_wp
        forc%rainfall_rate    = 0.0_wp
        forc%shortwave_down   = 0.0_wp
        forc%wind_speed       = 0.0_wp

        forc%q_sw_net  = 0.0_wp
        forc%q_lw_down = 0.0_wp
        forc%q_sh      = 0.0_wp
        forc%q_lh      = 0.0_wp

        forc%has_q_sw_net  = .TRUE.
        forc%has_q_lw_down = .TRUE.
        forc%has_q_sh      = .TRUE.
        forc%has_q_lh      = .TRUE.

        forc%relative_humidity     = 0.0_wp
        forc%has_relative_humidity = .FALSE.

        forc%air_pressure          = 101325.0_wp
        forc%prescribed_albedo     = 0.0_wp
        forc%has_prescribed_albedo = .FALSE.

        forc%latitude_deg        = 70.0_wp
        forc%day_of_year         = 1.0_wp
        forc%solar_longitude_deg = 0.0_wp

        return

    end subroutine quiet_forcing

    function sensible_energy(mass,temperature,n,c) result(e)

        implicit none

        real(wp),                intent(IN) :: mass(:)
        real(wp),                intent(IN) :: temperature(:)
        integer,                 intent(IN) :: n
        type(chion_const_class), intent(IN) :: c
        real(wp_acc) :: e

        ! Local variables
        integer :: k

        e = 0.0_wp_acc
        do k = 1, n
            e = e + real(mass(k),wp_acc)*real(c%ci,wp_acc)*real(temperature(k),wp_acc)
        end do

        return

    end function sensible_energy

    function random_uniform(seed) result(r)
        ! Reproducible linear congruential generator, so the Thomas comparison
        ! is deterministic across machines and compilers.

        implicit none

        integer, intent(INOUT) :: seed
        real(wp) :: r

        seed = mod(1103515245*seed + 12345, 2147483647)
        if (seed .lt. 0) seed = seed + 2147483647

        r = real(seed,wp)/2147483647.0_wp

        return

    end function random_uniform

    subroutine dense_solve(A,b,x,n)
        ! Gaussian elimination with partial pivoting. Independent reference
        ! implementation for the Thomas comparison.

        implicit none

        real(wp), intent(IN)  :: A(:,:)
        real(wp), intent(IN)  :: b(:)
        real(wp), intent(OUT) :: x(:)
        integer,  intent(IN)  :: n

        ! Local variables
        integer  :: i, j, k, ipiv
        real(wp) :: M(n,n+1), tmp(n+1), f, s

        M = 0.0_wp
        do i = 1, n
            do j = 1, n
                M(i,j) = A(i,j)
            end do
            M(i,n+1) = b(i)
        end do

        do k = 1, n-1
            ipiv = k
            do i = k+1, n
                if (abs(M(i,k)) .gt. abs(M(ipiv,k))) ipiv = i
            end do
            if (ipiv .ne. k) then
                tmp        = M(k,:)
                M(k,:)     = M(ipiv,:)
                M(ipiv,:)  = tmp
            end if
            do i = k+1, n
                f = M(i,k)/M(k,k)
                do j = k, n+1
                    M(i,j) = M(i,j) - f*M(k,j)
                end do
            end do
        end do

        x = 0.0_wp
        do i = n, 1, -1
            s = M(i,n+1)
            do j = i+1, n
                s = s - M(i,j)*x(j)
            end do
            x(i) = s/M(i,i)
        end do

        return

    end subroutine dense_solve

    subroutine check(label,condition,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        logical,          intent(IN)    :: condition
        integer,          intent(INOUT) :: nfail

        if (condition) then
            write(*,"(a,a)") "  ok   : ", trim(label)
        else
            write(*,"(a,a)") "  FAIL : ", trim(label)
            nfail = nfail + 1
        end if

        return

    end subroutine check

    subroutine check_val(label,value,expected,tol,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value
        real(wp),         intent(IN)    :: expected
        real(wp),         intent(IN)    :: tol
        integer,          intent(INOUT) :: nfail

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_val

    subroutine check_rel(label,value,expected,tol,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        real(wp_acc),     intent(IN)    :: tol
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp_acc) :: err

        err = abs(value-expected)/max(abs(expected),1.0e-30_wp_acc)

        if (err .le. tol) then
            write(*,"(a,a,a,g16.8)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_rel

    ! =====================================================================
    ! 8. CLIMBER-X SEMIX: the surface row only
    ! =====================================================================

    subroutine test_semix_surface_row(nfail)
        ! The SEMIX scheme swaps the flux formulas feeding rhs(1)/diag(1) and
        ! nothing else (docs/semix_port_scope.md, coupling decision alpha).
        ! Two things are checked: that q_const/q_lin pick up exactly the
        ! aerodynamic sensible term plus ebal's num_lh/denom_lh, and that the
        ! conduction machinery below row 1 is untouched -- verified by running
        ! the same column under both schemes with the turbulent terms
        ! PRESCRIBED, which must give bit-identical answers.

        implicit none

        integer, intent(INOUT) :: nfail

        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res_b, res_s
        type(semix_exchange_class)     :: sx
        type(semix_flux_lin_class)     :: lwc

        integer,  parameter :: n = 5
        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_b(Ntot), t_s(Ntot)
        real(wp) :: t_srf_b, t_srf_s, dt, Tn, h_snow
        real(wp) :: lw_const, lw_lin, sw_abs, expect_const, expect_lin
        real(wp) :: tol_wp
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- 8. seb_scheme = semix, turbulent_flux_scheme = climberx surface row ---"

        call chion_const_init(c)

        forc%air_temperature       = 263.15_wp
        forc%dt_days               = 1.0_wp/24.0_wp
        forc%snowfall_rate         = 0.0_wp
        forc%rainfall_rate         = 0.0_wp
        forc%shortwave_down        = 150.0_wp
        forc%wind_speed            = 5.0_wp
        forc%relative_humidity     = 0.75_wp
        forc%has_relative_humidity = .TRUE.
        forc%air_pressure          = 101325.0_wp

        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp

        do k = 1, n
            mass(k)        = 120.0_wp
            density(k)     = 320.0_wp + 20.0_wp*real(k,wp)
            temperature(k) = 258.0_wp - 1.0_wp*real(k,wp)
        end do

        Tn     = temperature(1)
        h_snow = semix_snow_depth(mass,density,n)
        dt     = 3600.0_wp

        ! --- the surface row under semix --------------------------------
        c%seb_scheme            = CHION_SEB_SEMIX
        c%turbulent_flux_scheme = CHION_TURB_CLIMBERX
        t_s(1:n)     = temperature(1:n)
        t_srf_s      = Tn

        call snow_energy_flux(mass,density,t_s,t_srf_s,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_s)

        sx = semix_turbulent_exchange(c,h_snow,forc%air_temperature,Tn, &
                                      forc%wind_speed,forc%air_pressure, &
                                      forc%relative_humidity,.TRUE.)

        ! Rebuild q_const/q_lin from the pieces. Shortwave is the untouched
        ! BESSI expression; longwave and both turbulent terms are SEMIX's.
        lwc = semix_longwave_linearized(c,semix_surface_emissivity(c,.TRUE.), &
                                        semix_longwave_down(c,forc%q_lw_down, &
                                                            forc%has_q_lw_down, &
                                                            forc%air_temperature),Tn)
        lw_const = lwc%constant
        lw_lin   = lwc%linear
        sw_abs   = (1.0_wp - 0.75_wp)*forc%shortwave_down

        expect_const = forc%air_temperature*sx%f_sh + lw_const + sw_abs &
                       - sx%f_lh*(sx%qsat - sx%dqsatdT*Tn - sx%q_air)
        expect_lin   = sx%f_sh + lw_lin + sx%f_lh*sx%dqsatdT

        ! The longwave the semix branch builds must absorb LESS than BESSI's,
        ! which takes the downwelling flux at face value.
        call check("semix LW constant < bessi LW constant", &
                   lw_const .lt. c%sigma_sb*(c%eps_air*forc%air_temperature**4 &
                                             + c%eps_snow*3.0_wp*Tn**4), nfail)

        ! check_val takes an ABSOLUTE tolerance. q_const is ~1000 W m-2 and
        ! num_lh is itself a difference of ~600 W m-2 terms, so 0.1 W m-2 is
        ! about what single precision can hold here; q_lin is ~14 W m-2 K-1
        ! with no such cancellation.
        call check_val("q_const = SH + LW + SW + ebal num_lh", &
                       res_s%surface_flux_constant, expect_const, 0.1_wp, nfail)
        call check_val("q_lin = f_sh + LW_lin + ebal denom_lh", &
                       res_s%surface_flux_linear, expect_lin, 1.0e-3_wp, nfail)

        ! The scheme must actually have changed something.
        c%seb_scheme            = CHION_SEB_BESSI
        c%turbulent_flux_scheme = CHION_TURB_BESSI
        t_b(1:n)     = temperature(1:n)
        t_srf_b      = Tn

        call snow_energy_flux(mass,density,t_b,t_srf_b,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_b)

        call check("semix and bessi surface rows differ", &
                   res_s%surface_flux_linear .ne. res_b%surface_flux_linear, nfail)

        ! --- nothing below row 1 moves ----------------------------------
        ! With both turbulent fluxes prescribed AND the emissivities equal,
        ! neither scheme's coefficients are consulted and the longwave terms
        ! coincide, so the whole solve must agree to round-off of wp.
        !
        ! eps_snow = 1 is what makes the longwave agree: it is the only value
        ! at which absorbing the downwelling flux with emissivity (semix) and
        ! absorbing it in full (bessi) are the same operation. Not to the bit:
        ! the two branches evaluate it in a different order
        ! (sigma*(eps_air*Ta^4 + 3*Tn^4) vs sigma*eps_air*Ta^4 + 3*sigma*Tn^4),
        ! which at sp may round one ulp apart (ifx -O2 -fp-model precise does).
        forc%has_q_sh = .TRUE.
        forc%q_sh     = -12.0_wp
        forc%has_q_lh = .TRUE.
        forc%q_lh     = -3.0_wp
        c%eps_snow    = 1.0_wp

        c%seb_scheme            = CHION_SEB_SEMIX
        c%turbulent_flux_scheme = CHION_TURB_CLIMBERX
        t_s(1:n)     = temperature(1:n)
        t_srf_s      = Tn
        call snow_energy_flux(mass,density,t_s,t_srf_s,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_s)

        c%seb_scheme            = CHION_SEB_BESSI
        c%turbulent_flux_scheme = CHION_TURB_BESSI
        t_b(1:n)     = temperature(1:n)
        t_srf_b      = Tn
        call snow_energy_flux(mass,density,t_b,t_srf_b,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_b)

        tol_wp = 4.0_wp*epsilon(1.0_wp)
        call check_val("prescribed fluxes: same q_const", res_s%surface_flux_constant, &
                       res_b%surface_flux_constant, &
                       tol_wp*abs(res_b%surface_flux_constant), nfail)
        call check_val("prescribed fluxes: same q_lin", res_s%surface_flux_linear, &
                       res_b%surface_flux_linear, &
                       tol_wp*abs(res_b%surface_flux_linear), nfail)
        call check_val("prescribed fluxes: same t_srf", t_srf_s, t_srf_b, &
                       tol_wp*t_srf_b, nfail)
        call check_val("prescribed fluxes: same column temperatures", &
                       maxval(abs(t_s(1:n)-t_b(1:n))), 0.0_wp, &
                       tol_wp*maxval(t_b(1:n)), nfail)

        return

    end subroutine test_semix_surface_row

end program test_energy
