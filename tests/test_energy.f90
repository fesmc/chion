program test_energy
    ! WP5 acceptance test: snow_energy (the implicit conduction solve).
    !
    ! Coverage:
    !   1. Thomas solver vs a dense Gaussian-elimination solve, random
    !      strictly-diagonally-dominant systems at several sizes.
    !   2. Pure diffusion: a uniform profile is preserved, and total sensible
    !      energy is conserved  (Chion.jl/docs/src/tests/test_energy_flux_analytical.md).
    !   3. Steady-state conduction under a constant surface flux: the converged
    !      profile reproduces the analytic conductive flux profile.
    !   4. Energy conservation in the non-melting case.
    !   5. The melting case: surface pinned exactly at T0, residual melt
    !      energy non-negative.
    !   6. The single-layer shortcut vs a 2-layer column whose second layer is
    !      thermally negligible.
    !   7. The n <= 0 and mass(1) <= 0 early exits write nothing.
    !   8. seb_scheme = "semix": the surface row picks up SEMIX's f_sh and the
    !      ebal num_lh/denom_lh decomposition, and nothing below row 1 moves.

    use chion_defs,   only : wp, wp_acc, chion_const_class, chion_step_forcing_class, &
                             chion_const_init, CHION_SEB_BESSI, CHION_SEB_SEMIX
    use snow_energy
    use snow_seb_semix, only : semix_exchange_class, semix_snow_depth, &
                               semix_turbulent_exchange

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
        ! The assembled stencil makes the interface flux 2*G_k*(T_k - T_k+1),
        ! so the analytic temperature drops are
        !     T_k - T_(k+1) = q_k/(2*G_k) = F*(1 - k/n)*dz/K       (uniform column)
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
        real(wp) :: t_srf, dt, F, m_layer, rho, dz, K_cond, G_face
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
        K_cond = c%Ki*(rho*1.0e-3_wp)**1.88_wp
        G_face = (K_cond*dz + K_cond*dz)/((dz + dz)**2)      ! = K/(2 dz)

        worst = 0.0_wp
        do k = 1, n-1
            drop        = temperature(k) - temperature(k+1)
            drop_expect = F*(1.0_wp - real(k,wp)/real(n,wp))/(2.0_wp*G_face)
            worst       = max(worst,abs(drop-drop_expect)/abs(drop_expect))
            write(*,"(a,i2,a,g14.6,a,g14.6)") "         drop across interface ", k, &
                                              " = ", drop, "  analytic ", drop_expect
        end do

        call check("converged drops match the analytic flux profile (1e-2 rel)", &
                   worst .lt. 1.0e-2_wp, nfail)
        write(*,"(a,g14.6)") "         worst relative deviation = ", worst

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

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_srf, dt
        integer  :: k
        logical  :: all_clamped

        write(*,*)
        write(*,"(a)") "--- Melting case (two-pass re-solve) ---"

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

        call snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,0.7_wp, &
                              0.0_wp,0.0_wp,dt,res)

        call check("needs_melt is flagged", res%needs_melt, nfail)
        call check("surface is pinned EXACTLY at T0", t_srf .eq. c%T0, nfail)
        call check("temperature(1) equals t_srf", temperature(1) .eq. c%T0, nfail)
        call check("melt_energy_available >= 0", &
                   res%melt_energy_available .ge. 0.0_wp_acc, nfail)
        call check("melt_energy_available > 0 under strong forcing", &
                   res%melt_energy_available .gt. 0.0_wp_acc, nfail)
        call check("heating == energy_to_melting in the melting branch", &
                   res%heating .eq. res%energy_to_melting, nfail)

        all_clamped = .TRUE.
        do k = 1, n
            if (temperature(k) .gt. c%T0) all_clamped = .FALSE.
        end do
        call check("all layers clamped to <= T0", all_clamped, nfail)

        write(*,"(a,g16.8)") "         energy_to_melting     = ", res%energy_to_melting
        write(*,"(a,g16.8)") "         melt_energy_available = ", res%melt_energy_available

        ! Single-layer melting branch, same expectations.
        mass        = 0.0_wp
        density     = 0.0_wp
        temperature = 0.0_wp
        mass(1)        = 100.0_wp
        density(1)     = 350.0_wp
        temperature(1) = 272.0_wp
        t_srf          = 272.0_wp

        call snow_energy_flux(mass,density,temperature,t_srf,1,c,forc,0.7_wp, &
                              0.0_wp,0.0_wp,dt,res)

        call check("single layer: needs_melt flagged", res%needs_melt, nfail)
        call check("single layer: surface pinned exactly at T0", t_srf .eq. c%T0, nfail)
        call check("single layer: melt_energy_available >= 0", &
                   res%melt_energy_available .ge. 0.0_wp_acc, nfail)

        return

    end subroutine test_melting

    ! =====================================================================
    ! 6. Single-layer shortcut vs a 2-layer column with a negligible layer 2
    ! =====================================================================

    subroutine test_single_layer_shortcut(nfail)
        ! The n == 1 branch is a separate closed form, not the solver. It must
        ! agree with the solver in the limit where the second layer carries no
        ! heat capacity and starts at the same temperature.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        type(snow_energy_result_class) :: res1, res2

        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_srf1, t_srf2, dt

        write(*,*)
        write(*,"(a)") "--- Single-layer shortcut vs 2-layer limit ---"

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
                   res%energy_to_melting .eq. 0.0_wp_acc .and. &
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
    ! Helpers
    ! =====================================================================

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
    ! 8. seb_scheme = "semix": the surface row only
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

        integer,  parameter :: n = 5
        real(wp) :: mass(Ntot), density(Ntot), temperature(Ntot)
        real(wp) :: t_b(Ntot), t_s(Ntot)
        real(wp) :: t_srf_b, t_srf_s, dt, Tn, h_snow
        real(wp) :: lw_const, lw_lin, sw_abs, expect_const, expect_lin
        integer  :: k

        write(*,*)
        write(*,"(a)") "--- 8. seb_scheme = semix surface row ---"

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
        c%seb_scheme = CHION_SEB_SEMIX
        t_s(1:n)     = temperature(1:n)
        t_srf_s      = Tn

        call snow_energy_flux(mass,density,t_s,t_srf_s,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_s)

        sx = semix_turbulent_exchange(c,h_snow,forc%air_temperature,Tn, &
                                      forc%wind_speed,forc%air_pressure, &
                                      forc%relative_humidity,.TRUE.)

        ! Rebuild q_const/q_lin from the pieces: longwave and shortwave are
        ! the untouched BESSI expressions, the two turbulent terms are SEMIX's.
        lw_const = c%sigma_sb*(c%eps_air*forc%air_temperature**4 &
                               + c%eps_snow*3.0_wp*Tn**4)
        lw_lin   = c%sigma_sb*c%eps_snow*4.0_wp*Tn**3
        sw_abs   = (1.0_wp - 0.75_wp)*forc%shortwave_down

        expect_const = forc%air_temperature*sx%f_sh + lw_const + sw_abs &
                       - sx%f_lh*(sx%qsat - sx%dqsatdT*Tn - sx%q_air)
        expect_lin   = sx%f_sh + lw_lin + sx%f_lh*sx%dqsatdT

        ! check_val takes an ABSOLUTE tolerance. q_const is ~1000 W m-2 and
        ! num_lh is itself a difference of ~600 W m-2 terms, so 0.1 W m-2 is
        ! about what single precision can hold here; q_lin is ~14 W m-2 K-1
        ! with no such cancellation.
        call check_val("q_const = SH + LW + SW + ebal num_lh", &
                       res_s%surface_flux_constant, expect_const, 0.1_wp, nfail)
        call check_val("q_lin = f_sh + LW_lin + ebal denom_lh", &
                       res_s%surface_flux_linear, expect_lin, 1.0e-3_wp, nfail)

        ! The scheme must actually have changed something.
        c%seb_scheme = CHION_SEB_BESSI
        t_b(1:n)     = temperature(1:n)
        t_srf_b      = Tn

        call snow_energy_flux(mass,density,t_b,t_srf_b,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_b)

        call check("semix and bessi surface rows differ", &
                   res_s%surface_flux_linear .ne. res_b%surface_flux_linear, nfail)

        ! --- nothing below row 1 moves ----------------------------------
        ! With both turbulent fluxes prescribed, neither scheme's coefficients
        ! are consulted, so the whole solve must agree to the last bit.
        forc%has_q_sh = .TRUE.
        forc%q_sh     = -12.0_wp
        forc%has_q_lh = .TRUE.
        forc%q_lh     = -3.0_wp

        c%seb_scheme = CHION_SEB_SEMIX
        t_s(1:n)     = temperature(1:n)
        t_srf_s      = Tn
        call snow_energy_flux(mass,density,t_s,t_srf_s,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_s)

        c%seb_scheme = CHION_SEB_BESSI
        t_b(1:n)     = temperature(1:n)
        t_srf_b      = Tn
        call snow_energy_flux(mass,density,t_b,t_srf_b,n,c,forc,0.75_wp, &
                              0.0_wp,0.0_wp,dt,res_b)

        call check("prescribed fluxes make the two schemes identical", &
                   all(t_s(1:n) .eq. t_b(1:n)) .and. t_srf_s .eq. t_srf_b .and. &
                   res_s%surface_flux_constant .eq. res_b%surface_flux_constant .and. &
                   res_s%surface_flux_linear .eq. res_b%surface_flux_linear, nfail)

        return

    end subroutine test_semix_surface_row

end program test_energy
