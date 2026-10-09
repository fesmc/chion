module snow_energy
    ! Implicit (backward-Euler) snow temperature solve: a Robin surface
    ! boundary (linearized surface energy balance closed against the top
    ! cell's half-thickness conductance), vertical heat conduction below,
    ! zero-flux bottom, solved with the Thomas algorithm, plus a melting
    ! re-solve with the interface held at T0 (Chion.jl 03bb445).
    !
    ! Port of Chion.jl/src/processes/energy_flux.jl:
    !   _go_energy_flux_resolved!            -> snow_energy_flux
    !   _thomas_forward! / _thomas_backward! -> solve_tridiagonal_thomas
    !   _snow_thermal_conductivity           -> snow_thermal_conductivity
    !   interface_conductance                -> interface_conductance
    !   _clamp_to_melt!                      -> inlined (clamp loop)
    !   shortwave_absorbed                   -> inlined (step 2; see note below)
    ! Prose reference: Chion.jl/docs/src/processes/energy.md.
    !
    ! This is the only implicit solve in the model. Order of operations matters
    ! more here than anywhere else in the port, so the assembly below follows
    ! energy_flux.jl expression for expression, including the parenthesization.
    !
    ! CALLING CONVENTION (docs/porting_notes.md D8): contiguous column slices
    ! mass(:,icol) etc. plus the active layer count n.
    !
    ! WORKSPACE (docs/PLAN.md section 4.1, trap 12): Julia's EnergyWorkspace
    ! holds 8 (Ntot,ncol) arrays, of which layer_thickness and
    ! thermal_conductivity are allocated but never used, and
    ! previous_temperature is not a temperature at all -- it is reused as the
    ! scratch copy of the main diagonal, which the Thomas forward sweep
    ! destroys. chion drops the workspace type entirely and uses stack-local
    ! AUTOMATIC arrays of size(mass): Ntot is small (default 15), and stack
    ! locals are automatically private per OpenMP thread.
    !
    ! Five work arrays are genuinely needed: lower, diag, upper, rhs and the
    ! destroyed copy solver_diag. Julia's sixth array (interface_conductance)
    ! is also unnecessary: each interface value feeds exactly two coefficients,
    ! upper(k) of row k and lower(k) of row k+1, so it is consumed immediately
    ! in the assembly loop and never needs to be stored. The arithmetic is
    ! unchanged -- the same expression, evaluated once, in the same order.
    !
    ! TRAP 2 (docs/PLAN.md section 5): the surface fluxes here are linearized
    ! about the OLD interface temperature Ts^n, while snow_surface_fluxes
    ! evaluates the same physics exactly at a known temperature (T0 for bare
    ! ice, Ts^{n+1} for the post-solve vapor mass). The energy and mass
    ! budgets therefore use slightly different latent fluxes. This is
    ! deliberate in Chion.jl and is NOT to be reconciled.
    !
    ! TRAP 6 (docs/PLAN.md section 5): the melting re-solve does not turn
    ! row 1 into a Dirichlet row: T0 is imposed on the INTERFACE, behind the
    ! half-cell conductance, and row 1 keeps its coupling to row 2 -- see
    ! step 5 below.

    use chion_defs, only : wp, wp_acc, CHION_SEB_SEMIX, &
                           chion_const_class, chion_step_forcing_class

    ! SEMIX aerodynamic surface scheme, selected by c%seb_scheme. It replaces
    ! the sensible and turbulent-latent coefficients of step 2 and nothing
    ! else: the conduction assembly, the melting-point re-solve and the whole
    ! firn column below row 1 are untouched (docs/semix_port_scope.md).
    use snow_seb_semix, only : semix_exchange_class, semix_flux_lin_class, &
                               semix_snow_depth, semix_turbulent_exchange, &
                               semix_surface_emissivity, semix_longwave_down, &
                               semix_longwave_linearized

    ! Vapor-pressure / turbulent-latent helpers shared with the unlinearized
    ! twin used for bare ice and post-solve vapor mass (WP5, other half).
    use snow_vapor, only : safe_positive, latent_vapor_flux_linearized, &
                           latent_vapor_flux_lin_class

    implicit none

    private

    public :: snow_energy_result_class
    public :: snow_energy_flux
    public :: solve_tridiagonal_thomas
    public :: snow_thermal_conductivity
    public :: interface_conductance

    type snow_energy_result_class
        ! Mirrors the named tuple built by _energy_flux_result
        ! (energy_flux.jl, 9ec6cc7), in the same order, minus the longwave and
        ! sensible components chion does not use. energy_to_melting is gone
        ! with the surface heat capacity (03bb445).
        !
        ! The two energies are wp_acc: they are differences of large numbers
        ! (Q*dt against the conducted Gs*(T0 - T1)*dt) and they feed the melt
        ! mass directly. See docs/PLAN.md section 3.1.

        logical      :: needs_melt                       ! interface reached T0 this step
        real(wp_acc) :: melt_energy_available            ! [J m-2] surface energy not conducted
        real(wp_acc) :: heating                          ! [J m-2] net energy into the column
        real(wp)     :: surface_flux_constant            ! [W m-2] Q_const
        real(wp)     :: surface_flux_linear              ! [W m-2 K-1] Q_lin
        real(wp)     :: latent_heat_linear_coefficient   ! [W m-2 K-1] echoed through
        real(wp)     :: latent_heat_constant_term        ! [W m-2] echoed through
    end type snow_energy_result_class

contains

    pure function snow_thermal_conductivity(rho,T,rho_i) result(K)
        ! Calonne et al. (2019), Eq. (5); Chion.jl/src/processes/energy_flux.jl:40-59
        ! (dev_nils 49990e6). A logistic blend, centred on 450 kg m-3, of a
        ! snow regression (quadratic in rho) and a firn regression (linear,
        ! reaching the ice value 2.107 at rho_i), each scaled by the
        ! temperature dependence of its conducting phases: ice
        ! k_ice(T) = 9.828 exp(-5.7e-3 T) and air
        ! k_air(T) = 2.334e-3 T^1.5/(164.54 + T), relative to their reference
        ! values 2.107 and 0.024 W m-1 K-1. The coefficients are the
        ! paper's, carried verbatim; the ice density is the model's.
        !
        ! Evaluated in Julia's order and grouping.

        implicit none

        real(wp), intent(IN) :: rho          ! [kg m-3] layer density
        real(wp), intent(IN) :: T            ! [K] layer temperature
        real(wp), intent(IN) :: rho_i        ! [kg m-3] ice density
        real(wp) :: K                        ! [W m-1 K-1]

        ! Local variables
        real(wp), parameter :: K_ICE_REF = 2.107_wp   ! [W m-1 K-1]
        real(wp), parameter :: K_AIR_REF = 0.024_wp   ! [W m-1 K-1]
        real(wp) :: transition, k_ice, k_air, k_snow, k_firn, snow_scale, firn_scale

        transition = 1.0_wp/(1.0_wp + exp(-0.04_wp*(rho - 450.0_wp)))
        k_ice      = 9.828_wp*exp(-5.7e-3_wp*T)
        k_air      = 2.334e-3_wp*T**1.5_wp/safe_positive(164.54_wp + T)
        k_snow     = K_AIR_REF - 1.23e-4_wp*rho + 2.5e-6_wp*(rho*rho)
        k_firn     = K_ICE_REF + 3.618e-3_wp*(rho - rho_i)
        snow_scale = k_ice*k_air/(K_ICE_REF*K_AIR_REF)
        firn_scale = k_ice/K_ICE_REF

        K = (1.0_wp - transition)*snow_scale*k_snow + transition*firn_scale*k_firn

        return

    end function snow_thermal_conductivity

    pure function interface_conductance(K_i,dz_i,K_j,dz_j) result(G)
        ! Chion.jl/src/processes/energy_flux.jl:179-184 (dev_nils 81034fa,
        ! d0146e1): the thermal resistances of the two half layers in series,
        !     G = 1/(dz_i/(2 K_i) + dz_j/(2 K_j)) = 2 K_i K_j/(K_j dz_i + K_i dz_j),
        ! the physical conductance between the two layer centres, so the
        ! interface flux is G*(T_i - T_j) (beta_scale carries no factor 2).
        ! For a uniform column (K, dz) this is K/dz, the same centre-to-centre
        ! flux as the arithmetic form it replaces; across a sharp K contrast
        ! the poorer conductor now controls the flux.
        !
        ! Evaluated as Julia does, ((2 K_i) K_j)/(...).

        implicit none

        real(wp), intent(IN) :: K_i          ! [W m-1 K-1] conductivity, layer i
        real(wp), intent(IN) :: dz_i         ! [m] thickness, layer i
        real(wp), intent(IN) :: K_j          ! [W m-1 K-1] conductivity, layer j
        real(wp), intent(IN) :: dz_j         ! [m] thickness, layer j
        real(wp) :: G                        ! [W m-2 K-1]

        G = 2.0_wp*K_i*K_j/safe_positive(K_j*dz_i + K_i*dz_j)

        return

    end function interface_conductance

    subroutine solve_tridiagonal_thomas(lower,diag,upper,rhs,n)
        ! Thomas algorithm, Chion.jl/src/processes/energy_flux.jl:113-189.
        !
        ! INDEXING CONVENTION, taken verbatim from Julia:
        !     lower(k) is the SUB-diagonal entry of row k+1
        !     upper(k) is the SUPER-diagonal entry of row k
        ! so a row k reads  lower(k-1)*x(k-1) + diag(k)*x(k) + upper(k)*x(k+1).
        ! lower(n) and upper(n) are never referenced.
        !
        ! No pivoting and no zero-check on the diagonal, exactly as in Julia:
        ! the assembled matrix is strictly diagonally dominant by construction
        ! (all off-diagonals are negative, the diagonal is 1 minus their sum,
        ! plus a non-negative surface term), so pivoting is unnecessary.
        !
        ! DESTRUCTIVE: diag is overwritten by the forward sweep and rhs is
        ! overwritten by the solution. Callers must copy diag before each solve.

        implicit none

        real(wp), intent(IN)    :: lower(:)  ! (n) sub-diagonal, lower(k) in row k+1
        real(wp), intent(INOUT) :: diag(:)   ! (n) main diagonal; DESTROYED
        real(wp), intent(IN)    :: upper(:)  ! (n) super-diagonal, upper(k) in row k
        real(wp), intent(INOUT) :: rhs(:)    ! (n) right-hand side; holds the solution
        integer,  intent(IN)    :: n

        ! Local variables
        integer  :: row
        real(wp) :: f

        ! Forward elimination
        do row = 2, n
            f = lower(row-1)/diag(row-1)
            diag(row) = diag(row) - f*upper(row-1)
            rhs(row)  = rhs(row)  - f*rhs(row-1)
        end do

        ! Back substitution
        rhs(n) = rhs(n)/diag(n)

        do row = n-1, 1, -1
            rhs(row) = (rhs(row) - upper(row)*rhs(row+1))/diag(row)
        end do

        return

    end subroutine solve_tridiagonal_thomas

    subroutine snow_energy_flux(mass,density,temperature,t_srf,n,c,forc,albedo, &
                                latent_heat_linear_coefficient, &
                                latent_heat_constant_term,dt_seconds,res)
        ! Advance the temperature profile of one column over one step.
        ! Chion.jl/src/processes/energy_flux.jl _go_energy_flux_resolved!
        ! (main 9ec6cc7 = 03bb445), without the ice substrate rows.
        !
        ! ROBIN SURFACE BOUNDARY (03bb445). t_srf is the physical
        ! atmosphere-snow interface temperature Ts, a state of its own; the
        ! first numerical temperature T1 = temperature(1) sits at the centre
        ! of the top finite-volume cell. The linearized surface energy balance
        ! Q(Ts) = q_const - q_lin*Ts is closed against the half-cell
        ! conductance Gs = 2 K1/dz1:
        !     Q(Ts) = Gs*(Ts - T1)  =>  Ts = a + b*T1,
        !     den = q_lin + Gs,  a = q_const/den,  b = Gs/den.
        ! Eliminating Ts leaves no surface heat capacity: row 1 is a regular
        ! cell whose top face carries Gs*(Ts - T1) = Q(Ts). One layer goes
        ! through the same matrix (no closed-form special case).
        !
        ! Melting: if Ts = a + b*T1^{n+1} > T0, the column is re-solved with
        ! Ts = T0 held behind Gs, Tsrf = T0 and
        !     melt energy = max((Q(T0) - Gs*(T0 - T1^{n+1}))*dt, 0),
        ! the surface energy not conducted into the snow.
        !
        ! Deviation from Julia, structural only: Julia flattens the forcing
        ! into fifteen positional scalars plus use_* flags because the kernel
        ! must be GPU-callable. chion passes chion_step_forcing_class, which is
        ! already the per-column, per-substep backend contract and carries the
        ! has_* flags with it (docs/PLAN.md section 1). Julia's mass_w argument
        ! is dropped: the routine never reads it.

        implicit none

        real(wp), intent(IN)    :: mass(:)        ! (Ntot) [kg m-2] solid mass
        real(wp), intent(IN)    :: density(:)     ! (Ntot) [kg m-3] layer density
        real(wp), intent(INOUT) :: temperature(:) ! (Ntot) [K] layer (cell-centre) temperature
        real(wp), intent(INOUT) :: t_srf          ! [K] interface temperature Ts
        integer,  intent(IN)    :: n              ! number of active layers

        type(chion_const_class),        intent(IN) :: c
        type(chion_step_forcing_class), intent(IN) :: forc

        real(wp), intent(IN) :: albedo            ! [1] diagnosed surface albedo
        real(wp), intent(IN) :: latent_heat_linear_coefficient  ! [W m-2 K-1] precip
        real(wp), intent(IN) :: latent_heat_constant_term       ! [W m-2] precip
        real(wp), intent(IN) :: dt_seconds        ! [s]

        type(snow_energy_result_class), intent(OUT) :: res

        ! Work arrays. Automatic, stack-local, OpenMP-private by construction.
        ! Five arrays; see the workspace note in the module header.
        real(wp), dimension(size(mass)) :: lower, diag, upper, rhs, solver_diag

        ! Local variables
        integer  :: k
        real(wp) :: m1, Ts_n, T1_n
        real(wp) :: Ts_sq, Ts_cube, Ts_fourth
        real(wp) :: sw_abs, lw_const, lw_lin, sh_const, sh_lin
        real(wp) :: lh_turb_const, lh_turb_lin, lh_const, lh_lin
        real(wp) :: q_const, q_lin
        real(wp) :: dz_prev, dz_k, K_prev, K_k, G_k
        real(wp) :: beta_scale, beta_1, beta_km1, beta_k
        real(wp) :: G_s, surface_den, surface_const, surface_coef, boundary_term
        real(wp) :: ts_new

        logical  :: uses_semix_seb

        type(latent_vapor_flux_lin_class) :: lh_coef
        type(semix_exchange_class)        :: sx
        type(semix_flux_lin_class)        :: lw_coef

        ! === Step 0: early exit ==============================================
        ! energy_flux.jl: n_snow = n if n > 0 and mass(1) > 0, else 0. NOTE the
        ! threshold is mass(1) <= 0, NOT TOL_EMPTY_LAYER as in surface_has_snow.
        ! The two guards differ deliberately -- see docs/PLAN.md section 5,
        ! item 1. Nothing is written: temperature and t_srf are left untouched,
        ! and the two precipitation latent coefficients are echoed through.

        res%needs_melt                     = .FALSE.
        res%melt_energy_available          = 0.0_wp_acc
        res%heating                        = 0.0_wp_acc
        res%surface_flux_constant          = 0.0_wp
        res%surface_flux_linear            = 0.0_wp
        res%latent_heat_linear_coefficient = latent_heat_linear_coefficient
        res%latent_heat_constant_term      = latent_heat_constant_term

        if (n .le. 0) return
        if (mass(1) .le. 0.0_wp) return

        ! === Step 1: surface scalars =========================================
        ! The fluxes are linearized about the previous INTERFACE temperature
        ! Ts^n; the top cell's conductivity is taken at its own T1^n.

        m1   = safe_positive(mass(1))
        Ts_n = t_srf
        T1_n = temperature(1)

        Ts_sq     = Ts_n*Ts_n
        Ts_cube   = Ts_sq*Ts_n
        Ts_fourth = Ts_sq*Ts_sq

        ! SEMIX exchange coefficients, built ONCE at the linearization point
        ! Ts^n, which is the temperature the whole of step 2 linearizes about.
        uses_semix_seb = (c%seb_scheme .eq. CHION_SEB_SEMIX)

        if (uses_semix_seb) then
            sx = semix_turbulent_exchange(c,semix_snow_depth(mass,density,n), &
                                          forc%air_temperature,Ts_n, &
                                          forc%wind_speed,forc%air_pressure, &
                                          forc%relative_humidity, &
                                          forc%has_relative_humidity)
        end if

        ! === Step 2: linearized surface energy balance =======================
        ! Each term is either taken from the forcing (has_* true) or
        ! parameterized internally.

        ! Shortwave. NOTE there is deliberately NO max(shortwave_down,0) here,
        ! while surface_fluxes.jl:74 (the bare-ice twin) does apply one. The
        ! asymmetry is real and is preserved: docs/PLAN.md section 4, WP5.
        if (forc%has_q_sw_net) then
            sw_abs = forc%q_sw_net
        else
            sw_abs = (1.0_wp - min(max(albedo,0.0_wp),1.0_wp))*forc%shortwave_down
        end if

        ! Longwave. The emitted term eps_snow*sigma*Ts^4 is linearized about
        ! Ts^n as Ts^4 ~ 4*Ts_n^3*Ts - 3*Ts_n^4, so +3*sigma*eps_snow*Ts_n^4 goes
        ! to the constant part and 4*sigma*eps_snow*Ts_n^3 to the linear part.
        !
        ! Under semix the same linearization applies but the DOWNWELLING flux
        ! is absorbed with the surface emissivity too (ebal num_lw/denom_lw),
        ! where BESSI absorbs it in full. Always the snow emissivity here:
        ! snow_energy_flux has already returned on a bare column.
        !
        ! The two branches are written out separately rather than factored
        ! through a shared lw_down: the BESSI else-branch groups its two terms
        ! inside one multiply by sigma, and refactoring that grouping would
        ! change the bessi answer in single precision.
        if (uses_semix_seb) then
            lw_coef = semix_longwave_linearized(c, &
                          semix_surface_emissivity(c,.TRUE.), &
                          semix_longwave_down(c,forc%q_lw_down,forc%has_q_lw_down, &
                                              forc%air_temperature), &
                          Ts_n)
            lw_const = lw_coef%constant
            lw_lin   = lw_coef%linear
        else if (forc%has_q_lw_down) then
            lw_const = forc%q_lw_down + c%sigma_sb*c%eps_snow*3.0_wp*Ts_fourth
        else
            lw_const = c%sigma_sb*(c%eps_air*forc%air_temperature**4 &
                                   + c%eps_snow*3.0_wp*Ts_fourth)
        end if

        if (.not. uses_semix_seb) lw_lin = c%sigma_sb*c%eps_snow*4.0_wp*Ts_cube

        ! Sensible heat. A prescribed flux still wins over either scheme.
        if (forc%has_q_sh) then
            sh_const = forc%q_sh
            sh_lin   = 0.0_wp
        else if (uses_semix_seb) then
            ! ebal num_sh/denom_sh (smb_ebal.f90:125-126). Same shape as the
            ! BESSI branch below, with the aerodynamic f_sh in place of D_sh.
            sh_const = forc%air_temperature*sx%f_sh
            sh_lin   = sx%f_sh
        else
            sh_const = forc%air_temperature*c%D_sh
            sh_lin   = c%D_sh
        end if

        ! Turbulent latent heat. Three-way, in this order (energy_flux.jl).
        if (forc%has_q_lh) then
            lh_turb_const = forc%q_lh
            lh_turb_lin   = 0.0_wp
        else if (uses_semix_seb) then
            ! ebal num_lh/denom_lh (smb_ebal.f90:122-123), which ARE chion's
            ! q_const/q_lin contributions -- coupling decision alpha. Note the
            ! sign flip: SEMIX counts the latent flux positive away from the
            ! surface, chion positive into it.
            !
            ! No third branch is needed for missing humidity forcing:
            ! semix_turbulent_exchange already zeroes f_lh in that case, which
            ! collapses both terms to zero exactly as the BESSI selection does.
            lh_turb_const = -sx%f_lh*(sx%qsat - sx%dqsatdT*Ts_n - sx%q_air)
            lh_turb_lin   =  sx%f_lh*sx%dqsatdT
        else if (forc%has_relative_humidity) then
            lh_coef       = latent_vapor_flux_linearized(Ts_n,c,forc%air_temperature, &
                                                         forc%relative_humidity, &
                                                         forc%air_pressure)
            lh_turb_const = lh_coef%constant
            lh_turb_lin   = lh_coef%linear
        else
            lh_turb_const = 0.0_wp
            lh_turb_lin   = 0.0_wp
        end if

        ! The precipitation heat coefficients passed in by the caller are ADDED
        ! to the turbulent term, not replaced by it: a prescribed q_lh is an
        ! ADDITIONAL flux (energy_flux.jl, docs energy.md).
        lh_const = latent_heat_constant_term + lh_turb_const
        lh_lin   = latent_heat_linear_coefficient + lh_turb_lin

        q_const = sh_const + lw_const + sw_abs + lh_const
        q_lin   = sh_lin + lw_lin + lh_lin

        res%surface_flux_constant = q_const
        res%surface_flux_linear   = q_lin

        ! === Step 3: matrix assembly =========================================
        !
        ! Julia stores the interface conductances in their own array and
        ! assembles in a second pass. Here each interface value is consumed
        ! immediately by the two coefficients it feeds -- upper(k-1) of row k-1
        ! and lower(k-1) of row k -- which removes one work array without
        ! changing a single floating-point operation.

        lower = 0.0_wp
        diag  = 0.0_wp
        upper = 0.0_wp
        rhs   = 0.0_wp

        ! -dt G/(ci m) off the diagonal: G is the full centre-to-centre
        ! conductance (energy_flux.jl, physical interface conductance).
        beta_scale = -dt_seconds/c%ci

        ! NOTE the surface layer thickness uses the SAFE-POSITIVE mass m1,
        ! while every other layer uses its raw mass. Each layer's conductivity
        ! is at its own start-of-step temperature.
        dz_prev = m1/safe_positive(density(1))
        K_prev  = snow_thermal_conductivity(density(1),T1_n,c%rho_i)

        ! Robin boundary terms from the top cell's half-thickness.
        G_s           = 2.0_wp*K_prev/safe_positive(dz_prev)
        surface_den   = safe_positive(q_lin + G_s)
        surface_const = q_const/surface_den
        surface_coef  = G_s/surface_den
        beta_1        = beta_scale/m1
        boundary_term = beta_1*G_s

        rhs(1) = T1_n - boundary_term*surface_const

        do k = 2, n

            dz_k = mass(k)/safe_positive(density(k))
            K_k  = snow_thermal_conductivity(density(k),temperature(k),c%rho_i)

            G_k = interface_conductance(K_prev,dz_prev,K_k,dz_k)

            beta_km1 = beta_scale/safe_positive(mass(k-1))
            beta_k   = beta_scale/safe_positive(mass(k))

            upper(k-1) = beta_km1*G_k      ! super-diagonal of row k-1
            lower(k-1) = beta_k*G_k        ! sub-diagonal   of row k

            rhs(k) = temperature(k)

            dz_prev = dz_k
            K_prev  = K_k

        end do

        ! Diagonal. Row 1 is 1 - bt*(1 - b), written as 1 - bt*(q_lin/den) so
        ! that a vanishing top cell (bt and b both huge/near 1) keeps an exact
        ! 1 when q_lin = 0 instead of the difference of two enormous terms.
        ! Row n is zero-flux at the bottom (no lower(n), no upper(n)).
        diag(1) = 1.0_wp - boundary_term*(q_lin/surface_den)

        if (n .gt. 1) then
            diag(1) = diag(1) - upper(1)
            diag(n) = 1.0_wp - lower(n-1)
        end if

        do k = 2, n-1
            diag(k) = 1.0_wp - lower(k-1) - upper(k)
        end do

        ! === Step 4: first solve =============================================

        solver_diag(1:n) = diag(1:n)
        call solve_tridiagonal_thomas(lower,solver_diag,upper,rhs,n)

        ts_new = surface_const + surface_coef*rhs(1)

        ! === Step 5: melting re-solve ========================================
        ! The interface is held at T0 as a Dirichlet value BEHIND the
        ! half-cell conductance: rhs(1) gains -bt*T0 and the diagonal loses
        ! the eliminated Robin coefficient bt*b, while row 1 KEEPS its
        ! conduction coupling to row 2.

        if (ts_new .gt. c%T0) then

            res%needs_melt = .TRUE.

            ! Rebuild the rhs from the ORIGINAL temperatures. temperature() has
            ! not been written yet, which is exactly why it is only updated at
            ! step 6.
            do k = 1, n
                rhs(k) = temperature(k)
            end do
            rhs(1) = rhs(1) - boundary_term*c%T0

            solver_diag(1:n) = diag(1:n)
            solver_diag(1)   = diag(1) - boundary_term*surface_coef

            call solve_tridiagonal_thomas(lower,solver_diag,upper,rhs,n)

            do k = 1, n
                if (rhs(k) .gt. c%T0) rhs(k) = c%T0
            end do

            ts_new = c%T0

            res%heating = real(dt_seconds,wp_acc) &
                          *(real(q_const,wp_acc) - real(q_lin,wp_acc)*real(c%T0,wp_acc))

            ! Surface energy at T0 not conducted into the top cell.
            res%melt_energy_available = &
                max((real(q_const,wp_acc) - real(q_lin,wp_acc)*real(c%T0,wp_acc) &
                     - real(G_s,wp_acc)*(real(c%T0,wp_acc) - real(rhs(1),wp_acc))) &
                    *real(dt_seconds,wp_acc), 0.0_wp_acc)

        else

            do k = 1, n
                if (rhs(k) .gt. c%T0) rhs(k) = c%T0
            end do

            res%heating = real(dt_seconds,wp_acc) &
                          *(real(q_const,wp_acc) - real(q_lin,wp_acc)*real(ts_new,wp_acc))

        end if

        ! === Step 6: write state =============================================

        do k = 1, n
            temperature(k) = rhs(k)
        end do

        t_srf = ts_new

        return

    end subroutine snow_energy_flux

end module snow_energy
