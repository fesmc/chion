module snow_energy
    ! Implicit (backward-Euler) snow temperature solve: linearized surface
    ! energy balance on the top row, vertical heat conduction below, zero-flux
    ! bottom, solved with the Thomas algorithm, plus a two-pass melting-point
    ! re-solve.
    !
    ! Port of Chion.jl/src/processes/energy_flux.jl:
    !   _go_energy_flux_resolved!            -> snow_energy_flux
    !   _thomas_forward! / _thomas_backward! -> solve_tridiagonal_thomas
    !   _snow_thermal_conductivity           -> snow_thermal_conductivity
    !   interface_conductance                -> interface_conductance
    !   _clamp_to_melt!                      -> inlined (clamp loop)
    !   _residual_melt_energy                -> inlined (step 6)
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
    ! about the OLD surface temperature T^n, while snow_surface_fluxes
    ! evaluates the same physics exactly at a known temperature (T0 for bare
    ! ice, T^{n+1} for the post-solve vapor mass). The energy and mass budgets
    ! therefore use slightly different latent fluxes. This is deliberate in
    ! Chion.jl and is NOT to be reconciled.
    !
    ! TRAP 6 (docs/PLAN.md section 5): the melting-point re-solve is NOT a
    ! Dirichlet row -- see the comment at step 5 below.

    use chion_defs, only : wp, wp_acc, chion_const_class, chion_step_forcing_class

    ! Vapor-pressure / turbulent-latent helpers shared with the unlinearized
    ! twin used for bare ice and post-solve vapor mass (WP5, other half).
    use snow_surface_fluxes, only : safe_positive, latent_vapor_flux_linearized, &
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
        ! (energy_flux.jl:197-239), field for field and in the same order.
        !
        ! The three energies are wp_acc: they are differences of large numbers
        ! (T0 - T)*ci*m and Q*dt, and they feed the melt mass directly.
        ! See docs/PLAN.md section 3.1.

        logical      :: needs_melt                       ! surface reached T0 this step
        real(wp_acc) :: energy_to_melting                ! [J m-2] spent reaching T0
        real(wp_acc) :: melt_energy_available            ! [J m-2] left over for melting
        real(wp_acc) :: heating                          ! [J m-2] net energy into the column
        real(wp)     :: surface_flux_constant            ! [W m-2] Q_const
        real(wp)     :: surface_flux_linear              ! [W m-2 K-1] Q_lin
        real(wp)     :: latent_heat_linear_coefficient   ! [W m-2 K-1] echoed through
        real(wp)     :: latent_heat_constant_term        ! [W m-2] echoed through
    end type snow_energy_result_class

contains

    pure function snow_thermal_conductivity(rho,Ki) result(K)
        ! Chion.jl/src/processes/energy_flux.jl:40-41:  K = Ki*(rho*1e-3)^1.88
        ! Note the exponent is applied to density in g cm-3, not kg m-3.

        implicit none

        real(wp), intent(IN) :: rho          ! [kg m-3] layer density
        real(wp), intent(IN) :: Ki           ! [W m-1 K-1] conductivity of ice
        real(wp) :: K

        K = Ki*(rho*1.0e-3_wp)**1.88_wp

        return

    end function snow_thermal_conductivity

    pure function interface_conductance(K_i,dz_i,K_j,dz_j) result(G)
        ! Chion.jl/src/processes/energy_flux.jl:104-105:
        !     G = (K_i*dz_i + K_j*dz_j) / safe_positive((dz_i + dz_j)^2)
        !
        ! NOTE the denominator is the SQUARE of the summed thickness. The
        ! 1/(dz_i+dz_j) that turns a conductivity into a conductance is already
        ! contained in it -- do not add another division by dz. For a uniform
        ! column (K, dz) this gives G = K/(2 dz), and the assembled interface
        ! flux 2*G*(T_i - T_j) = (K/dz)*(T_i - T_j) is the correct
        ! centre-to-centre conductive flux.

        implicit none

        real(wp), intent(IN) :: K_i          ! [W m-1 K-1] conductivity, layer i
        real(wp), intent(IN) :: dz_i         ! [m] thickness, layer i
        real(wp), intent(IN) :: K_j          ! [W m-1 K-1] conductivity, layer j
        real(wp), intent(IN) :: dz_j         ! [m] thickness, layer j
        real(wp) :: G                        ! [W m-2 K-1]

        G = (K_i*dz_i + K_j*dz_j)/safe_positive((dz_i + dz_j)**2)

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
        ! Chion.jl/src/processes/energy_flux.jl:315-526.
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
        real(wp), intent(INOUT) :: temperature(:) ! (Ntot) [K] layer temperature
        real(wp), intent(INOUT) :: t_srf          ! [K] surface temperature
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
        real(wp) :: m1, Tn, lambda
        real(wp) :: Tn_sq, Tn_cube, Tn_fourth
        real(wp) :: sw_abs, lw_const, lw_lin, sh_const, sh_lin
        real(wp) :: lh_turb_const, lh_turb_lin, lh_const, lh_lin
        real(wp) :: q_const, q_lin, rhs_surf, diag_surf
        real(wp) :: dz_prev, dz_k, K_prev, K_k, G_k
        real(wp) :: beta_scale, beta_km1, beta_k
        real(wp) :: t_new

        type(latent_vapor_flux_lin_class) :: lh_coef

        ! === Step 0: early exit ==============================================
        ! energy_flux.jl:343-355. NOTE the threshold is mass(1) <= 0, NOT
        ! TOL_EMPTY_LAYER as in surface_has_snow. The two guards differ
        ! deliberately -- see docs/PLAN.md section 5, item 1.
        ! Nothing is written: temperature and t_srf are left untouched, and the
        ! two precipitation latent coefficients are echoed through unchanged.

        res%needs_melt                     = .FALSE.
        res%energy_to_melting              = 0.0_wp_acc
        res%melt_energy_available          = 0.0_wp_acc
        res%heating                        = 0.0_wp_acc
        res%surface_flux_constant          = 0.0_wp
        res%surface_flux_linear            = 0.0_wp
        res%latent_heat_linear_coefficient = latent_heat_linear_coefficient
        res%latent_heat_constant_term      = latent_heat_constant_term

        if (n .le. 0) return
        if (mass(1) .le. 0.0_wp) return

        ! === Step 1: surface layer scalars ===================================

        m1     = safe_positive(mass(1))
        Tn     = temperature(1)
        lambda = dt_seconds/c%ci/m1

        Tn_sq     = Tn*Tn
        Tn_cube   = Tn_sq*Tn
        Tn_fourth = Tn_sq*Tn_sq

        ! === Step 2: linearized surface energy balance =======================
        ! energy_flux.jl:370-392. Each term is either taken from the forcing
        ! (has_* true) or parameterized internally.

        ! Shortwave. NOTE there is deliberately NO max(shortwave_down,0) here,
        ! while surface_fluxes.jl:74 (the bare-ice twin) does apply one. The
        ! asymmetry is real and is preserved: docs/PLAN.md section 4, WP5.
        if (forc%has_q_sw_net) then
            sw_abs = forc%q_sw_net
        else
            sw_abs = (1.0_wp - min(max(albedo,0.0_wp),1.0_wp))*forc%shortwave_down
        end if

        ! Longwave. The emitted term eps_snow*sigma*T^4 is linearized about Tn
        ! as T^4 ~ 4*Tn^3*T - 3*Tn^4, so +3*sigma*eps_snow*Tn^4 goes to the
        ! constant part and 4*sigma*eps_snow*Tn^3 to the linear part.
        if (forc%has_q_lw_down) then
            lw_const = forc%q_lw_down + c%sigma_sb*c%eps_snow*3.0_wp*Tn_fourth
        else
            lw_const = c%sigma_sb*(c%eps_air*forc%air_temperature**4 &
                                   + c%eps_snow*3.0_wp*Tn_fourth)
        end if
        lw_lin = c%sigma_sb*c%eps_snow*4.0_wp*Tn_cube

        ! Sensible heat
        if (forc%has_q_sh) then
            sh_const = forc%q_sh
            sh_lin   = 0.0_wp
        else
            sh_const = forc%air_temperature*c%D_sh
            sh_lin   = c%D_sh
        end if

        ! Turbulent latent heat. Three-way, in this order (energy_flux.jl:379).
        if (forc%has_q_lh) then
            lh_turb_const = forc%q_lh
            lh_turb_lin   = 0.0_wp
        else if (forc%has_relative_humidity) then
            lh_coef       = latent_vapor_flux_linearized(Tn,c,forc%air_temperature, &
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
        ! ADDITIONAL flux (energy_flux.jl:386-387, docs energy.md).
        lh_const = latent_heat_constant_term + lh_turb_const
        lh_lin   = latent_heat_linear_coefficient + lh_turb_lin

        q_const = sh_const + lw_const + sw_abs + lh_const
        q_lin   = sh_lin + lw_lin + lh_lin

        rhs_surf  = lambda*q_const
        diag_surf = lambda*q_lin

        res%surface_flux_constant = q_const
        res%surface_flux_linear   = q_lin

        ! === Step 3a: single-layer closed form ===============================
        ! energy_flux.jl:398-429. A separate branch, NOT the solver with n = 1.

        if (n .eq. 1) then

            t_new = (Tn + rhs_surf)/safe_positive(1.0_wp + diag_surf)

            if (t_new .gt. c%T0) then
                res%needs_melt        = .TRUE.
                res%energy_to_melting = (real(c%T0,wp_acc) - real(Tn,wp_acc)) &
                                        *real(c%ci,wp_acc)*real(m1,wp_acc)
                t_new                 = c%T0
                res%heating           = res%energy_to_melting
            else
                res%heating = real(dt_seconds,wp_acc) &
                              *(real(q_const,wp_acc) - real(q_lin,wp_acc)*real(t_new,wp_acc))
            end if

            t_new          = min(t_new,c%T0)
            temperature(1) = t_new
            t_srf          = t_new

            call residual_melt_energy(res,q_const,q_lin,t_srf,dt_seconds)

            return

        end if

        ! === Step 3b: multilayer matrix assembly =============================
        ! energy_flux.jl:432-477.
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

        beta_scale = -2.0_wp*dt_seconds/c%ci

        ! NOTE the surface layer thickness uses the SAFE-POSITIVE mass m1,
        ! while every other layer uses its raw mass (energy_flux.jl:433 vs 439).
        dz_prev = m1/safe_positive(density(1))
        K_prev  = snow_thermal_conductivity(density(1),c%Ki)

        rhs(1) = Tn + rhs_surf

        do k = 2, n

            dz_k = mass(k)/safe_positive(density(k))
            K_k  = snow_thermal_conductivity(density(k),c%Ki)

            G_k = interface_conductance(K_prev,dz_prev,K_k,dz_k)

            beta_km1 = beta_scale/safe_positive(mass(k-1))
            beta_k   = beta_scale/safe_positive(mass(k))

            upper(k-1) = beta_km1*G_k      ! super-diagonal of row k-1
            lower(k-1) = beta_k*G_k        ! sub-diagonal   of row k

            rhs(k) = temperature(k)

            dz_prev = dz_k
            K_prev  = K_k

        end do

        ! Diagonal. Row 1 carries the linearized surface flux; row n is
        ! zero-flux at the bottom (no lower(n), no upper(n)).
        diag(1) = 1.0_wp - upper(1) + diag_surf
        diag(n) = 1.0_wp - lower(n-1)

        do k = 2, n-1
            diag(k) = 1.0_wp - lower(k-1) - upper(k)
        end do

        ! === Step 4: first solve =============================================

        solver_diag(1:n) = diag(1:n)
        call solve_tridiagonal_thomas(lower,solver_diag,upper,rhs,n)

        ! === Step 5: melting-point re-solve ==================================
        ! energy_flux.jl:482-501.
        !
        ! TRAP 6 (docs/PLAN.md section 5): this is NOT a strict Dirichlet row.
        ! Row 1 keeps its conduction coupling to row 2 -- only the surface FLUX
        ! feedback diag_surf is removed from the diagonal, and the rhs is
        ! replaced by T0. Turning it into a true Dirichlet row (diag = 1,
        ! upper = 0) would change the physics and is explicitly not allowed
        ! without asking (docs/PLAN.md section 4.1).
        !
        ! energy_to_melting is accumulated in two contributions: what pass 1
        ! spent bringing the surface from Tn up to T0, plus what pass 2's
        ! conduction removes from the pinned surface.

        if (rhs(1) .gt. c%T0) then

            res%needs_melt        = .TRUE.
            res%energy_to_melting = (real(c%T0,wp_acc) - real(Tn,wp_acc)) &
                                    *real(c%ci,wp_acc)*real(m1,wp_acc)

            ! Rebuild the rhs from the ORIGINAL temperatures. temperature() has
            ! not been written yet, which is exactly why it is only updated at
            ! step 6.
            do k = 1, n
                rhs(k) = temperature(k)
            end do
            rhs(1) = c%T0

            solver_diag(1:n) = diag(1:n)
            solver_diag(1)   = solver_diag(1) - diag_surf

            call solve_tridiagonal_thomas(lower,solver_diag,upper,rhs,n)

            res%energy_to_melting = res%energy_to_melting &
                                    + (real(c%T0,wp_acc) - real(rhs(1),wp_acc)) &
                                      *real(c%ci,wp_acc)*real(m1,wp_acc)

            rhs(1) = c%T0

            do k = 1, n
                if (rhs(k) .gt. c%T0) rhs(k) = c%T0
            end do

            res%heating = res%energy_to_melting

        else

            do k = 1, n
                if (rhs(k) .gt. c%T0) rhs(k) = c%T0
            end do

            res%heating = real(dt_seconds,wp_acc) &
                          *(real(q_const,wp_acc) - real(q_lin,wp_acc)*real(rhs(1),wp_acc))

        end if

        ! === Step 6: write state and diagnose residual melt energy ===========

        do k = 1, n
            temperature(k) = rhs(k)
        end do

        t_srf = rhs(1)

        call residual_melt_energy(res,q_const,q_lin,t_srf,dt_seconds)

        return

    end subroutine snow_energy_flux

    subroutine residual_melt_energy(res,q_const,q_lin,t_surface,dt_seconds)
        ! Chion.jl/src/processes/energy_flux.jl:247-274.
        ! Energy left over once the surface has been brought to the melting
        ! point. Zero unless the step actually reached melting.

        implicit none

        type(snow_energy_result_class), intent(INOUT) :: res
        real(wp),                       intent(IN)    :: q_const     ! [W m-2]
        real(wp),                       intent(IN)    :: q_lin       ! [W m-2 K-1]
        real(wp),                       intent(IN)    :: t_surface   ! [K]
        real(wp),                       intent(IN)    :: dt_seconds  ! [s]

        if (.not. res%needs_melt) then
            res%melt_energy_available = 0.0_wp_acc
            return
        end if

        res%melt_energy_available = &
            max((real(q_const,wp_acc) - real(q_lin,wp_acc)*real(t_surface,wp_acc)) &
                *real(dt_seconds,wp_acc) - res%energy_to_melting, 0.0_wp_acc)

        return

    end subroutine residual_melt_energy

end module snow_energy
