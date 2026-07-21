module snow_refreezing
    ! Refreezing of retained liquid water inside a single snowpack column.
    !
    ! Port of Chion.jl/src/processes/refreezing.jl (_go_refreezing!, the kernel
    ! that step.jl calls directly).
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n. See docs/porting_notes.md D8.
    !
    ! SCHEME: every active layer is treated INDEPENDENTLY -- there is no
    ! vertical coupling of any kind, and no ordering dependence. A layer
    ! refreezes water only if it simultaneously has solid mass, liquid water,
    ! and is below the melting point. The available cold content
    !     Q_cold = (T0 - T)*ci*m_s
    ! is compared against the latent heat that would be released by freezing
    ! all the water,
    !     Q_lat  = m_w*Lm.
    !
    !   Q_cold <  Q_lat  (partial):  dm = Q_cold/Lm freezes, T is set to T0
    !                                exactly, and the remaining water stays.
    !   Q_cold >= Q_lat  (complete): all the water freezes and the layer warms
    !                                to the mixing temperature.
    !
    ! Q_cold MUST be evaluated in wp_acc. At 273 K the sp spacing is ~3e-5 K,
    ! so (T0 - T) is pure noise for layers close to melting. See
    ! docs/PLAN.md section 3.1 and docs/porting_notes.md D1.
    !
    ! TRAP (deliberate, do NOT "fix"): the density update is capped at rho_i
    ! while the solid mass still gains the FULL refrozen increment. Mass is
    ! therefore conserved but volume is not -- once a layer reaches ice
    ! density, further refreezing adds mass without adding thickness. This is
    ! Chion.jl behaviour and is listed as trap 7 in docs/PLAN.md section 5;
    ! changing it is explicitly a "not allowed without asking" item under the
    ! section 4.1 cleanup policy.

    use chion_defs, only : wp, wp_acc

    implicit none

    private

    public :: apply_refreezing

contains

    subroutine apply_refreezing(mass,mass_w,density,temperature,n,T0,ci,Lm,rho_i,refrozen)
        ! Chion.jl/src/processes/refreezing.jl:12-69.
        !
        ! refrozen is the mass frozen during this call [kg m-2] (NOT a running
        ! total). The released latent heat is refrozen*Lm; the caller computes
        ! that if it needs it, exactly as go_refreezing! does in Julia.

        implicit none

        real(wp),     intent(INOUT) :: mass(:)         ! (Ntot) [kg m-2] solid
        real(wp),     intent(INOUT) :: mass_w(:)       ! (Ntot) [kg m-2] liquid
        real(wp),     intent(INOUT) :: density(:)      ! (Ntot) [kg m-3]
        real(wp),     intent(INOUT) :: temperature(:)  ! (Ntot) [K]
        integer,      intent(IN)    :: n               ! active layer count
        real(wp),     intent(IN)    :: T0              ! [K] melting temperature
        real(wp),     intent(IN)    :: ci              ! [J kg-1 K-1]
        real(wp),     intent(IN)    :: Lm              ! [J kg-1]
        real(wp),     intent(IN)    :: rho_i           ! [kg m-3]
        real(wp_acc), intent(OUT)   :: refrozen        ! [kg m-2] this call only

        ! Local variables
        integer      :: k
        real(wp)     :: m_s, m_w, temp, rho
        real(wp_acc) :: q_cold, q_lat, dm, rho_new, temp_new

        refrozen = 0.0_wp_acc

        do k = 1, n

            ! Pre-update copies. The density update below uses the values from
            ! BEFORE mass and mass_w are overwritten, so these copies are
            ! load-bearing, not a convenience.
            m_s  = mass(k)
            m_w  = mass_w(k)
            temp = temperature(k)
            rho  = density(k)

            ! Strict triple entry condition: exact zeros, NOT TOL_TINY. This
            ! is one of the places where Chion.jl deliberately uses "> 0"
            ! rather than one of the two empty-layer tolerances.
            ! See docs/PLAN.md section 5, item 1.
            if (.not. (m_s .gt. 0.0_wp .and. m_w .gt. 0.0_wp .and. temp .lt. T0)) cycle

            q_cold = (real(T0,wp_acc) - real(temp,wp_acc)) &
                     *real(ci,wp_acc)*real(m_s,wp_acc)
            q_lat  = real(m_w,wp_acc)*real(Lm,wp_acc)

            if (q_cold .lt. q_lat) then

                ! --- Partial refreezing: cold content is the limit ---------
                dm = q_cold/real(Lm,wp_acc)

                temperature(k) = T0

                ! Density cap at rho_i; see the module header trap note.
                rho_new = real(rho,wp_acc)*(dm + real(m_s,wp_acc))/real(m_s,wp_acc)
                density(k) = real(min(rho_new,real(rho_i,wp_acc)),wp)

                mass(k)   = real(real(m_s,wp_acc) + dm,wp)
                mass_w(k) = real(real(m_w,wp_acc) - dm,wp)

                refrozen = refrozen + dm

            else

                ! --- Complete refreezing: all liquid water freezes ---------
                ! Mixing temperature from the PRE-update m_s, m_w and T.
                temp_new = (real(m_w,wp_acc)*real(Lm,wp_acc)/real(ci,wp_acc) &
                            + real(m_w,wp_acc)*real(T0,wp_acc)               &
                            + real(temp,wp_acc)*real(m_s,wp_acc))            &
                           /(real(m_w,wp_acc) + real(m_s,wp_acc))
                temperature(k) = real(temp_new,wp)

                rho_new = real(rho,wp_acc)*(real(m_w,wp_acc) + real(m_s,wp_acc)) &
                          /real(m_s,wp_acc)
                density(k) = real(min(rho_new,real(rho_i,wp_acc)),wp)

                mass(k)   = real(real(m_s,wp_acc) + real(m_w,wp_acc),wp)
                mass_w(k) = 0.0_wp

                refrozen = refrozen + real(m_w,wp_acc)

            end if

        end do

        return

    end subroutine apply_refreezing

end module snow_refreezing
