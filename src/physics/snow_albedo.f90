module snow_albedo
    ! Surface albedo: constant, dynamic and prescribed schemes.
    !
    ! Port of Chion.jl/src/processes/albedo.jl.
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n, per docs/porting_notes.md D8. The per-column albedo scalar is passed
    ! as intent(INOUT), mirroring Julia's in-place mutation of state.albedo.
    !
    ! PRESERVED QUIRKS (docs/PLAN.md section 5):
    !   * item 5 -- the dynamic aging law has NO dt. It decays once per CALL,
    !     so albedo_update must be called exactly once per model step. Do not
    !     "fix" this by scaling with dt; that is explicitly not allowed
    !     without asking (docs/PLAN.md section 4.1).
    !   * item 9 -- CHION_ALBEDO_PRESCRIBED takes the DYNAMIC code path inside
    !     this module. The caller (WP8) overrides the result afterwards, and
    !     only when forc%has_prescribed_albedo is true; when it is false the
    !     scheme silently behaves as dynamic. Every branch here therefore
    !     tests "is this the CONSTANT scheme?" and never "is this dynamic?".
    !   * item 1 -- two different empty-layer thresholds are used, one line
    !     apart: TOL_EMPTY_LAYER gates "no surface snow" (-> alpha_ice), while
    !     TOL_TINY gates the surface liquid-water content. Not interchangeable.
    !   * item 4 -- c%max_lwc_albedo is NOT the percolation max_lwc. They share
    !     a default of 0.1 and nothing else.
    !
    ! Note on magnitudes: a single snowfall event can brighten the surface by
    ! at most (alpha_dry - alpha_wet), because the brightening increment is
    ! (alpha_dry-alpha_wet)*(1-exp(-dm/3)) < (alpha_dry-alpha_wet).

    use chion_defs, only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER, &
                           chion_const_class, CHION_ALBEDO_CONSTANT
    use snow_column_utils, only : layer_lwc

    implicit none

    private

    ! Dynamic aging law coefficients (Chion.jl albedo.jl:114). Bare magic
    ! numbers upstream; named here per docs/PLAN.md section 4.1.
    real(wp), parameter, public :: ALBEDO_AGING_TEMP_COEFF = 1.35e-3_wp  ! [K-1] per call
    real(wp), parameter, public :: ALBEDO_AGING_OFFSET     = 0.0278_wp   ! [1]  per call

    ! Snowfall brightening e-folding mass (Chion.jl albedo.jl:77).
    real(wp), parameter, public :: ALBEDO_SNOWFALL_EFOLD_MASS = 3.0_wp   ! [kg m-2]

    public :: constant_surface_albedo
    public :: surface_liquid_water_content
    public :: albedo_refresh_from_snowfall
    public :: albedo_update

contains

    pure function constant_surface_albedo(mass,temperature,n,c) result(alb)
        ! Chion.jl/src/processes/albedo.jl:11-22.
        ! Memoryless: depends only on the current surface state, never on the
        ! previous albedo.

        implicit none

        real(wp),                intent(IN) :: mass(:)         ! (Ntot) [kg m-2] solid
        real(wp),                intent(IN) :: temperature(:)  ! (Ntot) [K]
        integer,                 intent(IN) :: n               ! active layers
        type(chion_const_class), intent(IN) :: c
        real(wp) :: alb

        ! No surface snow -> bare ice. Threshold TOL_EMPTY_LAYER.
        alb = c%alpha_ice
        if (n .le. 0) return
        if (mass(1) .le. TOL_EMPTY_LAYER) return

        if (temperature(1) .ge. c%T0) then
            alb = c%alpha_wet
        else
            alb = c%alpha_dry
        end if

        return

    end function constant_surface_albedo

    pure function surface_liquid_water_content(mass,mass_w,density,n,c) result(lwc)
        ! Chion.jl/src/processes/albedo.jl:30-53.
        !
        ! The pore-volume and LWC arithmetic itself is snow_column_utils'
        ! layer_lwc (computed in wp_acc; see docs/PLAN.md section 3.1). What is
        ! local to the albedo call site are its three extra guards:
        !   n <= 0 or mass(1) <= TOL_TINY        -> 0   (note: TOL_TINY, not
        !                                                TOL_EMPTY_LAYER)
        !   density(1) <= TOL_TINY               -> 0
        !   density(1) >= rho_i - TOL_TINY       -> 0   (ice-dense surface)
        ! The percolation call site guards differently; do not merge them.
        !
        ! Not clamped above 1: the wetness law below relies on lwc being able
        ! to exceed max_lwc_albedo, and clamps afterwards.

        implicit none

        real(wp),                intent(IN) :: mass(:)      ! (Ntot) [kg m-2] solid
        real(wp),                intent(IN) :: mass_w(:)    ! (Ntot) [kg m-2] liquid
        real(wp),                intent(IN) :: density(:)   ! (Ntot) [kg m-3]
        integer,                 intent(IN) :: n
        type(chion_const_class), intent(IN) :: c
        real(wp_acc) :: lwc

        lwc = 0.0_wp_acc

        if (n .le. 0) return
        if (real(mass(1),wp_acc) .le. TOL_TINY) return
        if (real(density(1),wp_acc) .le. TOL_TINY) return
        if (real(density(1),wp_acc) .ge. real(c%rho_i,wp_acc) - TOL_TINY) return

        lwc = layer_lwc(mass(1),mass_w(1),density(1),c%rho_i,c%rho_w)

        return

    end function surface_liquid_water_content

    subroutine albedo_refresh_from_snowfall(albedo,c,snowfall_mass)
        ! Chion.jl/src/processes/albedo.jl:61-81.
        ! Brightening by fresh snow, applied inside the accumulation step.
        !
        !   alpha <- min(alpha_dry, alpha + (alpha_dry-alpha_wet)*(1-exp(-dm/3)))
        !
        ! No-op below TOL_TINY of added mass. Under the CONSTANT scheme the
        ! albedo is simply reset to alpha_dry (and then recomputed from scratch
        ! by albedo_update, since the constant scheme is memoryless).
        !
        ! Trap 9: the test is on the CONSTANT scheme only, so PRESCRIBED lands
        ! in the dynamic branch here, exactly as in Julia.

        implicit none

        real(wp),                intent(INOUT) :: albedo         ! [1] column albedo
        type(chion_const_class), intent(IN)    :: c
        real(wp),                intent(IN)    :: snowfall_mass  ! [kg m-2] added this step

        if (real(snowfall_mass,wp_acc) .le. TOL_TINY) return

        if (c%albedo_scheme .eq. CHION_ALBEDO_CONSTANT) then
            albedo = c%alpha_dry
            return
        end if

        albedo = min(c%alpha_dry, albedo + (c%alpha_dry - c%alpha_wet) &
                     *(1.0_wp - exp(-snowfall_mass/ALBEDO_SNOWFALL_EFOLD_MASS)))

        return

    end subroutine albedo_refresh_from_snowfall

    subroutine albedo_update(mass,mass_w,density,temperature,n,c,albedo)
        ! Chion.jl/src/processes/albedo.jl:89-129.
        !
        ! Order of operations, all applied once per CALL (see trap 5):
        !   0. no surface snow (n<=0 or mass(1) <= TOL_EMPTY_LAYER)
        !                                       -> alpha_ice, return
        !   1. constant scheme                  -> constant_surface_albedo, return
        !   2. a = clamp(a_prev, alpha_wet, alpha_dry)
        !   3. aging   a = min(a, a - (1.35e-3*(Ts-T0) + 0.0278))
        !   4. floor   a = max(a, alpha_wet)
        !   5. wetness a = max(alpha_wet, min(a, a - (a-alpha_wet)*lwc/max_lwc))
        !   6. clamp   a = clamp(a, alpha_wet, alpha_dry)
        !
        ! Step 3's min() makes the law non-brightening whenever the bracket is
        ! non-negative, i.e. for Ts >= T0 - 0.0278/1.35e-3 ~ T0 - 20.6 K. Above
        ! that offset the bracket turns negative and the min() clamps the result
        ! to a_prev, so aging can never brighten. Preserved as written.

        implicit none

        real(wp),                intent(IN)    :: mass(:)         ! (Ntot) [kg m-2]
        real(wp),                intent(IN)    :: mass_w(:)       ! (Ntot) [kg m-2]
        real(wp),                intent(IN)    :: density(:)      ! (Ntot) [kg m-3]
        real(wp),                intent(IN)    :: temperature(:)  ! (Ntot) [K]
        integer,                 intent(IN)    :: n
        type(chion_const_class), intent(IN)    :: c
        real(wp),                intent(INOUT) :: albedo          ! [1]

        ! Local variables
        real(wp)     :: alb, t_srf
        real(wp_acc) :: lwc, alb_wet_adj

        ! Step 0: bare surface. Threshold TOL_EMPTY_LAYER (not TOL_TINY).
        if (n .le. 0) then
            albedo = c%alpha_ice
            return
        end if
        if (mass(1) .le. TOL_EMPTY_LAYER) then
            albedo = c%alpha_ice
            return
        end if

        ! Step 1: constant scheme is memoryless and ignores albedo entirely.
        if (c%albedo_scheme .eq. CHION_ALBEDO_CONSTANT) then
            albedo = constant_surface_albedo(mass,temperature,n,c)
            return
        end if

        ! Step 2: the incoming albedo may be out of range (e.g. alpha_ice left
        ! behind by a bare step), so clamp before aging.
        alb = min(max(albedo,c%alpha_wet),c%alpha_dry)

        ! Step 3: aging. NO dt -- see trap 5.
        t_srf = temperature(1)
        alb   = min(alb, alb - ((t_srf - c%T0)*ALBEDO_AGING_TEMP_COEFF + ALBEDO_AGING_OFFSET))

        ! Step 4
        alb = max(alb,c%alpha_wet)

        ! Step 5: wetness relaxation towards alpha_wet.
        lwc = surface_liquid_water_content(mass,mass_w,density,n,c)

        if (lwc .gt. 0.0_wp_acc .and. real(c%max_lwc_albedo,wp_acc) .gt. TOL_TINY) then
            alb_wet_adj = real(alb,wp_acc) &
                        - (real(alb,wp_acc) - real(c%alpha_wet,wp_acc)) &
                          *(lwc/real(c%max_lwc_albedo,wp_acc))
            alb = max(c%alpha_wet, min(alb,real(alb_wet_adj,wp)))
        end if

        ! Step 6
        albedo = min(max(alb,c%alpha_wet),c%alpha_dry)

        return

    end subroutine albedo_update

end module snow_albedo
