module snow_albedo
    ! Surface albedo: constant, dynamic, aging and prescribed schemes.
    !
    ! Port of Chion.jl/src/processes/albedo.jl.
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n, per docs/porting_notes.md D8. The per-column albedo scalar is passed
    ! as intent(INOUT), mirroring Julia's in-place mutation of state.albedo.
    !
    ! The dynamic aging decrement is scaled by dt_days (Chion.jl dev_nils
    ! 6d077c5), so diurnal substeps age by their fraction of the day and a
    ! daily step is unchanged (x*1 is exact). Before that the law decayed once
    ! per call (former trap 5, upstream defect 19).
    !
    ! PRESERVED QUIRKS (docs/PLAN.md section 5):
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
    !
    ! THIN SNOW (chion, docs/porting_notes.md D40): the schemes here compute
    ! the albedo of the SNOW (bessi's albedo_snow). What the surface energy
    ! balance sees is that blended with the background by a snow-cover
    ! fraction, thin_snow_albedo(snow_cover_fraction(SWE), alpha_snow,
    ! alpha_bg), assembled by the caller (snow_bessi).

    use chion_defs, only : wp, wp_acc, TOL_TINY, TOL_EMPTY_LAYER, &
                           chion_const_class, CHION_ALBEDO_CONSTANT, CHION_ALBEDO_AGING, &
                           ALBEDO_AGING_BINARY_REFRESH
    use snow_column_utils, only : layer_lwc

    implicit none

    private

    ! Dynamic aging law coefficients (Chion.jl albedo.jl:114). Bare magic
    ! numbers upstream; named here per docs/PLAN.md section 4.1.
    real(wp), parameter, public :: ALBEDO_AGING_TEMP_COEFF = 1.35e-3_wp  ! [K-1 d-1]
    real(wp), parameter, public :: ALBEDO_AGING_OFFSET     = 0.0278_wp   ! [d-1]

    ! Snowfall brightening e-folding mass (Chion.jl albedo.jl:77).
    real(wp), parameter, public :: ALBEDO_SNOWFALL_EFOLD_MASS = 3.0_wp   ! [kg m-2]

    public :: constant_surface_albedo
    public :: surface_liquid_water_content
    public :: albedo_refresh_from_snowfall
    public :: albedo_aging_rejuvenate
    public :: albedo_update
    public :: albedo_update_aging
    public :: snow_cover_fraction
    public :: thin_snow_albedo

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

    subroutine albedo_refresh_from_snowfall(albedo,snow_age_days,c,snowfall_mass,onto_bare)
        ! Chion.jl/src/processes/albedo.jl:61-81.
        ! Brightening by fresh snow, applied inside the accumulation step.
        !
        !   alpha <- min(alpha_dry, alpha + (alpha_dry-alpha_wet)*(1-exp(-dm/3)))
        !
        ! No-op below TOL_TINY of added mass. Under the CONSTANT scheme the
        ! albedo is simply reset to alpha_dry (and then recomputed from scratch
        ! by albedo_update, since the constant scheme is memoryless). Under the
        ! AGING scheme chion rejuvenates aged snow in proportion to the
        ! snowfall (albedo_aging_rejuvenate, D30); snow onto a bare surface
        ! (onto_bare) is all fresh, alpha_dry and age 0, however little falls:
        ! its thinness is the snow-cover fraction's business (D40). Chion.jl,
        ! and legacy_chion builds, reset to alpha_dry on any snowfall
        ! (6d06af6) and leave the snow age to albedo_update_aging.
        !
        ! Trap 9: the test is on the CONSTANT scheme only, so PRESCRIBED lands
        ! in the dynamic branch here, exactly as in Julia.

        implicit none

        real(wp),                intent(INOUT) :: albedo         ! [1] column albedo
        real(wp),                intent(INOUT) :: snow_age_days  ! [d] aging scheme only
        type(chion_const_class), intent(IN)    :: c
        real(wp),                intent(IN)    :: snowfall_mass  ! [kg m-2] added this step
        logical,                 intent(IN)    :: onto_bare      ! no surface snow before it

        if (real(snowfall_mass,wp_acc) .le. TOL_TINY) return

        if (c%albedo_scheme .eq. CHION_ALBEDO_CONSTANT) then
            albedo = c%alpha_dry
            return
        end if

        if (c%albedo_scheme .eq. CHION_ALBEDO_AGING) then
            if (ALBEDO_AGING_BINARY_REFRESH) then
                albedo = c%alpha_dry
            else if (onto_bare) then
                albedo        = c%alpha_dry
                snow_age_days = 0.0_wp
            else
                call albedo_aging_rejuvenate(albedo,snow_age_days,c,snowfall_mass)
            end if
            return
        end if

        albedo = min(c%alpha_dry, albedo + (c%alpha_dry - c%alpha_wet) &
                     *(1.0_wp - exp(-snowfall_mass/ALBEDO_SNOWFALL_EFOLD_MASS)))

        return

    end subroutine albedo_refresh_from_snowfall

    subroutine albedo_aging_rejuvenate(albedo,snow_age_days,c,snowfall_mass)
        ! Partial rejuvenation of aged snow by fresh snow (chion deviation,
        ! docs/porting_notes.md D30). With the refresh fraction
        !
        !   f = 1 - exp(-S/aging_snowfall_ref),   S = snowfall_mass,
        !
        ! the aging progress E = -ln((a - alpha_wet)/(alpha_dry - alpha_wet))
        ! (= sum of dt/tau, exact for the aging scheme) and the snow age are
        ! scaled by 1 - f = exp(-S/aging_snowfall_ref). In albedo space, with
        ! no extra state:
        !
        !   a   <- alpha_wet + (alpha_dry - alpha_wet)*x**(1-f),
        !          x = (a - alpha_wet)/(alpha_dry - alpha_wet)
        !   age <- (1 - f)*age
        !
        ! S = 0 changes nothing, a trace almost nothing, and refreshes compose
        ! exactly: exp(-S1/S_ref)*exp(-S2/S_ref) = exp(-(S1+S2)/S_ref), so a
        ! day's snowfall refreshes the same whether it falls in one step or in
        ! diurnal substeps (review Q14).
        !
        ! x = 0 (E infinite, also alpha_dry = alpha_wet, where x is undefined)
        ! stays alpha_wet: no partial refresh restores it. Snow onto a bare
        ! surface does not come here (albedo_refresh_from_snowfall).

        implicit none

        real(wp),                intent(INOUT) :: albedo         ! [1]
        real(wp),                intent(INOUT) :: snow_age_days  ! [d]
        type(chion_const_class), intent(IN)    :: c
        real(wp),                intent(IN)    :: snowfall_mass  ! [kg m-2] added this step

        ! Local variables
        real(wp) :: f_keep, alb, span

        ! 1 - f: the fraction of the aging progress that survives.
        f_keep = exp(-max(snowfall_mass,0.0_wp)/c%aging_snowfall_ref)

        snow_age_days = f_keep*snow_age_days

        alb  = min(max(albedo,c%alpha_wet),c%alpha_dry)
        span = c%alpha_dry - c%alpha_wet

        if (alb .gt. c%alpha_wet) then
            albedo = c%alpha_wet + span*((alb - c%alpha_wet)/span)**f_keep
        else
            albedo = c%alpha_wet
        end if

        return

    end subroutine albedo_aging_rejuvenate

    subroutine albedo_update(mass,mass_w,density,temperature,n,c,dt_days,albedo)
        ! Chion.jl/src/processes/albedo.jl:89-129 (dev_nils 6d077c5).
        !
        ! Order of operations:
        !   0. no surface snow (n<=0 or mass(1) <= TOL_EMPTY_LAYER)
        !                                       -> alpha_ice, return
        !   1. constant scheme                  -> constant_surface_albedo, return
        !   2. a = clamp(a_prev, alpha_wet, alpha_dry)
        !   3. aging   a = min(a, a - (1.35e-3*(Ts-T0) + 0.0278)*dt_days)
        !   4. floor   a = max(a, alpha_wet)
        !   5. wetness, if dt_days > 0 (Chion.jl 03bb445):
        !        r = clamp(lwc/max_lwc, 0, 1)
        !        a = max(alpha_wet, min(a, alpha_wet + (a-alpha_wet)*(1-r)**dt_days))
        !      The one-day relaxation of the former linear law
        !      a - (a-alpha_wet)*r, but composing exactly across substeps at
        !      fixed wetness; saturated snow reaches alpha_wet directly.
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
        real(wp),                intent(IN)    :: dt_days         ! [d] step length
        real(wp),                intent(INOUT) :: albedo          ! [1]

        ! Local variables
        real(wp)     :: alb, t_srf
        real(wp_acc) :: lwc, wet_fraction, retention, alb_wet_adj

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

        ! Step 3: aging, scaled by the step length.
        t_srf = temperature(1)
        alb   = min(alb, alb - ((t_srf - c%T0)*ALBEDO_AGING_TEMP_COEFF + ALBEDO_AGING_OFFSET)*dt_days)

        ! Step 4
        alb = max(alb,c%alpha_wet)

        ! Step 5: wetness relaxation towards alpha_wet, dt-consistent.
        lwc = surface_liquid_water_content(mass,mass_w,density,n,c)

        if (dt_days .gt. 0.0_wp .and. lwc .gt. 0.0_wp_acc .and. &
            real(c%max_lwc_albedo,wp_acc) .gt. TOL_TINY) then
            wet_fraction = min(max(lwc/real(c%max_lwc_albedo,wp_acc),0.0_wp_acc),1.0_wp_acc)
            retention    = (1.0_wp_acc - wet_fraction)**real(dt_days,wp_acc)
            alb_wet_adj  = real(c%alpha_wet,wp_acc) &
                         + (real(alb,wp_acc) - real(c%alpha_wet,wp_acc))*retention
            alb = max(c%alpha_wet, min(alb,real(alb_wet_adj,wp)))
        end if

        ! Step 6
        albedo = min(max(alb,c%alpha_wet),c%alpha_dry)

        return

    end subroutine albedo_update

    subroutine albedo_update_aging(mass,temperature,n,c,snowfall_rate,dt_days, &
                                   albedo,snow_age_days)
        ! Chion.jl/src/processes/albedo.jl _update_aging_surface_albedo_arrays!
        ! (6d06af6). Snow albedo from the time since the latest snowfall:
        !
        !   0. no surface snow                  -> alpha_ice, age 0, return
        !   1. snowfall this step (rate > 0)    -> alpha_dry, age 0, return
        !   2. age = max(age, 0) + dt_days
        !   3. a = alpha_wet + (clamp(a_prev, alpha_wet, alpha_dry) - alpha_wet)
        !          *exp(-dt_days/tau)
        !      tau = aging_melting_timescale_days if Ts >= T0,
        !            aging_cold_timescale_days    otherwise
        !   4. clamp   a = clamp(a, alpha_wet, alpha_dry)
        !
        ! Step 1 tests the RATE, so any snowfall at all fully rejuvenates the
        ! surface, however little mass it adds. That is Chion.jl's form, kept
        ! under legacy_chion only (ALBEDO_AGING_BINARY_REFRESH). chion instead
        ! rejuvenates in proportion to the snowfall in the accumulation step
        ! (albedo_aging_rejuvenate, D30), so step 1 is skipped and fresh snow
        ! ages from the step it fell. Liquid water plays no role.

        implicit none

        real(wp),                intent(IN)    :: mass(:)         ! (Ntot) [kg m-2]
        real(wp),                intent(IN)    :: temperature(:)  ! (Ntot) [K]
        integer,                 intent(IN)    :: n
        type(chion_const_class), intent(IN)    :: c
        real(wp),                intent(IN)    :: snowfall_rate   ! [kg m-2 s-1]
        real(wp),                intent(IN)    :: dt_days         ! [d] step length
        real(wp),                intent(INOUT) :: albedo          ! [1]
        real(wp),                intent(INOUT) :: snow_age_days   ! [d]

        ! Local variables
        real(wp) :: alb, tau

        ! Step 0: bare surface. Threshold TOL_EMPTY_LAYER, as albedo_update.
        if (n .le. 0) then
            albedo        = c%alpha_ice
            snow_age_days = 0.0_wp
            return
        end if
        if (mass(1) .le. TOL_EMPTY_LAYER) then
            albedo        = c%alpha_ice
            snow_age_days = 0.0_wp
            return
        end if

        ! Step 1: fresh snow (Chion.jl's binary reset; legacy_chion only).
        if (ALBEDO_AGING_BINARY_REFRESH) then
            if (snowfall_rate .gt. 0.0_wp) then
                albedo        = c%alpha_dry
                snow_age_days = 0.0_wp
                return
            end if
        end if

        ! Step 2
        snow_age_days = max(snow_age_days,0.0_wp) + dt_days

        ! Step 3: exponential relaxation towards alpha_wet.
        if (temperature(1) .ge. c%T0) then
            tau = c%aging_melting_timescale_days
        else
            tau = c%aging_cold_timescale_days
        end if

        alb = min(max(albedo,c%alpha_wet),c%alpha_dry)
        alb = c%alpha_wet + (alb - c%alpha_wet)*exp(-dt_days/tau)

        ! Step 4
        albedo = min(max(alb,c%alpha_wet),c%alpha_dry)

        return

    end subroutine albedo_update_aging

    pure function snow_cover_fraction(swe,swe_crit) result(f)
        ! Snow-cover fraction of the thin-snow albedo (chion, D40):
        !
        !   f = min(1, SWE/swe_crit)     swe_crit > 0
        !   f = 1                        swe_crit = 0 (blend off)
        !
        ! SWE is the column's snow water equivalent [kg m-2], solid plus
        ! liquid over the snow layers (not the ice substrate): a column of
        ! firn covers fully, fresh snow on bare ice partly. Whether there is
        ! surface snow at all (f = 0 without) is the caller's test.

        implicit none

        real(wp), intent(IN) :: swe        ! [kg m-2]
        real(wp), intent(IN) :: swe_crit   ! [kg m-2] >= 0
        real(wp) :: f

        if (swe_crit .gt. 0.0_wp) then
            f = min(1.0_wp, max(swe,0.0_wp)/swe_crit)
        else
            f = 1.0_wp
        end if

        return

    end function snow_cover_fraction

    pure function thin_snow_albedo(f,alb_snow,alb_bg) result(alb)
        ! The albedo the surface energy balance sees (D40):
        !
        !   alpha = f*alpha_snow + (1 - f)*alpha_bg
        !
        ! Written as two products so f = 1 returns alpha_snow and f = 0
        ! alpha_bg exactly.

        implicit none

        real(wp), intent(IN) :: f          ! [1] snow-cover fraction
        real(wp), intent(IN) :: alb_snow   ! [1] snow albedo
        real(wp), intent(IN) :: alb_bg     ! [1] background albedo
        real(wp) :: alb

        alb = f*alb_snow + (1.0_wp - f)*alb_bg

        return

    end function thin_snow_albedo

end module snow_albedo
