module snow_layers
    ! Layer-structure updates and basal depletion for a single snowpack column.
    !
    ! Port of Chion.jl/src/processes/layer_structure.jl (branch main), in full:
    !     _reset_layer_at_index!                  -> reset_layer_at_index
    !     _split_surface_layer!                   -> split_surface_layer
    !     _merge_surface_layer!                   -> merge_surface_layer
    !     _merge_bottom_layer!                    -> merge_bottom_layer
    !     _remove_surface_layer!                  -> remove_surface_layer
    !     _remove_depleted_surface_and_route_water! -> remove_depleted_surface_and_route_water
    !     _continuous_bottom_deplete!             -> continuous_bottom_deplete
    !     _free_slot_for_surface_split!           -> free_slot_for_surface_split
    !     _enforce_snow_depth_cap!                -> enforce_snow_depth_cap
    !
    ! CALLING CONVENTION (docs/porting_notes.md D8): every routine takes
    ! contiguous column slices -- mass(:,icol), mass_w(:,icol), ... -- together
    ! with the active layer count n. Routines that change the layer count take
    ! n as intent(INOUT). The Julia (array, idx) indirection exists only for GPU
    ! dispatch and has no Fortran analogue.
    !
    ! PRECISION (docs/PLAN.md section 3.1, porting_notes.md D1):
    !   - mass, mass_w, density, temperature, t_srf, albedo are wp (sp).
    !   - mass_base, smb_ice, runoff are per-column CUMULATIVE accumulators and
    !     are therefore wp_acc (dp) scalars.
    !   - Every sum over layers, every running remainder, and every difference
    !     of near-equal numbers (the depth accumulation, the depletion
    !     remainder d_m, the rho_i excess export) is carried in wp_acc locals
    !     and converted back to wp only on store.
    !
    ! SHORT CIRCUITING (docs/porting_notes.md D10): Julia's `||`/`&&`
    ! short-circuit; Fortran's `.or.`/`.and.` do not. Two guards in this file
    ! rely on that in Julia and are written as nested ifs here:
    !   - the tail-trim loop of continuous_bottom_deplete
    !     (`n > 0 && mass(n) <= EPS_EMPTY_LAYER`), which would read mass(0);
    !   - the split guard is safe either way but is left as a single .and.
    !     because both operands are unconditionally evaluable.
    !
    ! THRESHOLDS: the three "empty" tests are used deliberately and differently
    ! within this file, exactly as in Julia. Do not unify them
    ! (docs/PLAN.md section 5, item 1):
    !   > 0                free_slot_for_surface_split bottom-mass test,
    !                      depth-cap layer inclusion, mass-weighted-mean guard
    !   > TOL_TINY         depletion remainder loop, depth-cap density guard
    !   < TOL_EMPTY_LAYER  merge_surface_layer single-layer collapse,
    !                      continuous_bottom_deplete empty-layer skip and trim
    !
    ! DEPTH CAP: enforce_snow_depth_cap uses the hard-coded
    ! BESSI_REFERENCE_LAYER_COUNT = 15 and BESSI_REFERENCE_DEPTH_DENSITY = 300
    ! with a 1.5 factor, INDEPENDENT of the configured Ntot. Preserved
    ! deliberately -- see docs/PLAN.md section 5, item 11. Making the cap
    ! respect Ntot is listed in section 4.1 as "not allowed without asking".

    use chion_defs, only : wp, wp_acc, io_unit_err, TOL_TINY, TOL_EMPTY_LAYER, &
                           BESSI_REFERENCE_LAYER_COUNT, BESSI_REFERENCE_DEPTH_DENSITY, &
                           chion_const_class

    implicit none

    private

    public :: reset_layer_at_index
    public :: split_surface_layer
    public :: merge_surface_layer
    public :: merge_bottom_layer
    public :: remove_surface_layer
    public :: remove_depleted_surface_and_route_water
    public :: continuous_bottom_deplete
    public :: free_slot_for_surface_split
    public :: enforce_snow_depth_cap

contains

    ! =====================================================================
    ! Helpers (Chion.jl layer_structure.jl:10-20)
    ! =====================================================================

    pure function safe_nonnegative(x) result(y)
        ! Chion.jl _safe_nonnegative: clamp to zero from below.

        implicit none

        real(wp_acc), intent(IN) :: x
        real(wp_acc) :: y

        if (x .gt. 0.0_wp_acc) then
            y = x
        else
            y = 0.0_wp_acc
        end if

        return

    end function safe_nonnegative

    pure function mass_weighted_mean(m1,x1,m2,x2) result(xbar)
        ! Chion.jl _mass_weighted_mean: mass-weighted mean of two layer
        ! properties, returning zero when the combined mass is non-positive.
        !
        ! Evaluated in wp_acc: this is a sum-then-divide over two layers whose
        ! masses can differ by orders of magnitude just after a split.

        implicit none

        real(wp), intent(IN) :: m1, x1
        real(wp), intent(IN) :: m2, x2
        real(wp) :: xbar

        ! Local variables
        real(wp_acc) :: total_mass

        total_mass = real(m1,wp_acc) + real(m2,wp_acc)

        if (total_mass .gt. 0.0_wp_acc) then
            xbar = real((real(m1,wp_acc)*real(x1,wp_acc) &
                       + real(m2,wp_acc)*real(x2,wp_acc)) / total_mass, wp)
        else
            xbar = 0.0_wp
        end if

        return

    end function mass_weighted_mean

    ! =====================================================================
    ! reset_layer_at_index
    ! Chion.jl layer_structure.jl:28-42
    ! =====================================================================

    subroutine reset_layer_at_index(mass,mass_w,density,temperature,k,c)
        ! Reset one layer to an empty state. Note the temperature is reset to
        ! c%T0, NOT to zero -- an empty layer still carries the freezing point
        ! so that later mass-weighted merges stay well-behaved.

        implicit none

        real(wp), intent(INOUT) :: mass(:)          ! (Ntot) [kg m-2] solid
        real(wp), intent(INOUT) :: mass_w(:)        ! (Ntot) [kg m-2] liquid
        real(wp), intent(INOUT) :: density(:)       ! (Ntot) [kg m-3]
        real(wp), intent(INOUT) :: temperature(:)   ! (Ntot) [K]
        integer,  intent(IN)    :: k                ! layer index to reset
        type(chion_const_class), intent(IN) :: c

        mass(k)        = 0.0_wp
        mass_w(k)      = 0.0_wp
        density(k)     = 0.0_wp
        temperature(k) = c%T0

        return

    end subroutine reset_layer_at_index

    ! =====================================================================
    ! split_surface_layer
    ! Chion.jl layer_structure.jl:51-91
    ! =====================================================================

    subroutine split_surface_layer(mass,mass_w,density,temperature,n,Ntot,mass_max,mass_split)
        ! Split the surface layer in two when it exceeds mass_max, provided a
        ! free slot exists (n < Ntot). Layer 2 receives exactly mass_split and
        ! layer 1 keeps the remainder; density and temperature are copied
        ! unchanged into both halves, and liquid water is partitioned by the
        ! solid-mass fraction.
        !
        ! No-ops unless BOTH conditions hold. The caller is responsible for
        ! freeing a slot first (free_slot_for_surface_split) when n == Ntot.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n                ! active layer count
        integer,  intent(IN)    :: Ntot             ! maximum layer count
        real(wp), intent(IN)    :: mass_max         ! [kg m-2] split trigger
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] target layer mass

        ! Local variables
        integer  :: k, new_n
        real(wp) :: surface_mass, surface_mass_w
        real(wp) :: surface_density, surface_temperature
        real(wp) :: water_fraction

        surface_mass = mass(1)

        ! Julia: if !(n < Ntot && surface_mass > mass_max) return
        ! Both operands are safe to evaluate, so a single .and. is fine here.
        if (.not. (n .lt. Ntot .and. surface_mass .gt. mass_max)) return

        surface_mass_w      = mass_w(1)
        surface_density     = density(1)
        surface_temperature = temperature(1)

        new_n = n + 1
        n     = new_n

        ! Shift layers 2..new_n-1 down by one, leaving slot 2 free.
        do k = new_n, 3, -1
            mass(k)        = mass(k-1)
            mass_w(k)      = mass_w(k-1)
            density(k)     = density(k-1)
            temperature(k) = temperature(k-1)
        end do

        mass(2) = mass_split
        mass(1) = surface_mass - mass_split

        ! surface_mass > mass_max >= 0, so this division is safe.
        water_fraction = mass_split / surface_mass
        mass_w(2)      = surface_mass_w * water_fraction
        mass_w(1)      = surface_mass_w * (1.0_wp - water_fraction)

        density(1)     = surface_density
        density(2)     = surface_density
        temperature(1) = surface_temperature
        temperature(2) = surface_temperature

        return

    end subroutine split_surface_layer

    ! =====================================================================
    ! merge_surface_layer
    ! Chion.jl layer_structure.jl:99-172
    ! =====================================================================

    subroutine merge_surface_layer(mass,mass_w,density,temperature,n,mass_split,mass_min,c)
        ! Rebalance or merge the top two layers when the surface layer falls
        ! below mass_min. Three outcomes:
        !   1. n == 1 and the surface layer is essentially empty
        !      (< TOL_EMPTY_LAYER): the column collapses to zero layers.
        !   2. combined mass > 2*mass_split: PARTIAL TRANSFER. The surface layer
        !      is topped back up to exactly mass_split from layer 2; the layer
        !      count is unchanged.
        !   3. otherwise: FULL MERGE. Layers 1 and 2 are combined into layer 1
        !      and everything below shifts up by one.
        !
        ! PORTING NOTE: the Julia signature also carries Ntot, which its body
        ! never uses. Dropped here rather than carried as an unused dummy.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] target layer mass
        real(wp), intent(IN)    :: mass_min         ! [kg m-2] merge trigger
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer  :: k, new_n
        real(wp) :: surface_mass, subsurface_mass
        real(wp) :: transferred_to_surface, transferred_water
        real(wp) :: dens_1, dens_2, temp_1, temp_2
        real(wp_acc) :: combined_mass

        if (n .gt. 0) then
            surface_mass = mass(1)
        else
            surface_mass = 0.0_wp
        end if

        ! Threshold here is TOL_EMPTY_LAYER (strict <), not TOL_TINY.
        if (n .eq. 1 .and. real(surface_mass,wp_acc) .lt. TOL_EMPTY_LAYER) then
            n = 0
            call reset_layer_at_index(mass,mass_w,density,temperature,1,c)
            return
        else if (n .le. 1 .or. surface_mass .ge. mass_min) then
            return
        end if

        subsurface_mass = mass(2)
        combined_mass   = real(surface_mass,wp_acc) + real(subsurface_mass,wp_acc)

        dens_1 = density(1)
        dens_2 = density(2)
        temp_1 = temperature(1)
        temp_2 = temperature(2)

        if (combined_mass .gt. 2.0_wp_acc*real(mass_split,wp_acc)) then

            ! --- Partial transfer: top the surface layer back up to mass_split.
            ! subsurface_mass = combined - surface > 2*mass_split - mass_min,
            ! which is strictly positive for any sane parameter set, so the
            ! division below is safe. Julia does not guard it either.
            transferred_to_surface = mass_split - surface_mass
            transferred_water      = transferred_to_surface / subsurface_mass * mass_w(2)

            mass(1)   = mass_split
            mass(2)   = real(combined_mass - real(mass_split,wp_acc),wp)
            mass_w(1) = mass_w(1) + transferred_water
            mass_w(2) = mass_w(2) - transferred_water

            density(1)     = mass_weighted_mean(surface_mass,dens_1, &
                                                transferred_to_surface,dens_2)
            temperature(1) = mass_weighted_mean(surface_mass,temp_1, &
                                                transferred_to_surface,temp_2)
            return

        end if

        ! --- Full merge of layers 1 and 2.
        mass(1)   = real(combined_mass,wp)
        mass_w(1) = mass_w(1) + mass_w(2)

        density(1)     = mass_weighted_mean(surface_mass,dens_1,subsurface_mass,dens_2)
        temperature(1) = mass_weighted_mean(surface_mass,temp_1,subsurface_mass,temp_2)

        new_n = n - 1
        n     = new_n

        do k = 2, new_n
            mass(k)        = mass(k+1)
            mass_w(k)      = mass_w(k+1)
            density(k)     = density(k+1)
            temperature(k) = temperature(k+1)
        end do

        call reset_layer_at_index(mass,mass_w,density,temperature,new_n+1,c)

        return

    end subroutine merge_surface_layer

    ! =====================================================================
    ! merge_bottom_layer
    ! Chion.jl layer_structure.jl:180-232
    ! =====================================================================

    subroutine merge_bottom_layer(mass,mass_w,density,temperature,n, &
                                  mass_base,smb_ice,c)
        ! Merge the two deepest active layers, freeing one slot.
        !
        ! If the combined density would exceed pure ice, the density is pinned
        ! at rho_i and the mass in excess of what rho_i can hold at that
        ! thickness is EXPORTED to the basal accumulators (mass_base and
        ! smb_ice) rather than being kept in the column. That export is the
        ! only mass sink in this routine, and it is what keeps the column mass
        ! budget closed: mass leaving the layers reappears in mass_base.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp_acc), intent(INOUT) :: mass_base    ! [kg m-2] cumulative basal mass
        real(wp_acc), intent(INOUT) :: smb_ice      ! [kg m-2] cumulative ice SMB
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer      :: n_old, n_new
        real(wp)     :: combined_density, combined_temperature
        real(wp_acc) :: combined_mass, combined_mass_w
        real(wp_acc) :: mass_limited_to_ice_density, exported_excess

        if (n .lt. 2) return

        n_old = n
        n_new = n - 1
        n     = n_new

        combined_mass   = real(mass(n_new),wp_acc)   + real(mass(n_old),wp_acc)
        combined_mass_w = real(mass_w(n_new),wp_acc) + real(mass_w(n_old),wp_acc)

        combined_density     = mass_weighted_mean(mass(n_new),density(n_new), &
                                                  mass(n_old),density(n_old))
        combined_temperature = mass_weighted_mean(mass(n_new),temperature(n_new), &
                                                  mass(n_old),temperature(n_old))

        if (combined_density .gt. c%rho_i) then

            ! Difference of near-equal numbers -- carried in wp_acc.
            mass_limited_to_ice_density = combined_mass &
                                        * (real(c%rho_i,wp_acc)/real(combined_density,wp_acc))
            exported_excess = combined_mass - mass_limited_to_ice_density

            mass_base = mass_base + exported_excess
            smb_ice   = smb_ice   + exported_excess

            mass(n_new)        = real(mass_limited_to_ice_density,wp)
            mass_w(n_new)      = real(combined_mass_w,wp)
            density(n_new)     = c%rho_i
            temperature(n_new) = combined_temperature

        else

            mass(n_new)        = real(combined_mass,wp)
            mass_w(n_new)      = real(combined_mass_w,wp)
            density(n_new)     = combined_density
            temperature(n_new) = combined_temperature

        end if

        call reset_layer_at_index(mass,mass_w,density,temperature,n_old,c)

        return

    end subroutine merge_bottom_layer

    ! =====================================================================
    ! remove_surface_layer
    ! Chion.jl layer_structure.jl:240-267
    ! =====================================================================

    subroutine remove_surface_layer(mass,mass_w,density,temperature,n,c)
        ! Drop the top active layer and shift everything below it up by one.
        ! Any mass still held in layer 1 is DISCARDED -- callers must have
        ! emptied it first (see remove_depleted_surface_and_route_water).

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer :: k

        if (n .le. 0) then
            return
        else if (n .eq. 1) then
            call reset_layer_at_index(mass,mass_w,density,temperature,1,c)
            n = 0
            return
        end if

        do k = 1, n-1
            mass(k)        = mass(k+1)
            mass_w(k)      = mass_w(k+1)
            density(k)     = density(k+1)
            temperature(k) = temperature(k+1)
        end do

        call reset_layer_at_index(mass,mass_w,density,temperature,n,c)
        n = n - 1

        return

    end subroutine remove_surface_layer

    ! =====================================================================
    ! remove_depleted_surface_and_route_water
    ! Chion.jl layer_structure.jl:275-294
    ! =====================================================================

    subroutine remove_depleted_surface_and_route_water(mass,mass_w,density,temperature,n, &
                                                       runoff,c)
        ! Remove an exhausted surface layer, first moving its residual liquid
        ! water into layer 2 if one exists, or to runoff if it does not.
        !
        ! Julia does not guard n == 0 here: with no layers, mass_w(1) (which is
        ! zero after any reset) is added to runoff and remove_surface_layer
        ! no-ops. Preserved.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp_acc), intent(INOUT) :: runoff       ! [kg m-2] cumulative runoff
        type(chion_const_class), intent(IN) :: c

        if (n .gt. 1) then
            mass_w(2) = mass_w(2) + mass_w(1)
        else
            runoff = runoff + real(mass_w(1),wp_acc)
        end if

        mass_w(1) = 0.0_wp

        call remove_surface_layer(mass,mass_w,density,temperature,n,c)

        return

    end subroutine remove_depleted_surface_and_route_water

    ! =====================================================================
    ! continuous_bottom_deplete
    ! Chion.jl layer_structure.jl:302-366
    ! =====================================================================

    subroutine continuous_bottom_deplete(mass,mass_w,density,temperature,n, &
                                         mass_base,smb_ice,runoff,t_srf,albedo, &
                                         d_m_in,c,ice_to_base,runoff_out)
        ! Remove d_m_in of solid mass from the BOTTOM of the column, consuming
        ! whole layers until the remainder fits inside one. Removed solid mass
        ! goes to mass_base and smb_ice; the liquid water that travelled with it
        ! goes to runoff. This is the routine through which all basal export
        ! passes (depth cap, slot freeing, and the public wrapper).
        !
        ! d_m_in is wp_acc: it is a running remainder repeatedly reduced by
        ! layer masses of comparable magnitude, and the loop guard tests it
        ! against TOL_TINY, which is only meaningful in dp (porting_notes D1b).
        !
        ! The routine also refreshes the two instantaneous surface diagnostics:
        ! t_srf tracks layer 1 (or T0 when the column empties) and the albedo
        ! is forced to bare-ice when the column empties.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp_acc), intent(INOUT) :: mass_base
        real(wp_acc), intent(INOUT) :: smb_ice
        real(wp_acc), intent(INOUT) :: runoff
        real(wp), intent(INOUT) :: t_srf            ! [K] surface temperature diagnostic
        real(wp), intent(INOUT) :: albedo           ! [1] dynamic albedo
        real(wp_acc), intent(IN) :: d_m_in          ! [kg m-2] solid mass to remove
        type(chion_const_class), intent(IN) :: c
        real(wp_acc), optional, intent(OUT) :: ice_to_base   ! [kg m-2] solid removed
        real(wp_acc), optional, intent(OUT) :: runoff_out    ! [kg m-2] water released

        ! Local variables
        integer      :: nn, tail
        logical      :: keep_trimming
        real(wp_acc) :: d_m, d_lw, layer_mass, runoff_mass
        real(wp_acc) :: ice_total, runoff_total

        d_m          = safe_nonnegative(d_m_in)
        ice_total    = 0.0_wp_acc
        runoff_total = 0.0_wp_acc

        do while (d_m .gt. TOL_TINY .and. n .gt. 0)

            nn         = n
            layer_mass = real(mass(nn),wp_acc)

            if (layer_mass .le. TOL_EMPTY_LAYER) then
                ! An already-empty bottom layer is simply dropped, without
                ! consuming any of the remaining demand.
                call reset_layer_at_index(mass,mass_w,density,temperature,nn,c)
                n = nn - 1
                cycle
            end if

            if (d_m .gt. layer_mass) then

                ! Consume the whole layer and carry the remainder downward.
                d_m         = d_m - layer_mass
                ice_total   = ice_total + layer_mass
                runoff_mass = real(mass_w(nn),wp_acc)

                runoff_total = runoff_total + runoff_mass
                mass_base    = mass_base + layer_mass
                smb_ice      = smb_ice   + layer_mass
                runoff       = runoff    + runoff_mass

                call reset_layer_at_index(mass,mass_w,density,temperature,nn,c)
                n = nn - 1

            else

                ! Partial consumption: liquid water leaves in proportion to the
                ! solid mass removed. Density and temperature are unchanged.
                d_lw = d_m * real(mass_w(nn),wp_acc) / layer_mass

                mass(nn)   = real(layer_mass - d_m,wp)
                mass_w(nn) = real(real(mass_w(nn),wp_acc) - d_lw,wp)

                ice_total    = ice_total + d_m
                runoff_total = runoff_total + d_lw
                mass_base    = mass_base + d_m
                smb_ice      = smb_ice   + d_m
                runoff       = runoff    + d_lw

                d_m = 0.0_wp_acc

            end if

        end do

        ! Trim any residual empty layers off the bottom.
        ! Julia relies on `&&` short-circuiting here; mass(n) with n == 0 would
        ! be out of bounds in Fortran, hence the nested form (porting_notes D10).
        keep_trimming = .TRUE.
        do while (keep_trimming)
            keep_trimming = .FALSE.
            if (n .gt. 0) then
                if (real(mass(n),wp_acc) .le. TOL_EMPTY_LAYER) then
                    tail = n
                    call reset_layer_at_index(mass,mass_w,density,temperature,tail,c)
                    n = tail - 1
                    keep_trimming = .TRUE.
                end if
            end if
        end do

        if (n .gt. 0) then
            t_srf = temperature(1)
        else
            t_srf = c%T0
        end if

        if (n .eq. 0) albedo = c%alpha_ice

        if (present(ice_to_base)) ice_to_base = ice_total
        if (present(runoff_out))  runoff_out  = runoff_total

        return

    end subroutine continuous_bottom_deplete

    ! =====================================================================
    ! free_slot_for_surface_split
    ! Chion.jl layer_structure.jl:374-434
    ! =====================================================================

    subroutine free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
                                           mass_base,smb_ice,runoff,t_srf,albedo, &
                                           Ntot,mass_max,c)
        ! Make room for a surface split when the column is already at Ntot.
        !
        ! Ntot == 1 SPECIAL CASE: there is no slot to free, so instead the
        ! surface layer's own overflow beyond mass_max is exported basally.
        ! This is the only path by which a single-layer column sheds mass at
        ! the base during accumulation.
        !
        ! Otherwise the deepest layer is depleted in full (routing its water to
        ! runoff), or simply reset if it carries no mass.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp_acc), intent(INOUT) :: mass_base
        real(wp_acc), intent(INOUT) :: smb_ice
        real(wp_acc), intent(INOUT) :: runoff
        real(wp), intent(INOUT) :: t_srf
        real(wp), intent(INOUT) :: albedo
        integer,  intent(IN)    :: Ntot
        real(wp), intent(IN)    :: mass_max
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        real(wp_acc) :: overflow, bottom_mass

        if (Ntot .eq. 1) then

            overflow = max(real(mass(1),wp_acc) - real(mass_max,wp_acc), 0.0_wp_acc)

            if (overflow .gt. 0.0_wp_acc) then
                call continuous_bottom_deplete(mass,mass_w,density,temperature,n, &
                                               mass_base,smb_ice,runoff,t_srf,albedo, &
                                               overflow,c)
            end if

            return

        end if

        ! Julia indexes mass at _n_active(...) without checking it is >= 1;
        ! with n == 0 that is an out-of-bounds read. chion makes the
        ! precondition explicit (docs/PLAN.md section 4.1, "fix outright bugs").
        if (n .lt. 1) then
            write(io_unit_err,*) "free_slot_for_surface_split:: Error: &
                                 &called on a column with no active layers."
            write(io_unit_err,*) "n, Ntot = ", n, Ntot
            stop "Program stopped."
        end if

        bottom_mass = real(mass(n),wp_acc)

        ! Threshold here is exactly zero, not TOL_EMPTY_LAYER.
        if (bottom_mass .gt. 0.0_wp_acc) then
            call continuous_bottom_deplete(mass,mass_w,density,temperature,n, &
                                           mass_base,smb_ice,runoff,t_srf,albedo, &
                                           bottom_mass,c)
        else
            call reset_layer_at_index(mass,mass_w,density,temperature,n,c)
            n = n - 1
        end if

        return

    end subroutine free_slot_for_surface_split

    ! =====================================================================
    ! enforce_snow_depth_cap
    ! Chion.jl layer_structure.jl:442-513
    ! =====================================================================

    subroutine enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                      mass_base,smb_ice,runoff,t_srf,albedo, &
                                      mass_split,c)
        ! Cap the total active snow depth after accumulation, exporting the
        ! excess from the base of the column.
        !
        ! THE CAP IGNORES THE CONFIGURED Ntot. reference_depth is built from
        ! the hard-coded BESSI_REFERENCE_LAYER_COUNT = 15 and
        ! BESSI_REFERENCE_DEPTH_DENSITY = 300 with a factor 1.5, so a column
        ! configured with Ntot = 2 is still capped at the depth 15 reference
        ! layers would occupy. This is trap 11 in docs/PLAN.md section 5;
        ! changing it is explicitly listed in section 4.1 as a modelling
        ! decision that may not be taken inside a work package.
        !
        ! PORTING NOTE: the Julia signature also carries Ntot and dt_seconds,
        ! neither of which its body uses. Both dropped rather than carried as
        ! unused dummies.
        !
        ! Layer inclusion in the depth sum uses mass > 0 AND density > TOL_TINY.
        ! The reverse pass that converts an excess DEPTH back into an excess
        ! MASS treats a zero-density layer differently again: its mass is
        ! exported in full without reducing the remaining depth demand.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp_acc), intent(INOUT) :: mass_base
        real(wp_acc), intent(INOUT) :: smb_ice
        real(wp_acc), intent(INOUT) :: runoff
        real(wp), intent(INOUT) :: t_srf
        real(wp), intent(INOUT) :: albedo
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] target layer mass
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer      :: k
        real(wp_acc) :: total_active_snow_depth, reference_depth, excess_depth
        real(wp_acc) :: excess_basal_mass, remaining_excess_depth
        real(wp_acc) :: layer_mass, layer_density, layer_depth

        ! --- Total active solid depth. Sum over layers -> wp_acc.
        total_active_snow_depth = 0.0_wp_acc

        do k = 1, n
            layer_mass    = real(mass(k),wp_acc)
            layer_density = real(density(k),wp_acc)
            if (layer_mass .gt. 0.0_wp_acc .and. layer_density .gt. TOL_TINY) then
                total_active_snow_depth = total_active_snow_depth + layer_mass/layer_density
            end if
        end do

        reference_depth = real(BESSI_REFERENCE_LAYER_COUNT,wp_acc) &
                        * real(mass_split,wp_acc) * 1.5_wp_acc &
                        / real(BESSI_REFERENCE_DEPTH_DENSITY,wp_acc)

        excess_depth = total_active_snow_depth - reference_depth

        if (excess_depth .le. 0.0_wp_acc) return

        ! --- Convert the excess depth into an excess mass, from the base up.
        excess_basal_mass      = 0.0_wp_acc
        remaining_excess_depth = excess_depth

        do k = n, 1, -1

            layer_mass    = real(mass(k),wp_acc)
            layer_density = real(density(k),wp_acc)

            if (layer_mass .le. 0.0_wp_acc) then
                cycle
            else if (layer_density .le. TOL_TINY) then
                ! Massive but density-less: contributes no depth, so its mass
                ! is exported without reducing the remaining demand.
                excess_basal_mass = excess_basal_mass + layer_mass
                cycle
            end if

            layer_depth = layer_mass/layer_density

            if (layer_depth .gt. remaining_excess_depth) then
                excess_basal_mass = excess_basal_mass + remaining_excess_depth*layer_density
                exit
            end if

            excess_basal_mass      = excess_basal_mass + layer_mass
            remaining_excess_depth = remaining_excess_depth - layer_depth

        end do

        call continuous_bottom_deplete(mass,mass_w,density,temperature,n, &
                                       mass_base,smb_ice,runoff,t_srf,albedo, &
                                       excess_basal_mass,c)

        return

    end subroutine enforce_snow_depth_cap

end module snow_layers
