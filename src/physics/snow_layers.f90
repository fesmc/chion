module snow_layers
    ! Layer-structure updates and basal depletion for a single snowpack column.
    !
    ! Port of Chion.jl/src/processes/layer_structure.jl (branch main), in full:
    !     _reset_layer_at_index!                  -> reset_layer_at_index
    !     _split_surface_layer!                   -> split_layer (k = 1)
    !     _merge_surface_layer!                   -> merge_layer (k = 1)
    !     _merge_bottom_layer!                    -> merge_bottom_layer
    !     _remove_surface_layer!                  -> remove_surface_layer
    !     _remove_depleted_surface_and_route_water! -> remove_depleted_surface_and_route_water
    !     _continuous_bottom_deplete!             -> continuous_bottom_deplete
    !     _free_slot_for_surface_split!           -> free_slot_for_split (k = 1)
    !     _enforce_snow_depth_cap!                -> enforce_snow_depth_cap
    !     _cap_near_surface_layer_thicknesses!    -> cap_near_surface_layer_thicknesses
    !     _fill_near_surface_layer_thicknesses!   -> fill_near_surface_layer_thicknesses
    !     _remesh_near_surface_layers!            -> remesh_near_surface_layers
    ! plus chion's rebalance_layer, the split and merge loops of
    ! _apply_accumulation_resolved! (accumulation.jl) for one layer k.
    !
    ! THE MASS LAYER k: Chion.jl splits and merges the surface layer only.
    ! split_layer, merge_layer and free_slot_for_split act on layer k, the
    ! top MASS layer; k = 1 is Chion.jl's surface layer. With fine
    ! near-surface layers chion applies them also to the first layer below
    ! the fine ones (docs/porting_notes.md D32), see remesh_near_surface_layers.
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
    !   > 0                free_slot_for_split bottom-mass test,
    !                      depth-cap layer inclusion, mass-weighted-mean guard
    !   > TOL_TINY         depletion remainder loop, depth-cap density guard
    !   < TOL_EMPTY_LAYER  merge_layer single-layer collapse,
    !                      continuous_bottom_deplete empty-layer skip and trim
    !
    ! DEPTH CAP: enforce_snow_depth_cap caps the total solid depth at the
    ! constant BESSI_REFERENCE_SNOW_DEPTH_M = 22.5 m (Chion.jl 03bb445),
    ! INDEPENDENT of the configured Ntot and of mass_split.
    !
    ! FINE NEAR-SURFACE LAYERS (Chion.jl 03bb445): the remesh holds the top
    ! NEAR_SURFACE_LAYERS layers at maximum thicknesses h_max(k) [m] (0 = no
    ! limit, Julia's Inf; the limited layers are a leading block,
    ! bessi_par_validate). It conserves solid mass, liquid water, volume and
    ! sensible enthalpy (temperature by mass_weighted_mean, D31). chion then
    ! splits and merges the first layer below them by mass (D32; not under
    ! legacy_chion).

    use chion_defs, only : wp, wp_acc, io_unit_err, TOL_TINY, TOL_EMPTY_LAYER, &
                           BESSI_REFERENCE_SNOW_DEPTH_M, &
                           NEAR_SURFACE_SPLIT_MERGE_BELOW, &
                           chion_const_class

    implicit none

    private

    public :: reset_layer_at_index
    public :: split_layer
    public :: merge_layer
    public :: merge_bottom_layer
    public :: remove_surface_layer
    public :: remove_depleted_surface_and_route_water
    public :: continuous_bottom_deplete
    public :: free_slot_for_split
    public :: rebalance_layer
    public :: enforce_snow_depth_cap
    public :: near_surface_layer_count
    public :: cap_near_surface_layer_thicknesses
    public :: fill_near_surface_layer_thicknesses
    public :: remesh_near_surface_layers

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

    pure function safe_positive(x) result(y)
        ! Chion.jl _safe_positive (src/processes/energy_flux.jl:5): floor at
        ! EPS_TINY, used to protect divisions. wp_acc, like the remesh locals
        ! it guards (snow_vapor has the wp version).

        implicit none

        real(wp_acc), intent(IN) :: x
        real(wp_acc) :: y

        if (x .gt. TOL_TINY) then
            y = x
        else
            y = TOL_TINY
        end if

        return

    end function safe_positive

    pure function mass_weighted_mean(m1,x1,m2,x2) result(xbar)
        ! Chion.jl _mass_weighted_mean: mass-weighted mean of two layer
        ! properties, returning zero when the combined mass is non-positive.
        !
        ! Written as x1 + w2*(x2 - x1), w2 = m2/(m1 + m2), not Chion.jl's
        ! (m1*x1 + m2*x2)/(m1 + m2) (docs/porting_notes.md D31): equal values
        ! mix to exactly that value (two layers at T0 stay at T0, not T0 - 1
        ! ulp), and only the difference is weighted, so it is better
        ! conditioned. Evaluated in wp_acc, since the two masses can differ by
        ! orders of magnitude just after a split.

        implicit none

        real(wp), intent(IN) :: m1, x1
        real(wp), intent(IN) :: m2, x2
        real(wp) :: xbar

        ! Local variables
        real(wp_acc) :: total_mass, w2

        total_mass = real(m1,wp_acc) + real(m2,wp_acc)

        if (total_mass .gt. 0.0_wp_acc) then
            w2   = real(m2,wp_acc)/total_mass
            xbar = real(real(x1,wp_acc) + w2*(real(x2,wp_acc) - real(x1,wp_acc)), wp)
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
    ! split_layer
    ! Chion.jl layer_structure.jl:51-91 (_split_surface_layer!, k = 1)
    ! =====================================================================

    subroutine split_layer(mass,mass_w,density,temperature,n,k,Ntot,mass_max,mass_split)
        ! Split mass layer k in two when it exceeds mass_max, provided a free
        ! slot exists (n < Ntot). The new layer k+1 receives exactly
        ! mass_split and layer k keeps the remainder; density and temperature
        ! are copied unchanged into both halves, and liquid water is
        ! partitioned by the solid-mass fraction. Layers below k shift down.
        !
        ! No-ops unless BOTH conditions hold. The caller is responsible for
        ! freeing a slot first (free_slot_for_split) when n == Ntot.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n                ! active layer count
        integer,  intent(IN)    :: k                ! mass layer to split, 1 <= k <= Ntot
        integer,  intent(IN)    :: Ntot             ! maximum layer count
        real(wp), intent(IN)    :: mass_max         ! [kg m-2] split trigger
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] target layer mass

        ! Local variables
        integer  :: j, new_n
        real(wp) :: layer_mass, layer_mass_w
        real(wp) :: layer_density, layer_temperature
        real(wp) :: water_fraction

        layer_mass = mass(k)

        ! Julia: if !(n < Ntot && surface_mass > mass_max) return
        ! Both operands are safe to evaluate, so a single .and. is fine here.
        ! An inactive slot k > n holds no mass, so it never splits.
        if (.not. (n .lt. Ntot .and. layer_mass .gt. mass_max)) return

        layer_mass_w      = mass_w(k)
        layer_density     = density(k)
        layer_temperature = temperature(k)

        new_n = n + 1
        n     = new_n

        ! Shift layers k+1..new_n-1 down by one, leaving slot k+1 free.
        do j = new_n, k+2, -1
            mass(j)        = mass(j-1)
            mass_w(j)      = mass_w(j-1)
            density(j)     = density(j-1)
            temperature(j) = temperature(j-1)
        end do

        mass(k+1) = mass_split
        mass(k)   = layer_mass - mass_split

        ! layer_mass > mass_max >= 0, so this division is safe.
        water_fraction = mass_split / layer_mass
        mass_w(k+1)    = layer_mass_w * water_fraction
        mass_w(k)      = layer_mass_w * (1.0_wp - water_fraction)

        density(k)       = layer_density
        density(k+1)     = layer_density
        temperature(k)   = layer_temperature
        temperature(k+1) = layer_temperature

        return

    end subroutine split_layer

    ! =====================================================================
    ! merge_layer
    ! Chion.jl layer_structure.jl:99-172 (_merge_surface_layer!, k = 1)
    ! =====================================================================

    subroutine merge_layer(mass,mass_w,density,temperature,n,k,mass_split,mass_min,c)
        ! Rebalance or merge mass layer k with layer k+1 when layer k falls
        ! below mass_min. Three outcomes:
        !   1. n == k and layer k is essentially empty (< TOL_EMPTY_LAYER): it
        !      is dropped (k = 1: the column collapses to zero layers).
        !   2. combined mass > 2*mass_split: PARTIAL TRANSFER. Layer k is
        !      topped back up to exactly mass_split from layer k+1; the layer
        !      count is unchanged.
        !   3. otherwise: FULL MERGE. Layers k and k+1 are combined into layer
        !      k and everything below shifts up by one.
        !
        ! PORTING NOTE: the Julia signature also carries Ntot, which its body
        ! never uses. Dropped here rather than carried as an unused dummy.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        integer,  intent(IN)    :: k                ! mass layer to merge, k >= 1
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] target layer mass
        real(wp), intent(IN)    :: mass_min         ! [kg m-2] merge trigger
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer  :: j, new_n
        real(wp) :: layer_mass, below_mass
        real(wp) :: transferred_to_layer, transferred_water
        real(wp) :: dens_1, dens_2, temp_1, temp_2
        real(wp_acc) :: combined_mass

        if (n .ge. k) then
            layer_mass = mass(k)
        else
            layer_mass = 0.0_wp
        end if

        ! Threshold here is TOL_EMPTY_LAYER (strict <), not TOL_TINY.
        if (n .eq. k .and. real(layer_mass,wp_acc) .lt. TOL_EMPTY_LAYER) then
            n = k - 1
            call reset_layer_at_index(mass,mass_w,density,temperature,k,c)
            return
        else if (n .le. k .or. layer_mass .ge. mass_min) then
            return
        end if

        below_mass    = mass(k+1)
        combined_mass = real(layer_mass,wp_acc) + real(below_mass,wp_acc)

        dens_1 = density(k)
        dens_2 = density(k+1)
        temp_1 = temperature(k)
        temp_2 = temperature(k+1)

        if (combined_mass .gt. 2.0_wp_acc*real(mass_split,wp_acc)) then

            ! --- Partial transfer: top layer k back up to mass_split.
            ! below_mass = combined - layer > 2*mass_split - mass_min, which
            ! is strictly positive for any sane parameter set, so the
            ! division below is safe. Julia does not guard it either.
            transferred_to_layer = mass_split - layer_mass
            transferred_water    = transferred_to_layer / below_mass * mass_w(k+1)

            mass(k)     = mass_split
            mass(k+1)   = real(combined_mass - real(mass_split,wp_acc),wp)
            mass_w(k)   = mass_w(k)   + transferred_water
            mass_w(k+1) = mass_w(k+1) - transferred_water

            density(k)     = mass_weighted_mean(layer_mass,dens_1, &
                                                transferred_to_layer,dens_2)
            temperature(k) = mass_weighted_mean(layer_mass,temp_1, &
                                                transferred_to_layer,temp_2)
            return

        end if

        ! --- Full merge of layers k and k+1.
        mass(k)   = real(combined_mass,wp)
        mass_w(k) = mass_w(k) + mass_w(k+1)

        density(k)     = mass_weighted_mean(layer_mass,dens_1,below_mass,dens_2)
        temperature(k) = mass_weighted_mean(layer_mass,temp_1,below_mass,temp_2)

        new_n = n - 1
        n     = new_n

        do j = k+1, new_n
            mass(j)        = mass(j+1)
            mass_w(j)      = mass_w(j+1)
            density(j)     = density(j+1)
            temperature(j) = temperature(j+1)
        end do

        call reset_layer_at_index(mass,mass_w,density,temperature,new_n+1,c)

        return

    end subroutine merge_layer

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
    ! free_slot_for_split
    ! Chion.jl layer_structure.jl:374-434 (_free_slot_for_surface_split!, k = 1)
    ! =====================================================================

    subroutine free_slot_for_split(mass,mass_w,density,temperature,n,k, &
                                   mass_base,smb_ice,runoff,t_srf,albedo, &
                                   Ntot,mass_max,c)
        ! Make room for a split of mass layer k when the column is already at
        ! Ntot.
        !
        ! Ntot == k SPECIAL CASE: layer k is the last slot, so there is none
        ! to free; instead its own overflow beyond mass_max is exported
        ! basally. For k = 1 (Ntot == 1) this is the only path by which a
        ! single-layer column sheds mass at the base during accumulation.
        !
        ! Otherwise the deepest layer is depleted in full (routing its water to
        ! runoff), or simply reset if it carries no mass.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        integer,  intent(IN)    :: k                ! mass layer to be split
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

        if (Ntot .eq. k) then

            overflow = max(real(mass(k),wp_acc) - real(mass_max,wp_acc), 0.0_wp_acc)

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
            write(io_unit_err,*) "free_slot_for_split:: Error: &
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

    end subroutine free_slot_for_split

    ! =====================================================================
    ! rebalance_layer
    ! Chion.jl accumulation.jl:110-131 (the split and merge loops of
    ! _apply_accumulation_resolved!, k = 1)
    ! =====================================================================

    subroutine rebalance_layer(mass,mass_w,density,temperature,n,k, &
                               mass_base,smb_ice,runoff,t_srf,albedo, &
                               Ntot,mass_max,mass_split,mass_min,c)
        ! Keep mass layer k within [mass_min, mass_max]:
        !   split loop  while mass(k) > mass_max: split off mass_split below
        !               it, first freeing a slot at Ntot -- by depleting the
        !               bottom layer when at most two mass layers fit below
        !               the k-1 above (k = 1: Ntot <= 2), else by merging the
        !               two deepest layers;
        !   merge loop  while n > k and mass(k) < mass_min: merge or top up
        !               from layer k+1.
        ! k = 1 is Chion.jl's surface rebalance in accumulation.

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        integer,  intent(IN)    :: k                ! mass layer, 1 <= k <= Ntot
        real(wp_acc), intent(INOUT) :: mass_base
        real(wp_acc), intent(INOUT) :: smb_ice
        real(wp_acc), intent(INOUT) :: runoff
        real(wp), intent(INOUT) :: t_srf
        real(wp), intent(INOUT) :: albedo
        integer,  intent(IN)    :: Ntot
        real(wp), intent(IN)    :: mass_max         ! [kg m-2] split trigger
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] mass left below
        real(wp), intent(IN)    :: mass_min         ! [kg m-2] merge trigger
        type(chion_const_class), intent(IN) :: c

        ! --- Split loop. mass(k) is in bounds for any n (k <= Ntot), so the
        !     single .and. is safe.
        do while (n .ge. k .and. mass(k) .gt. mass_max)

            if (n .eq. Ntot) then

                if (Ntot - k .le. 1) then
                    call free_slot_for_split(mass,mass_w,density,temperature,n,k, &
                                             mass_base,smb_ice,runoff,t_srf,albedo, &
                                             Ntot,mass_max,c)
                else
                    call merge_bottom_layer(mass,mass_w,density,temperature,n, &
                                            mass_base,smb_ice,c)
                end if

                ! Re-test after freeing a slot: either layer k is gone or it
                ! is already back under mass_max.
                if (n .lt. k) exit
                if (mass(k) .le. mass_max) exit

            end if

            call split_layer(mass,mass_w,density,temperature,n,k,Ntot,mass_max,mass_split)

        end do

        ! --- Merge loop.
        do while (n .gt. k .and. mass(k) .lt. mass_min)
            call merge_layer(mass,mass_w,density,temperature,n,k,mass_split,mass_min,c)
        end do

        return

    end subroutine rebalance_layer

    ! =====================================================================
    ! enforce_snow_depth_cap
    ! Chion.jl layer_structure.jl:442-513
    ! =====================================================================

    subroutine enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                      mass_base,smb_ice,runoff,t_srf,albedo,c)
        ! Cap the total active snow depth after accumulation, exporting the
        ! excess from the base of the column.
        !
        ! THE CAP IGNORES THE CONFIGURED Ntot AND mass_split: reference_depth
        ! is the constant BESSI_REFERENCE_SNOW_DEPTH_M = 22.5 m (Chion.jl
        ! 03bb445). It replaced 15*mass_split*1.5/300 (trap 11 in
        ! docs/PLAN.md section 5), which equals 22.5 m at the default
        ! mass_split = 300.
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

        reference_depth = real(BESSI_REFERENCE_SNOW_DEPTH_M,wp_acc)

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

    ! =====================================================================
    ! Fine near-surface layers
    ! Chion.jl layer_structure.jl (03bb445) _cap_near_surface_layer_thicknesses!,
    ! _fill_near_surface_layer_thicknesses!, _remesh_near_surface_layers!
    ! =====================================================================

    pure function near_surface_layer_count(h_max) result(n_fine)
        ! Number of thickness-limited near-surface layers: the leading entries
        ! of h_max that are > 0 (0 = no limit). bessi_par_validate requires the
        ! limited layers to be a leading block, so this is also the index of
        ! the deepest one; 0 = no fine layers.

        implicit none

        real(wp), intent(IN) :: h_max(:)     ! (NEAR_SURFACE_LAYERS) [m]
        integer :: n_fine

        ! Local variables
        integer :: k

        n_fine = 0

        do k = 1, size(h_max)
            if (h_max(k) .le. 0.0_wp) exit
            n_fine = k
        end do

        return

    end function near_surface_layer_count

    subroutine cap_near_surface_layer_thicknesses(mass,mass_w,density,temperature,n, &
                                                  Ntot,h_max,c)
        ! The downward half of the remesh. For each limited layer k (top
        ! down, k <= n), mass above rho_k*h_max(k) moves into layer k+1 with
        ! its share of liquid water. Layer k+1 takes the combined mass, the
        ! volume-conserving density m/(m_excess/rho_k + m_below/rho_below)
        ! and the mass-weighted temperature. Density and temperature of
        ! layer k are unchanged.
        !
        ! If k is the deepest active layer, a new layer n+1 is opened for the
        ! excess -- unless n == Ntot, in which case layer k keeps it and the
        ! pass ends (Julia: "a full column keeps its deepest near-surface
        ! layer as it is").
        !
        ! Precision: the excess is taken against the capped mass AS STORED,
        ! so the pair (k, k+1) is conserved up to the one rounding of
        ! mass(k+1) (identical to Julia when wp = dp).

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        integer,  intent(IN)    :: Ntot
        real(wp), intent(IN)    :: h_max(:)          ! (NEAR_SURFACE_LAYERS) [m], 0 = no limit
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer      :: k
        real(wp)     :: mass_cap, excess_wp
        real(wp_acc) :: layer_mass, max_mass, excess_mass, excess_water
        real(wp_acc) :: below_mass, combined_mass, combined_volume

        do k = 1, size(h_max)

            if (k .gt. n) exit
            if (h_max(k) .le. 0.0_wp) cycle

            layer_mass = real(mass(k),wp_acc)
            max_mass   = real(density(k),wp_acc)*real(h_max(k),wp_acc)

            if (layer_mass .le. max_mass) cycle

            if (k .eq. n) then
                if (k .eq. Ntot) exit
                n = n + 1
                call reset_layer_at_index(mass,mass_w,density,temperature,n,c)
            end if

            mass_cap     = real(max_mass,wp)
            excess_mass  = layer_mass - real(mass_cap,wp_acc)
            excess_water = real(mass_w(k),wp_acc)*excess_mass/safe_positive(layer_mass)

            below_mass      = real(mass(k+1),wp_acc)
            combined_mass   = excess_mass + below_mass
            combined_volume = excess_mass/safe_positive(real(density(k),wp_acc)) &
                            + below_mass/safe_positive(real(density(k+1),wp_acc))

            excess_wp = real(excess_mass,wp)
            temperature(k+1) = mass_weighted_mean(excess_wp,temperature(k), &
                                                  mass(k+1),temperature(k+1))
            density(k+1)     = real(combined_mass/safe_positive(combined_volume),wp)

            mass(k)     = mass_cap
            mass_w(k)   = real(real(mass_w(k),wp_acc)   - excess_water,wp)
            mass(k+1)   = real(combined_mass,wp)
            mass_w(k+1) = real(real(mass_w(k+1),wp_acc) + excess_water,wp)

        end do

        return

    end subroutine cap_near_surface_layer_thicknesses

    subroutine fill_near_surface_layer_thicknesses(mass,mass_w,density,temperature,n, &
                                                   h_target,c)
        ! The upward half of the remesh. Each limited layer k < n thinner
        ! than h_target(k) pulls min(m_donor, deficit*rho_donor) from layer
        ! k+1, with liquid water in proportion; layer k takes the
        ! volume-conserving density and the mass-weighted temperature. A
        ! donor left with <= TOL_EMPTY_LAYER is removed (its residue with
        ! it, as in Julia) and the next layer becomes the donor. A shallow
        ! column is not padded: its deepest layer may stay thin.
        !
        ! Julia loops while the recomputed deficit exceeds EPS_TINY. After a
        ! transfer that leaves the donor non-empty the receiver is full by
        ! construction: in dp the recomputed deficit is ~1e-18 m and Julia
        ! exits, but in sp it is a few ulp of h (~1e-9 m) and the loop would
        ! keep moving round-off. So the loop ends there instead; it continues
        ! only after a donor is exhausted (identical to Julia when wp = dp).

        implicit none

        real(wp), intent(INOUT) :: mass(:)
        real(wp), intent(INOUT) :: mass_w(:)
        real(wp), intent(INOUT) :: density(:)
        real(wp), intent(INOUT) :: temperature(:)
        integer,  intent(INOUT) :: n
        real(wp), intent(IN)    :: h_target(:)       ! (NEAR_SURFACE_LAYERS) [m], 0 = no limit
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer      :: k, kd, j
        real(wp_acc) :: receiver_mass, receiver_volume, missing_volume
        real(wp_acc) :: donor_mass, donor_density, donor_water
        real(wp_acc) :: transferred_mass, transferred_water
        real(wp_acc) :: combined_mass, combined_volume, remaining_donor_mass

        do k = 1, size(h_target)

            if (k .ge. n) exit
            if (h_target(k) .le. 0.0_wp) cycle

            do while (k .lt. n)

                receiver_mass   = real(mass(k),wp_acc)
                receiver_volume = receiver_mass/safe_positive(real(density(k),wp_acc))
                missing_volume  = real(h_target(k),wp_acc) - receiver_volume

                if (missing_volume .le. TOL_TINY) exit

                kd = k + 1

                donor_mass       = real(mass(kd),wp_acc)
                donor_density    = real(density(kd),wp_acc)
                transferred_mass = min(donor_mass,missing_volume*donor_density)

                if (transferred_mass .le. TOL_TINY) exit

                donor_water       = real(mass_w(kd),wp_acc)
                transferred_water = donor_water*transferred_mass/safe_positive(donor_mass)

                combined_mass   = receiver_mass + transferred_mass
                combined_volume = receiver_volume + transferred_mass/safe_positive(donor_density)

                temperature(k) = mass_weighted_mean(mass(k),temperature(k), &
                                                    real(transferred_mass,wp),temperature(kd))
                mass(k)        = real(combined_mass,wp)
                mass_w(k)      = real(real(mass_w(k),wp_acc) + transferred_water,wp)
                density(k)     = real(combined_mass/safe_positive(combined_volume),wp)

                remaining_donor_mass = donor_mass - transferred_mass
                mass(kd)   = real(remaining_donor_mass,wp)
                mass_w(kd) = real(donor_water - transferred_water,wp)

                ! Receiver full, donor left: done with layer k (see header).
                if (remaining_donor_mass .gt. TOL_EMPTY_LAYER) exit

                ! Donor exhausted: shift the layers below it up by one.
                do j = kd, n-1
                    mass(j)        = mass(j+1)
                    mass_w(j)      = mass_w(j+1)
                    density(j)     = density(j+1)
                    temperature(j) = temperature(j+1)
                end do
                call reset_layer_at_index(mass,mass_w,density,temperature,n,c)
                n = n - 1

            end do

        end do

        return

    end subroutine fill_near_surface_layer_thicknesses

    subroutine remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                          mass_base,smb_ice,runoff,t_srf,albedo, &
                                          Ntot,mass_max,mass_split,mass_min,h_max,c)
        ! Conservative remesh of the fine near-surface layers: cap downward,
        ! then fill upward. Called twice per (sub)step by the BESSI kernel,
        ! after accumulation and after refreezing. A no-op without limits.
        !
        ! DEVIATION (docs/porting_notes.md D32; not under legacy_chion): the
        ! first layer below the fine ones, k0, is then kept within
        ! [mass_min, mass_max] by the mass-based split and merge that Chion.jl
        ! applies to the surface layer only (rebalance_layer at k0). In
        ! Chion.jl everything the cap pushes below the fine layers stays in
        ! layer k0, which is never split or merged, so a firn column ends as
        ! the fine layers over one cell of up to the 22.5 m depth cap. The
        ! accumulators are touched only when the split frees a slot at Ntot
        ! (bottom merge or depletion, as in accumulation).

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
        real(wp), intent(IN)    :: mass_max         ! [kg m-2] split trigger
        real(wp), intent(IN)    :: mass_split       ! [kg m-2] mass left below
        real(wp), intent(IN)    :: mass_min         ! [kg m-2] merge trigger
        real(wp), intent(IN)    :: h_max(:)         ! (NEAR_SURFACE_LAYERS) [m], 0 = no limit
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer :: k0

        call cap_near_surface_layer_thicknesses(mass,mass_w,density,temperature,n,Ntot,h_max,c)
        call fill_near_surface_layer_thicknesses(mass,mass_w,density,temperature,n,h_max,c)

        if (.not. NEAR_SURFACE_SPLIT_MERGE_BELOW) return

        k0 = near_surface_layer_count(h_max) + 1

        if (k0 .gt. 1 .and. k0 .le. Ntot) then
            call rebalance_layer(mass,mass_w,density,temperature,n,k0, &
                                 mass_base,smb_ice,runoff,t_srf,albedo, &
                                 Ntot,mass_max,mass_split,mass_min,c)
        end if

        return

    end subroutine remesh_near_surface_layers

end module snow_layers
