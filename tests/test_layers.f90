program test_layers
    ! WP4 acceptance test: snow_layers.
    !
    ! Drives accumulation/melt-like sequences through split, surface merge,
    ! bottom merge, basal depletion and the depth cap, and asserts after every
    ! operation that
    !
    !     sum(mass) + sum(mass_w) + mass_base + runoff
    !
    ! is conserved to 1e-6 relative. The layered state is sp, so 1e-6 is the
    ! right target (docs/PLAN.md section 3.1); the check itself accumulates in
    ! wp_acc because mass_base and runoff are wp_acc.
    !
    ! Covered explicitly: Ntot == 1 and Ntot == 2, a column driven to zero
    ! layers and back, and the combined_density > rho_i basal-ice export path
    ! in merge_bottom_layer.
    !
    ! Fine near-surface layers (Chion.jl 03bb445, plan C4): the remesh
    ! (cap down, fill up) conserves solid mass, liquid water, volume and
    ! sensible enthalpy sum(m*T); Chion.jl's test_geometric_accumulation
    ! (zero-input accumulation + remesh is the identity) is ported. Below the
    ! fine layers chion splits and merges the first layer by mass (C4b, D32);
    ! legacy_chion builds check Chion.jl's behaviour there instead.

    use chion_defs,  only : wp, wp_acc, chion_const_class, chion_const_init, &
                            BESSI_REFERENCE_SNOW_DEPTH_M, NEAR_SURFACE_LAYERS, &
                            NEAR_SURFACE_SPLIT_MERGE_BELOW
    use snow_layers
    use snow_accumulation, only : apply_accumulation

    implicit none

    integer, parameter :: NMAX = 20

    type(chion_const_class) :: c
    integer  :: nfail

    ! Column state, sized generously; the "Ntot" in use is passed per test.
    !
    ! Main-program variables are implicitly SAVE (F2008 5.3.16). The attribute
    ! is stated on the arrays because ifx 2023.2 with -qopenmp otherwise hands
    ! the internal procedures below (column_total, clear_column, ...) a
    ! separate, all-zero copy of them, at any -O level, so every conservation
    ! check compared 0 with 0 or failed. A bare SAVE statement does not avoid
    ! it; gfortran is unaffected.
    real(wp), save :: mass(NMAX), mass_w(NMAX), density(NMAX), temperature(NMAX)
    integer  :: n

    real(wp_acc) :: mass_base, smb_ice, runoff
    real(wp)     :: t_srf, albedo

    real(wp_acc) :: total_ref, total_now
    real(wp_acc) :: ice_out, run_out

    integer  :: k, istep
    real(wp) :: mass_max, mass_split, mass_min

    ! Fine near-surface layers
    real(wp)     :: h_fine(NEAR_SURFACE_LAYERS), h_off(NEAR_SURFACE_LAYERS)
    real(wp)     :: mass_before(NMAX), temp_before(NMAX), snow_age
    real(wp_acc) :: vol_ref, ent_ref, water_ref
    logical      :: at_target

    nfail = 0

    call chion_const_init(c)

    mass_max   = 500.0_wp
    mass_split = 300.0_wp
    mass_min   = 100.0_wp

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP4 acceptance test: snow_layers"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! =====================================================================
    ! reset_layer_at_index
    ! =====================================================================
    write(*,"(a)") "--- reset_layer_at_index ---"

    call clear_column()
    mass(3) = 111.0_wp ; mass_w(3) = 22.0_wp
    density(3) = 400.0_wp ; temperature(3) = 260.0_wp
    n = 3

    call reset_layer_at_index(mass,mass_w,density,temperature,3,c)

    call check("mass zeroed",        mass(3)    .eq. 0.0_wp, nfail)
    call check("mass_w zeroed",      mass_w(3)  .eq. 0.0_wp, nfail)
    call check("density zeroed",     density(3) .eq. 0.0_wp, nfail)
    call check_val("temperature set to T0", temperature(3), c%T0, nfail)

    ! =====================================================================
    ! split_layer
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- split_layer ---"

    call clear_column()
    n = 2
    mass(1) = 700.0_wp ; mass_w(1) = 70.0_wp
    density(1) = 350.0_wp ; temperature(1) = 265.0_wp
    mass(2) = 250.0_wp ; mass_w(2) = 5.0_wp
    density(2) = 500.0_wp ; temperature(2) = 268.0_wp

    total_ref = column_total()

    call split_layer(mass,mass_w,density,temperature,n,1,5,mass_max,mass_split)

    call check("n incremented to 3", n .eq. 3, nfail)
    call check_val("layer 2 holds mass_split",  mass(2), mass_split, nfail)
    call check_val("layer 1 holds remainder",   mass(1), 700.0_wp-mass_split, nfail)
    call check_val("old layer 2 pushed to 3",   mass(3), 250.0_wp, nfail)
    call check_val("old layer 2 water pushed",  mass_w(3), 5.0_wp, nfail)
    call check_val("water split by mass fraction", mass_w(2), &
                   70.0_wp*(mass_split/700.0_wp), nfail)
    call check_val("density copied to layer 2", density(2), 350.0_wp, nfail)
    call check_val("temperature copied to layer 2", temperature(2), 265.0_wp, nfail)
    call check_conserve("split conserves mass", total_ref, nfail)

    ! No free slot -> no-op.
    call clear_column()
    n = 3
    mass(1:3)   = [700.0_wp, 300.0_wp, 300.0_wp]
    density(1:3) = 350.0_wp
    total_ref = column_total()
    call split_layer(mass,mass_w,density,temperature,n,1,3,mass_max,mass_split)
    call check("no split when n == Ntot", n .eq. 3 .and. mass(1) .eq. 700.0_wp, nfail)

    ! Below mass_max -> no-op.
    call clear_column()
    n = 1
    mass(1) = 400.0_wp ; density(1) = 350.0_wp
    call split_layer(mass,mass_w,density,temperature,n,1,5,mass_max,mass_split)
    call check("no split below mass_max", n .eq. 1 .and. mass(1) .eq. 400.0_wp, nfail)

    ! =====================================================================
    ! merge_layer: partial-transfer branch
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- merge_layer (partial transfer) ---"

    call clear_column()
    n = 3
    mass(1) = 50.0_wp  ; mass_w(1) = 5.0_wp
    density(1) = 300.0_wp ; temperature(1) = 260.0_wp
    mass(2) = 700.0_wp ; mass_w(2) = 70.0_wp
    density(2) = 500.0_wp ; temperature(2) = 270.0_wp
    mass(3) = 400.0_wp ; mass_w(3) = 0.0_wp
    density(3) = 600.0_wp ; temperature(3) = 272.0_wp

    total_ref = column_total()

    ! combined = 750 > 2*mass_split = 600 -> partial transfer, n unchanged.
    call merge_layer(mass,mass_w,density,temperature,n,1,mass_split,mass_min,c)

    call check("n unchanged on partial transfer", n .eq. 3, nfail)
    call check_val("surface topped up to mass_split", mass(1), mass_split, nfail)
    call check_val("layer 2 keeps the remainder", mass(2), 750.0_wp-mass_split, nfail)
    call check_val("water transferred in proportion", mass_w(1), &
                   5.0_wp + (mass_split-50.0_wp)/700.0_wp*70.0_wp, nfail)
    call check_val("surface density is mass-weighted", density(1), &
                   (50.0_wp*300.0_wp + (mass_split-50.0_wp)*500.0_wp) &
                   /(50.0_wp + (mass_split-50.0_wp)), nfail)
    call check_conserve("partial transfer conserves mass", total_ref, nfail)

    ! =====================================================================
    ! merge_layer: full-merge branch
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- merge_layer (full merge) ---"

    call clear_column()
    n = 3
    mass(1) = 50.0_wp  ; mass_w(1) = 5.0_wp
    density(1) = 300.0_wp ; temperature(1) = 260.0_wp
    mass(2) = 200.0_wp ; mass_w(2) = 20.0_wp
    density(2) = 500.0_wp ; temperature(2) = 270.0_wp
    mass(3) = 400.0_wp ; mass_w(3) = 8.0_wp
    density(3) = 600.0_wp ; temperature(3) = 272.0_wp

    total_ref = column_total()

    ! combined = 250 <= 600 -> full merge, n drops to 2.
    call merge_layer(mass,mass_w,density,temperature,n,1,mass_split,mass_min,c)

    call check("n decremented on full merge", n .eq. 2, nfail)
    call check_val("merged surface mass", mass(1), 250.0_wp, nfail)
    call check_val("merged surface water", mass_w(1), 25.0_wp, nfail)
    call check_val("merged density is mass-weighted", density(1), &
                   (50.0_wp*300.0_wp+200.0_wp*500.0_wp)/250.0_wp, nfail)
    call check_val("merged temperature is mass-weighted", temperature(1), &
                   (50.0_wp*260.0_wp+200.0_wp*270.0_wp)/250.0_wp, nfail)
    call check_val("layer 3 shifted up to layer 2", mass(2), 400.0_wp, nfail)
    call check("vacated layer 3 is reset", mass(3) .eq. 0.0_wp &
               .and. density(3) .eq. 0.0_wp, nfail)
    call check_val("vacated layer 3 temperature is T0", temperature(3), c%T0, nfail)
    call check_conserve("full merge conserves mass", total_ref, nfail)

    ! Single nearly-empty layer -> column collapses.
    call clear_column()
    n = 1
    mass(1) = 1.0e-11_wp     ! below TOL_EMPTY_LAYER = 1e-10
    density(1) = 300.0_wp
    call merge_layer(mass,mass_w,density,temperature,n,1,mass_split,mass_min,c)
    call check("n=1 below TOL_EMPTY_LAYER collapses to n=0", n .eq. 0, nfail)

    call clear_column()
    n = 1
    mass(1) = 1.0e-9_wp      ! above TOL_EMPTY_LAYER
    density(1) = 300.0_wp
    call merge_layer(mass,mass_w,density,temperature,n,1,mass_split,mass_min,c)
    call check("n=1 above TOL_EMPTY_LAYER survives", n .eq. 1, nfail)

    ! Surface already at/above mass_min -> no-op.
    call clear_column()
    n = 2
    mass(1) = 150.0_wp ; mass(2) = 150.0_wp
    density(1:2) = 300.0_wp
    call merge_layer(mass,mass_w,density,temperature,n,1,mass_split,mass_min,c)
    call check("no merge when surface >= mass_min", n .eq. 2, nfail)

    ! Two layers at T0 merge to exactly T0 (D31). With these masses the
    ! sum-then-divide form gives T0 - 1 ulp at dp, which flips the aging
    ! albedo's T >= T0 timescale test.
    call clear_column()
    n = 2
    mass(1) = 94.2237843373378_wp ; mass(2) = 305.6260503473118_wp
    density(1) = 350.0_wp ; density(2) = 350.0_wp
    temperature(1:2) = c%T0
    call merge_layer(mass,mass_w,density,temperature,n,1,mass_split,mass_min,c)
    call check("full merge of two layers at T0 is exactly T0", &
               n .eq. 1 .and. temperature(1) .eq. c%T0, nfail)
    call check("full merge of two equal densities is exact", density(1) .eq. 350.0_wp, nfail)

    ! =====================================================================
    ! merge_bottom_layer: ordinary branch
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- merge_bottom_layer (below ice density) ---"

    call clear_column()
    call clear_accum()
    n = 3
    mass(1) = 300.0_wp ; density(1) = 350.0_wp ; temperature(1) = 260.0_wp
    mass(2) = 300.0_wp ; density(2) = 500.0_wp ; temperature(2) = 265.0_wp
    mass(3) = 100.0_wp ; density(3) = 700.0_wp ; temperature(3) = 270.0_wp
    mass_w(2) = 10.0_wp ; mass_w(3) = 4.0_wp

    total_ref = column_total()

    call merge_bottom_layer(mass,mass_w,density,temperature,n,mass_base,smb_ice,c)

    call check("n decremented", n .eq. 2, nfail)
    call check_val("combined bottom mass",  mass(2), 400.0_wp, nfail)
    call check_val("combined bottom water", mass_w(2), 14.0_wp, nfail)
    call check_val("combined bottom density", density(2), &
                   (300.0_wp*500.0_wp+100.0_wp*700.0_wp)/400.0_wp, nfail)
    call check_acc("no basal export below rho_i", mass_base, 0.0_wp_acc, nfail)
    call check_conserve("bottom merge conserves mass", total_ref, nfail)

    ! =====================================================================
    ! merge_bottom_layer: combined_density > rho_i export path
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- merge_bottom_layer (combined_density > rho_i export) ---"

    call clear_column()
    call clear_accum()
    n = 2
    mass(1) = 200.0_wp ; density(1) = 950.0_wp ; temperature(1) = 270.0_wp
    mass(2) = 200.0_wp ; density(2) = 990.0_wp ; temperature(2) = 272.0_wp
    mass_w(1) = 3.0_wp ; mass_w(2) = 1.0_wp

    total_ref = column_total()

    call merge_bottom_layer(mass,mass_w,density,temperature,n,mass_base,smb_ice,c)

    call check("n decremented to 1", n .eq. 1, nfail)
    call check_val("density pinned at rho_i", density(1), c%rho_i, nfail)
    call check_val("mass limited to what rho_i can hold", mass(1), &
                   400.0_wp*(c%rho_i/970.0_wp), nfail)
    call check_acc("excess exported to mass_base", mass_base, &
                   400.0_wp_acc*(1.0_wp_acc - real(c%rho_i,wp_acc)/970.0_wp_acc), nfail)
    call check_acc("smb_ice credited identically", smb_ice, mass_base, nfail)
    call check_val("liquid water fully retained", mass_w(1), 4.0_wp, nfail)
    call check_conserve("rho_i export conserves mass (via mass_base)", total_ref, nfail)

    ! =====================================================================
    ! remove_surface_layer / remove_depleted_surface_and_route_water
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- remove_surface_layer / route water ---"

    call clear_column()
    call clear_accum()
    n = 3
    mass(1:3)   = [10.0_wp, 200.0_wp, 300.0_wp]
    density(1:3) = [300.0_wp, 400.0_wp, 500.0_wp]
    temperature(1:3) = [260.0_wp, 265.0_wp, 270.0_wp]

    call remove_surface_layer(mass,mass_w,density,temperature,n,c)
    call check("n decremented", n .eq. 2, nfail)
    call check_val("layers shifted up", mass(1), 200.0_wp, nfail)
    call check_val("second layer shifted up", mass(2), 300.0_wp, nfail)
    call check("vacated layer 3 reset", mass(3) .eq. 0.0_wp, nfail)

    ! Water routed into layer 2 when one exists.
    call clear_column()
    call clear_accum()
    n = 2
    mass(1) = 0.0_wp   ; mass_w(1) = 7.0_wp
    mass(2) = 200.0_wp ; mass_w(2) = 3.0_wp ; density(2) = 400.0_wp

    total_ref = column_total()

    call remove_depleted_surface_and_route_water(mass,mass_w,density,temperature,n,runoff,c)
    call check("n decremented", n .eq. 1, nfail)
    call check_val("water pushed into the layer below", mass_w(1), 10.0_wp, nfail)
    call check_acc("no runoff generated", runoff, 0.0_wp_acc, nfail)
    call check_conserve("routing to layer below conserves mass", total_ref, nfail)

    ! Water goes to runoff when the surface layer is the last one.
    call clear_column()
    call clear_accum()
    n = 1
    mass(1) = 0.0_wp ; mass_w(1) = 7.0_wp

    total_ref = column_total()

    call remove_depleted_surface_and_route_water(mass,mass_w,density,temperature,n,runoff,c)
    call check("column emptied", n .eq. 0, nfail)
    call check_acc("water routed to runoff", runoff, 7.0_wp_acc, nfail)
    call check_conserve("routing to runoff conserves mass", total_ref, nfail)

    ! =====================================================================
    ! continuous_bottom_deplete
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- continuous_bottom_deplete ---"

    call clear_column()
    call clear_accum()
    n = 3
    mass(1:3)   = [300.0_wp, 300.0_wp, 300.0_wp]
    mass_w(1:3) = [ 10.0_wp,  20.0_wp,  30.0_wp]
    density(1:3) = 500.0_wp
    temperature(1:3) = [260.0_wp, 265.0_wp, 270.0_wp]

    total_ref = column_total()

    ! Remove 450: all of layer 3 (300) plus half of layer 2.
    call continuous_bottom_deplete(mass,mass_w,density,temperature,n, &
                                   mass_base,smb_ice,runoff,t_srf,albedo, &
                                   450.0_wp_acc,c,ice_out,run_out)

    call check("n decremented to 2", n .eq. 2, nfail)
    call check_val("partial bottom layer left", mass(2), 150.0_wp, nfail)
    call check_val("water removed in proportion", mass_w(2), 10.0_wp, nfail)
    call check_acc("ice_to_base = requested mass", ice_out, 450.0_wp_acc, nfail)
    call check_acc("mass_base = requested mass",   mass_base, 450.0_wp_acc, nfail)
    call check_acc("smb_ice = requested mass",     smb_ice, 450.0_wp_acc, nfail)
    call check_acc("runoff = 30 + 10",             runoff, 40.0_wp_acc, nfail)
    call check_val("t_srf tracks layer 1", t_srf, 260.0_wp, nfail)
    call check_conserve("depletion conserves mass", total_ref, nfail)

    ! Deplete more than the column holds: column empties, albedo -> alpha_ice.
    call clear_column()
    call clear_accum()
    n = 2
    mass(1:2)   = [100.0_wp, 100.0_wp]
    mass_w(1:2) = [  5.0_wp,   5.0_wp]
    density(1:2) = 500.0_wp
    albedo = 0.81_wp

    total_ref = column_total()

    call continuous_bottom_deplete(mass,mass_w,density,temperature,n, &
                                   mass_base,smb_ice,runoff,t_srf,albedo, &
                                   1000.0_wp_acc,c,ice_out,run_out)

    call check("column driven to zero layers", n .eq. 0, nfail)
    call check_acc("all solid mass exported", mass_base, 200.0_wp_acc, nfail)
    call check_acc("all liquid water to runoff", runoff, 10.0_wp_acc, nfail)
    call check_val("albedo forced to bare ice", albedo, c%alpha_ice, nfail)
    call check_val("t_srf reset to T0", t_srf, c%T0, nfail)
    call check_conserve("over-depletion conserves mass", total_ref, nfail)

    ! ...and back again: rebuild the column from zero layers.
    write(*,*)
    write(*,"(a)") "--- zero layers and back ---"

    n = 1
    mass(1) = 120.0_wp ; density(1) = 300.0_wp ; temperature(1) = 265.0_wp
    total_ref = column_total()
    call check("column rebuilt from n=0", n .eq. 1, nfail)

    do istep = 1, 6
        mass(1) = mass(1) + 200.0_wp
        total_ref = total_ref + 200.0_wp_acc
        do while (n .lt. 5 .and. mass(1) .gt. mass_max)
            call split_layer(mass,mass_w,density,temperature,n,1,5,mass_max,mass_split)
        end do
        call check_conserve("accumulate+split step conserves mass", total_ref, nfail)
    end do
    call check("column grew to several layers", n .ge. 3, nfail)

    ! Melt it away layer by layer, back to zero.
    do while (n .gt. 0)
        mass(1)   = 0.0_wp
        density(1) = 0.0_wp
        call remove_depleted_surface_and_route_water(mass,mass_w,density,temperature,n, &
                                                     runoff,c)
    end do
    call check("column returned to zero layers", n .eq. 0, nfail)

    ! =====================================================================
    ! free_slot_for_split: Ntot == 1 special case
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- free_slot_for_split (Ntot == 1) ---"

    call clear_column()
    call clear_accum()
    n = 1
    mass(1) = 800.0_wp ; mass_w(1) = 16.0_wp
    density(1) = 400.0_wp ; temperature(1) = 265.0_wp

    total_ref = column_total()

    call free_slot_for_split(mass,mass_w,density,temperature,n,1, &
                             mass_base,smb_ice,runoff,t_srf,albedo, &
                             1,mass_max,c)

    call check("still one layer", n .eq. 1, nfail)
    call check_val("surface trimmed back to mass_max", mass(1), mass_max, nfail)
    call check_acc("overflow exported basally", mass_base, 300.0_wp_acc, nfail)
    call check_acc("water removed in proportion", runoff, &
                   16.0_wp_acc*300.0_wp_acc/800.0_wp_acc, nfail)
    call check_conserve("Ntot=1 slot freeing conserves mass", total_ref, nfail)

    ! Below mass_max: nothing happens.
    call clear_accum()
    total_ref = column_total()
    call free_slot_for_split(mass,mass_w,density,temperature,n,1, &
                             mass_base,smb_ice,runoff,t_srf,albedo, &
                             1,mass_max,c)
    call check_acc("no export when below mass_max", mass_base, 0.0_wp_acc, nfail)
    call check_conserve("Ntot=1 no-op conserves mass", total_ref, nfail)

    ! =====================================================================
    ! free_slot_for_split: Ntot > 1
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- free_slot_for_split (Ntot == 2) ---"

    call clear_column()
    call clear_accum()
    n = 2
    mass(1) = 600.0_wp ; mass_w(1) = 6.0_wp ; density(1) = 350.0_wp
    mass(2) = 250.0_wp ; mass_w(2) = 5.0_wp ; density(2) = 500.0_wp

    total_ref = column_total()

    call free_slot_for_split(mass,mass_w,density,temperature,n,1, &
                             mass_base,smb_ice,runoff,t_srf,albedo, &
                             2,mass_max,c)

    call check("bottom layer consumed, n = 1", n .eq. 1, nfail)
    call check_acc("bottom mass exported basally", mass_base, 250.0_wp_acc, nfail)
    call check_acc("bottom water to runoff", runoff, 5.0_wp_acc, nfail)
    call check_conserve("Ntot=2 slot freeing conserves mass", total_ref, nfail)

    ! ...and the split can then proceed.
    call split_layer(mass,mass_w,density,temperature,n,1,2,mass_max,mass_split)
    call check("split now succeeds with the freed slot", n .eq. 2, nfail)
    call check_val("layer 2 holds mass_split", mass(2), mass_split, nfail)

    ! Massless bottom layer -> reset, no export.
    call clear_column()
    call clear_accum()
    n = 3
    mass(1:2)   = [600.0_wp, 200.0_wp]
    density(1:2) = 400.0_wp
    mass(3) = 0.0_wp ; density(3) = 0.0_wp

    total_ref = column_total()

    call free_slot_for_split(mass,mass_w,density,temperature,n,1, &
                             mass_base,smb_ice,runoff,t_srf,albedo, &
                             3,mass_max,c)
    call check("massless bottom layer just dropped", n .eq. 2, nfail)
    call check_acc("nothing exported", mass_base, 0.0_wp_acc, nfail)
    call check_conserve("massless slot freeing conserves mass", total_ref, nfail)

    ! =====================================================================
    ! enforce_snow_depth_cap
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- enforce_snow_depth_cap ---"

    ! reference_depth is a constant 22.5 m (Chion.jl 03bb445)
    call check_acc("reference depth is 22.5 m", &
                   real(BESSI_REFERENCE_SNOW_DEPTH_M,wp_acc), 22.5_wp_acc, nfail)

    ! Below the cap -> no-op.
    call clear_column()
    call clear_accum()
    n = 4
    mass(1:4)    = 300.0_wp
    density(1:4) = 300.0_wp      ! 4 m total
    total_ref = column_total()
    call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                mass_base,smb_ice,runoff,t_srf,albedo,c)
    call check("shallow column untouched", n .eq. 4, nfail)
    call check_acc("no basal export below the cap", mass_base, 0.0_wp_acc, nfail)
    call check_conserve("depth cap no-op conserves mass", total_ref, nfail)

    ! Above the cap: 12 layers of 300 kg m-2 at 300 kg m-3 = 12 m each... use
    ! 10 layers x 300/100 = 3 m each = 30 m, cap 22.5 m, excess 7.5 m at
    ! 100 kg m-3 = 750 kg m-2 to export.
    call clear_column()
    call clear_accum()
    n = 10
    mass(1:10)    = 300.0_wp
    mass_w(1:10)  = 5.0_wp
    density(1:10) = 100.0_wp
    temperature(1:10) = 265.0_wp

    total_ref = column_total()

    call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                mass_base,smb_ice,runoff,t_srf,albedo,c)

    call check_acc("exported exactly the excess mass", mass_base, 750.0_wp_acc, nfail)
    call check_acc("smb_ice credited identically", smb_ice, mass_base, nfail)
    call check("two full layers plus part of a third removed", n .eq. 8, nfail)
    call check_acc("remaining depth equals the cap", column_depth(), 22.5_wp_acc, nfail)
    call check_conserve("depth cap conserves mass", total_ref, nfail)

    ! Zero-density layer is exported in full without reducing the depth demand.
    call clear_column()
    call clear_accum()
    n = 9
    mass(1:9)    = 300.0_wp
    density(1:9) = 100.0_wp      ! 3 m each -> 27 m, excess 4.5 m
    density(9)   = 0.0_wp        ! contributes no depth, exported in full
    temperature(1:9) = 265.0_wp

    total_ref = column_total()

    call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                mass_base,smb_ice,runoff,t_srf,albedo,c)

    ! Depth is 8 * 3 = 24 m, cap 22.5 m -> excess 1.5 m.
    ! Reverse pass: layer 9 (rho = 0) exports 300 in full, then layer 8 covers
    ! 1.5 m at 100 kg m-3 = 150. Total 450.
    call check_acc("zero-density layer exported without reducing demand", &
                   mass_base, 450.0_wp_acc, nfail)
    call check_conserve("zero-density depth cap conserves mass", total_ref, nfail)

    ! =====================================================================
    ! Ntot == 1 and Ntot == 2 driven sequences
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- driven sequence, Ntot = 1 ---"

    call clear_column()
    call clear_accum()
    n = 0
    total_ref = 0.0_wp_acc
    albedo = c%alpha_dry

    do istep = 1, 40

        ! Accumulate.
        if (n .eq. 0) then
            n = 1
            mass(1) = 0.0_wp ; mass_w(1) = 0.0_wp
            density(1) = 300.0_wp ; temperature(1) = 265.0_wp
        end if
        mass(1)   = mass(1) + 60.0_wp
        mass_w(1) = mass_w(1) + 2.0_wp
        total_ref = total_ref + 62.0_wp_acc

        ! With Ntot = 1 the split can never happen; free_slot handles overflow.
        if (n .ge. 1 .and. mass(1) .gt. mass_max) then
            call free_slot_for_split(mass,mass_w,density,temperature,n,1, &
                                     mass_base,smb_ice,runoff,t_srf,albedo, &
                                     1,mass_max,c)
        end if

        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo,c)

        ! Melt-like conversion of solid to liquid every third step. Mass moves
        ! between mass and mass_w, so the conserved total is unchanged.
        if (mod(istep,3) .eq. 0 .and. n .ge. 1) then
            if (mass(1) .gt. 150.0_wp) then
                mass(1)   = mass(1) - 150.0_wp
                mass_w(1) = mass_w(1) + 150.0_wp
            end if
            call merge_layer(mass,mass_w,density,temperature,n,1, &
                             mass_split,mass_min,c)
        end if

        call check_conserve_quiet("Ntot=1 sequence conserves mass", total_ref, nfail, istep)

    end do
    call check("Ntot=1 sequence completed", .TRUE., nfail)
    call check_conserve("Ntot=1 final state conserves mass", total_ref, nfail)

    write(*,*)
    write(*,"(a)") "--- driven sequence, Ntot = 2 ---"

    call clear_column()
    call clear_accum()
    n = 0
    total_ref = 0.0_wp_acc
    albedo = c%alpha_dry

    do istep = 1, 60

        if (n .eq. 0) then
            n = 1
            mass(1) = 0.0_wp ; mass_w(1) = 0.0_wp
            density(1) = 320.0_wp ; temperature(1) = 264.0_wp
        end if
        mass(1)   = mass(1) + 80.0_wp
        mass_w(1) = mass_w(1) + 1.0_wp
        total_ref = total_ref + 81.0_wp_acc

        ! Accumulation's split loop, Ntot <= 2 variant: free a slot rather than
        ! merging the bottom pair (there is no pair to merge).
        do while (mass(1) .gt. mass_max)
            if (n .ge. 2) then
                call free_slot_for_split(mass,mass_w,density,temperature,n,1, &
                                         mass_base,smb_ice,runoff,t_srf,albedo, &
                                         2,mass_max,c)
            end if
            if (n .lt. 2) then
                call split_layer(mass,mass_w,density,temperature,n,1,2, &
                                 mass_max,mass_split)
            else
                exit
            end if
        end do

        do while (n .gt. 1 .and. mass(1) .lt. mass_min)
            call merge_layer(mass,mass_w,density,temperature,n,1, &
                             mass_split,mass_min,c)
        end do

        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo,c)

        call check_conserve_quiet("Ntot=2 sequence conserves mass", total_ref, nfail, istep)

    end do
    call check("Ntot=2 sequence completed", .TRUE., nfail)
    call check_conserve("Ntot=2 final state conserves mass", total_ref, nfail)
    call check("Ntot=2 column reached capacity", n .eq. 2, nfail)
    call check("Ntot=2 exported mass at the base", mass_base .gt. 0.0_wp_acc, nfail)

    ! =====================================================================
    ! Full accumulation/melt sequence at Ntot = 15
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- driven accumulation/melt sequence, Ntot = 15 ---"

    call clear_column()
    call clear_accum()
    n = 0
    total_ref = 0.0_wp_acc
    albedo = c%alpha_dry

    do istep = 1, 200

        ! --- Accumulation of fresh snow onto the surface layer.
        if (n .eq. 0) then
            n = 1
            mass(1) = 0.0_wp ; mass_w(1) = 0.0_wp
            density(1) = 315.0_wp ; temperature(1) = 262.0_wp
        end if
        mass(1)   = mass(1) + 90.0_wp
        total_ref = total_ref + 90.0_wp_acc

        ! --- Split loop, mirroring accumulation.jl's structure.
        do while (mass(1) .gt. mass_max)
            if (n .ge. 15) then
                call merge_bottom_layer(mass,mass_w,density,temperature,n,mass_base,smb_ice,c)
            end if
            if (n .lt. 15) then
                call split_layer(mass,mass_w,density,temperature,n,1,15, &
                                 mass_max,mass_split)
            else
                exit
            end if
        end do

        ! --- Merge loop.
        do while (n .gt. 1 .and. mass(1) .lt. mass_min)
            call merge_layer(mass,mass_w,density,temperature,n,1, &
                             mass_split,mass_min,c)
        end do

        ! --- Densify a little, so the depth cap and the rho_i export can bite.
        do k = 1, n
            density(k) = min(density(k)*1.02_wp,c%rho_i*1.05_wp)
        end do

        ! --- Depth cap.
        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo,c)

        ! --- Melt-like water production every 5th step: convert solid to
        !     liquid in place (mass leaves `mass`, appears in `mass_w`), which
        !     keeps the total budget unchanged.
        if (mod(istep,5) .eq. 0 .and. n .ge. 1) then
            if (mass(1) .gt. 40.0_wp) then
                mass(1)   = mass(1) - 40.0_wp
                mass_w(1) = mass_w(1) + 40.0_wp
            end if
            call merge_layer(mass,mass_w,density,temperature,n,1, &
                             mass_split,mass_min,c)
        end if

        call check_conserve_quiet("Ntot=15 sequence conserves mass", total_ref, nfail, istep)

    end do

    call check("Ntot=15 sequence reached capacity", n .eq. 15, nfail)
    call check("Ntot=15 sequence exported basal mass", mass_base .gt. 0.0_wp_acc, nfail)
    call check_conserve("Ntot=15 final state conserves mass", total_ref, nfail)

    total_now = column_total()
    write(*,"(a,g16.8)") "    final total (layers+base+runoff) = ", total_now
    write(*,"(a,g16.8)") "    expected                        = ", total_ref
    write(*,"(a,g16.8)") "    relative error                  = ", &
                         abs(total_now-total_ref)/max(abs(total_ref),1.0_wp_acc)

    ! =====================================================================
    ! Fine near-surface layers: remesh (Chion.jl 03bb445)
    ! =====================================================================
    write(*,"(a)") "--- near-surface remesh: cap down, fill up ---"

    h_fine = [0.02_wp, 0.05_wp, 0.10_wp, 0.30_wp]
    h_off  = 0.0_wp

    call check("near_surface_layer_count: all four limited", &
               near_surface_layer_count(h_fine) .eq. 4, nfail)
    call check("near_surface_layer_count: none", near_surface_layer_count(h_off) .eq. 0, nfail)
    call check("near_surface_layer_count: top layer only", &
               near_surface_layer_count([0.02_wp,0.0_wp,0.0_wp,0.0_wp]) .eq. 1, nfail)

    ! --- Chion.jl test_geometric_accumulation: a column already at its
    !     target thicknesses is left alone by zero-input accumulation (with
    !     the surface merge off, mass_min = 0) followed by the remesh.
    call clear_column()
    call clear_accum()
    n = 5
    mass(1:5)        = [6.0_wp, 15.0_wp, 30.0_wp, 90.0_wp, 300.0_wp]
    density(1:5)     = 300.0_wp
    temperature(1:5) = [250.0_wp, 255.0_wp, 260.0_wp, 265.0_wp, 270.0_wp]
    mass_before = mass
    temp_before = temperature
    snow_age    = 0.0_wp
    call apply_accumulation(mass,mass_w,density,temperature,n, &
                            mass_base,smb_ice,runoff,t_srf,albedo,snow_age, &
                            c,8,mass_max,mass_split,0.0_wp, &
                            0.0_wp,0.0_wp,86400.0_wp,250.0_wp,5.0_wp)
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    8,mass_max,mass_split,mass_min,h_fine,c)
    call check("geometric accumulation: n unchanged", n .eq. 5, nfail)
    call check("geometric accumulation: masses unchanged (rel 1e-6)", &
               all(abs(mass(1:5) - mass_before(1:5)) .le. 1.0e-6_wp*mass_before(1:5)), nfail)
    call check("geometric accumulation: temperatures unchanged (rel 1e-6)", &
               all(abs(temperature(1:5) - temp_before(1:5)) .le. 1.0e-6_wp*temp_before(1:5)), nfail)

    ! --- Cap down: a 0.5 m fresh surface layer over 1 m of firn is pushed
    !     down through all four fine layers, opening layers 3-5 as it goes.
    call clear_column()
    call clear_accum()
    n = 2
    mass(1:2)        = [50.0_wp, 400.0_wp]
    mass_w(1:2)      = [5.0_wp, 2.0_wp]
    density(1:2)     = [100.0_wp, 400.0_wp]
    temperature(1:2) = [260.0_wp, 270.0_wp]
    total_ref = column_total()
    vol_ref   = column_volume()
    ent_ref   = column_enthalpy()
    water_ref = sum(real(mass_w,wp_acc))

    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_fine,c)

    call check("cap: four fine layers plus one below", n .eq. 5, nfail)
    at_target = .TRUE.
    do k = 1, 4
        at_target = at_target .and. abs(mass(k)/density(k) - h_fine(k)) .le. 1.0e-5_wp*h_fine(k)
    end do
    call check("cap: layers 1-4 at their target thickness", at_target, nfail)
    call check("cap: surface layer keeps its density and temperature", &
               density(1) .eq. 100.0_wp .and. temperature(1) .eq. 260.0_wp, nfail)
    call check("cap: water moves with the mass (layer 1 keeps 5*2/50)", &
               abs(mass_w(1) - 0.2_wp) .le. 1.0e-5_wp, nfail)
    call check_conserve("cap conserves mass + water", total_ref, nfail)
    call check_acc("cap conserves water",    sum(real(mass_w,wp_acc)), water_ref, nfail)
    call check_acc("cap conserves volume",   column_volume(),   vol_ref, nfail)
    call check_acc("cap conserves sum(m*T)", column_enthalpy(), ent_ref, nfail)

    ! --- Fill up: melt-like loss of most of layer 1 (solid -> liquid in
    !     place) is replaced from below, layer by layer, down to layer 5.
    mass(1)   = mass(1) - 1.5_wp
    mass_w(1) = mass_w(1) + 1.5_wp
    total_ref = column_total()
    vol_ref   = column_volume()
    ent_ref   = column_enthalpy()
    water_ref = sum(real(mass_w,wp_acc))
    mass_before = mass

    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_fine,c)

    call check("fill: layer count unchanged", n .eq. 5, nfail)
    at_target = .TRUE.
    do k = 1, 4
        at_target = at_target .and. abs(mass(k)/density(k) - h_fine(k)) .le. 1.0e-5_wp*h_fine(k)
    end do
    call check("fill: layers 1-4 back at their target thickness", at_target, nfail)
    call check("fill: layer 5 gave the mass", mass(5) .lt. mass_before(5), nfail)
    call check_conserve("fill conserves mass + water", total_ref, nfail)
    call check_acc("fill conserves water",    sum(real(mass_w,wp_acc)), water_ref, nfail)
    call check_acc("fill conserves volume",   column_volume(),   vol_ref, nfail)
    call check_acc("fill conserves sum(m*T)", column_enthalpy(), ent_ref, nfail)

    ! --- Fill exhausts a donor: layer 2 is emptied into layer 1 and removed,
    !     and layer 1 goes on filling from the former layer 3.
    call clear_column()
    call clear_accum()
    n = 3
    mass(1:3)        = [0.2_wp, 1.0_wp, 100.0_wp]
    mass_w(1:3)      = [0.0_wp, 0.5_wp, 1.0_wp]
    density(1:3)     = 300.0_wp
    temperature(1:3) = [265.0_wp, 268.0_wp, 271.0_wp]
    total_ref = column_total()
    ent_ref   = column_enthalpy()

    call fill_near_surface_layer_thicknesses(mass,mass_w,density,temperature,n,h_fine,c)

    call check("fill: exhausted donor removed (n = 2)", n .eq. 2, nfail)
    call check_val("fill: layer 1 at 0.02 m", mass(1)/density(1), 0.02_wp, nfail)
    call check("fill: vacated slot 3 reset", mass(3) .eq. 0.0_wp .and. mass_w(3) .eq. 0.0_wp, nfail)
    call check_conserve("fill with exhausted donor conserves mass + water", total_ref, nfail)
    call check_acc("fill with exhausted donor conserves sum(m*T)", column_enthalpy(), ent_ref, nfail)

    ! --- A shallow column is not padded: layer 1 takes all of layer 2 and
    !     stays thinner than its target.
    call clear_column()
    n = 2
    mass(1:2)    = [1.0_wp, 2.0_wp]
    density(1:2) = 300.0_wp
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_fine,c)
    call check("shallow column: one thin layer left", &
               n .eq. 1 .and. abs(mass(1) - 3.0_wp) .le. 1.0e-6_wp, nfail)

    ! --- A full column keeps the excess in its deepest fine layer.
    call clear_column()
    n = 3
    mass(1:3)    = [6.0_wp, 15.0_wp, 200.0_wp]
    density(1:3) = 300.0_wp
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    3,mass_max,mass_split,mass_min,h_fine,c)
    call check("full column (n = Ntot): deepest fine layer keeps its excess", &
               n .eq. 3 .and. mass(3) .eq. 200.0_wp, nfail)

    ! --- No limits: the remesh is a no-op, bit for bit.
    call clear_column()
    n = 2
    mass(1:2)        = [50.0_wp, 1.0_wp]
    density(1:2)     = [100.0_wp, 400.0_wp]
    temperature(1:2) = [260.0_wp, 270.0_wp]
    mass_before = mass
    temp_before = temperature
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_off,c)
    call check("no limits: remesh is a no-op", n .eq. 2 .and. &
               all(mass .eq. mass_before) .and. all(temperature .eq. temp_before), nfail)

    write(*,*)

    ! =====================================================================
    ! Below the fine layers: split and merge by mass (C4b, D32)
    ! =====================================================================
    if (NEAR_SURFACE_SPLIT_MERGE_BELOW) then
        write(*,"(a)") "--- below the fine layers: split/merge of layer 5 (C4b) ---"
    else
        write(*,"(a)") "--- below the fine layers: layer 5 left alone (legacy_chion) ---"
    end if

    ! --- Split: layer 5 above mass_max.
    call fine_column([800.0_wp], 400.0_wp)
    total_ref = column_total()
    ent_ref   = column_enthalpy()
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_fine,c)
    if (NEAR_SURFACE_SPLIT_MERGE_BELOW) then
        call check("split below: layer 5 split (n = 6, 500 + 300)", n .eq. 6 .and. &
                   near(mass(5),500.0_wp) .and. near(mass(6),mass_split), nfail)
        call check("split below: both halves keep density and temperature", &
                   density(6) .eq. 400.0_wp .and. temperature(6) .eq. temperature(5), nfail)
    else
        call check("legacy: layer 5 not split", n .eq. 5 .and. near(mass(5),800.0_wp), nfail)
    end if
    call check("split below: fine layers untouched", near(mass(1),6.0_wp) .and. &
               near(mass(2),15.0_wp) .and. near(mass(3),30.0_wp) .and. near(mass(4),90.0_wp), nfail)
    call check_conserve("split below conserves mass + water", total_ref, nfail)
    call check_acc("split below conserves sum(m*T)", column_enthalpy(), ent_ref, nfail)

    ! --- Full merge: layer 5 below mass_min, layer 6 small.
    call fine_column([50.0_wp, 300.0_wp], 400.0_wp)
    total_ref = column_total()
    ent_ref   = column_enthalpy()
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_fine,c)
    if (NEAR_SURFACE_SPLIT_MERGE_BELOW) then
        call check("merge below: layers 5 and 6 merged (n = 5, 350)", &
                   n .eq. 5 .and. near(mass(5),350.0_wp) .and. mass(6) .eq. 0.0_wp, nfail)
    else
        call check("legacy: layer 5 not merged", n .eq. 6 .and. near(mass(5),50.0_wp), nfail)
    end if
    call check_conserve("merge below conserves mass + water", total_ref, nfail)
    call check_acc("merge below conserves sum(m*T)", column_enthalpy(), ent_ref, nfail)

    ! --- Partial transfer: layer 5 below mass_min, layer 6 large.
    call fine_column([50.0_wp, 700.0_wp], 400.0_wp)
    total_ref = column_total()
    ent_ref   = column_enthalpy()
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    10,mass_max,mass_split,mass_min,h_fine,c)
    if (NEAR_SURFACE_SPLIT_MERGE_BELOW) then
        call check("top-up below: layer 5 back to mass_split from layer 6", &
                   n .eq. 6 .and. near(mass(5),mass_split) .and. near(mass(6),450.0_wp), nfail)
    else
        call check("legacy: layer 5 not topped up", n .eq. 6 .and. near(mass(5),50.0_wp), nfail)
    end if
    call check_conserve("top-up below conserves mass + water", total_ref, nfail)
    call check_acc("top-up below conserves sum(m*T)", column_enthalpy(), ent_ref, nfail)

    ! --- At capacity (Ntot = 7): each split first merges the two deepest
    !     layers, until layer 5 is within bounds.
    call fine_column([900.0_wp, 300.0_wp, 300.0_wp], 400.0_wp)
    total_ref = column_total()
    call remesh_near_surface_layers(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    7,mass_max,mass_split,mass_min,h_fine,c)
    if (NEAR_SURFACE_SPLIT_MERGE_BELOW) then
        call check("at Ntot: layer 5 split after bottom merges (300, 300, 900)", &
                   n .eq. 7 .and. near(mass(5),300.0_wp) .and. near(mass(6),300.0_wp) &
                   .and. near(mass(7),900.0_wp), nfail)
    else
        call check("legacy: full column untouched", n .eq. 7 .and. near(mass(5),900.0_wp), nfail)
    end if
    call check_conserve("at Ntot conserves mass + water", total_ref, nfail)

    write(*,*)

    ! =====================================================================
    ! Summary
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP4: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP4: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    pure function near(x,y) result(ok)
        ! |x - y| <= 1e-5 |y|: equal up to sp round-off of the remesh.

        implicit none

        real(wp), intent(IN) :: x, y
        logical :: ok

        ok = abs(x - y) .le. 1.0e-5_wp*abs(y)

        return

    end function near

    subroutine fine_column(m_below,rho_below)
        ! Four fine layers exactly at their targets (density 300: 6, 15, 30,
        ! 90 kg m-2) over the given mass layers at density rho_below. Some
        ! water and a temperature gradient, so conservation checks bite.

        implicit none

        real(wp), intent(IN) :: m_below(:)          ! [kg m-2] layers 5, 6, ...
        real(wp), intent(IN) :: rho_below           ! [kg m-3]

        ! Local variables
        integer :: kk

        call clear_column()
        call clear_accum()

        n = 4 + size(m_below)
        mass(1:4)    = [6.0_wp, 15.0_wp, 30.0_wp, 90.0_wp]
        density(1:4) = 300.0_wp
        do kk = 1, size(m_below)
            mass(4+kk)    = m_below(kk)
            density(4+kk) = rho_below
        end do
        do kk = 1, n
            mass_w(kk)      = 0.01_wp*mass(kk)
            temperature(kk) = 250.0_wp + 2.0_wp*real(kk,wp)
        end do

        return

    end subroutine fine_column

    subroutine clear_column()

        implicit none

        mass        = 0.0_wp
        mass_w      = 0.0_wp
        density     = 0.0_wp
        temperature = c%T0
        n           = 0

        return

    end subroutine clear_column

    subroutine clear_accum()

        implicit none

        mass_base = 0.0_wp_acc
        smb_ice   = 0.0_wp_acc
        runoff    = 0.0_wp_acc
        t_srf     = c%T0
        albedo    = c%alpha_dry

        return

    end subroutine clear_accum

    function column_total() result(total)
        ! The conserved quantity: solid + liquid still in the column, plus
        ! everything exported to the base and to runoff. Accumulated in wp_acc.

        implicit none

        real(wp_acc) :: total

        ! Local variables
        integer :: kk

        total = mass_base + runoff

        do kk = 1, NMAX
            total = total + real(mass(kk),wp_acc) + real(mass_w(kk),wp_acc)
        end do

        return

    end function column_total

    function column_volume() result(vol)
        ! Total volume of the active layers [m], in wp_acc.

        implicit none

        real(wp_acc) :: vol

        ! Local variables
        integer :: kk

        vol = 0.0_wp_acc

        do kk = 1, n
            vol = vol + real(mass(kk),wp_acc)/real(density(kk),wp_acc)
        end do

        return

    end function column_volume

    function column_enthalpy() result(ent)
        ! Sensible enthalpy of the solid, up to the factor ci: sum(m*T) over
        ! the active layers [kg K m-2], in wp_acc.

        implicit none

        real(wp_acc) :: ent

        ! Local variables
        integer :: kk

        ent = 0.0_wp_acc

        do kk = 1, n
            ent = ent + real(mass(kk),wp_acc)*real(temperature(kk),wp_acc)
        end do

        return

    end function column_enthalpy

    function column_depth() result(depth)

        implicit none

        real(wp_acc) :: depth

        ! Local variables
        integer :: kk

        depth = 0.0_wp_acc

        do kk = 1, n
            if (mass(kk) .gt. 0.0_wp .and. density(kk) .gt. 0.0_wp) then
                depth = depth + real(mass(kk),wp_acc)/real(density(kk),wp_acc)
            end if
        end do

        return

    end function column_depth

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

    subroutine check_val(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp),         intent(IN)    :: value
        real(wp),         intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp) :: tol

        tol = max(abs(expected)*1.0e-6_wp,1.0e-20_wp)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g14.6)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g14.6,a,g14.6)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_val

    subroutine check_acc(label,value,expected,nfail)

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp_acc) :: tol

        ! Inputs originate in sp, so sp accuracy is the achievable target.
        tol = max(abs(expected)*1.0e-6_wp_acc,1.0e-20_wp_acc)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g16.8)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_acc

    subroutine check_conserve(label,expected,nfail)
        ! sum(mass) + sum(mass_w) + mass_base + runoff conserved to 1e-6
        ! relative, accumulated in wp_acc.

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp_acc) :: total, rel

        total = column_total()
        rel   = abs(total-expected)/max(abs(expected),1.0_wp_acc)

        if (rel .le. 1.0e-6_wp_acc) then
            write(*,"(a,a,a,es10.2)") "  ok   : ", trim(label), "   rel err = ", rel
        else
            write(*,"(a,a,a,es10.2,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                    "   rel err = ", rel, "  got ", total, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_conserve

    subroutine check_conserve_quiet(label,expected,nfail,istep)
        ! As check_conserve, but only reports on failure -- used inside the
        ! long driven sequences so the output stays readable.

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: expected
        integer,          intent(INOUT) :: nfail
        integer,          intent(IN)    :: istep

        ! Local variables
        real(wp_acc) :: total, rel

        total = column_total()
        rel   = abs(total-expected)/max(abs(expected),1.0_wp_acc)

        if (rel .gt. 1.0e-6_wp_acc) then
            write(*,"(a,a,a,i0,a,es10.2,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                    "  step ", istep, "  rel err = ", rel, &
                    "  got ", total, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_conserve_quiet

end program test_layers
