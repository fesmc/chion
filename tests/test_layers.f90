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

    use chion_defs,  only : wp, wp_acc, chion_const_class, chion_const_init, &
                            BESSI_REFERENCE_LAYER_COUNT, BESSI_REFERENCE_DEPTH_DENSITY
    use snow_layers

    implicit none

    integer, parameter :: NMAX = 20

    type(chion_const_class) :: c
    integer  :: nfail

    ! Column state, sized generously; the "Ntot" in use is passed per test.
    real(wp) :: mass(NMAX), mass_w(NMAX), density(NMAX), temperature(NMAX)
    integer  :: n

    real(wp_acc) :: mass_base, smb_ice, runoff
    real(wp)     :: t_srf, albedo

    real(wp_acc) :: total_ref, total_now
    real(wp_acc) :: ice_out, run_out

    integer  :: k, istep
    real(wp) :: mass_max, mass_split, mass_min

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
    ! split_surface_layer
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- split_surface_layer ---"

    call clear_column()
    n = 2
    mass(1) = 700.0_wp ; mass_w(1) = 70.0_wp
    density(1) = 350.0_wp ; temperature(1) = 265.0_wp
    mass(2) = 250.0_wp ; mass_w(2) = 5.0_wp
    density(2) = 500.0_wp ; temperature(2) = 268.0_wp

    total_ref = column_total()

    call split_surface_layer(mass,mass_w,density,temperature,n,5,mass_max,mass_split)

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
    call split_surface_layer(mass,mass_w,density,temperature,n,3,mass_max,mass_split)
    call check("no split when n == Ntot", n .eq. 3 .and. mass(1) .eq. 700.0_wp, nfail)

    ! Below mass_max -> no-op.
    call clear_column()
    n = 1
    mass(1) = 400.0_wp ; density(1) = 350.0_wp
    call split_surface_layer(mass,mass_w,density,temperature,n,5,mass_max,mass_split)
    call check("no split below mass_max", n .eq. 1 .and. mass(1) .eq. 400.0_wp, nfail)

    ! =====================================================================
    ! merge_surface_layer: partial-transfer branch
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- merge_surface_layer (partial transfer) ---"

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
    call merge_surface_layer(mass,mass_w,density,temperature,n,mass_split,mass_min,c)

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
    ! merge_surface_layer: full-merge branch
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- merge_surface_layer (full merge) ---"

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
    call merge_surface_layer(mass,mass_w,density,temperature,n,mass_split,mass_min,c)

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
    call merge_surface_layer(mass,mass_w,density,temperature,n,mass_split,mass_min,c)
    call check("n=1 below TOL_EMPTY_LAYER collapses to n=0", n .eq. 0, nfail)

    call clear_column()
    n = 1
    mass(1) = 1.0e-9_wp      ! above TOL_EMPTY_LAYER
    density(1) = 300.0_wp
    call merge_surface_layer(mass,mass_w,density,temperature,n,mass_split,mass_min,c)
    call check("n=1 above TOL_EMPTY_LAYER survives", n .eq. 1, nfail)

    ! Surface already at/above mass_min -> no-op.
    call clear_column()
    n = 2
    mass(1) = 150.0_wp ; mass(2) = 150.0_wp
    density(1:2) = 300.0_wp
    call merge_surface_layer(mass,mass_w,density,temperature,n,mass_split,mass_min,c)
    call check("no merge when surface >= mass_min", n .eq. 2, nfail)

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
            call split_surface_layer(mass,mass_w,density,temperature,n,5,mass_max,mass_split)
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
    ! free_slot_for_surface_split: Ntot == 1 special case
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- free_slot_for_surface_split (Ntot == 1) ---"

    call clear_column()
    call clear_accum()
    n = 1
    mass(1) = 800.0_wp ; mass_w(1) = 16.0_wp
    density(1) = 400.0_wp ; temperature(1) = 265.0_wp

    total_ref = column_total()

    call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
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
    call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
                                     mass_base,smb_ice,runoff,t_srf,albedo, &
                                     1,mass_max,c)
    call check_acc("no export when below mass_max", mass_base, 0.0_wp_acc, nfail)
    call check_conserve("Ntot=1 no-op conserves mass", total_ref, nfail)

    ! =====================================================================
    ! free_slot_for_surface_split: Ntot > 1
    ! =====================================================================
    write(*,*)
    write(*,"(a)") "--- free_slot_for_surface_split (Ntot == 2) ---"

    call clear_column()
    call clear_accum()
    n = 2
    mass(1) = 600.0_wp ; mass_w(1) = 6.0_wp ; density(1) = 350.0_wp
    mass(2) = 250.0_wp ; mass_w(2) = 5.0_wp ; density(2) = 500.0_wp

    total_ref = column_total()

    call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
                                     mass_base,smb_ice,runoff,t_srf,albedo, &
                                     2,mass_max,c)

    call check("bottom layer consumed, n = 1", n .eq. 1, nfail)
    call check_acc("bottom mass exported basally", mass_base, 250.0_wp_acc, nfail)
    call check_acc("bottom water to runoff", runoff, 5.0_wp_acc, nfail)
    call check_conserve("Ntot=2 slot freeing conserves mass", total_ref, nfail)

    ! ...and the split can then proceed.
    call split_surface_layer(mass,mass_w,density,temperature,n,2,mass_max,mass_split)
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

    call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
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

    ! reference_depth = 15 * 300 * 1.5 / 300 = 22.5 m
    call check_acc("reference depth is 22.5 m for mass_split = 300", &
                   real(BESSI_REFERENCE_LAYER_COUNT,wp_acc)*300.0_wp_acc*1.5_wp_acc &
                   /real(BESSI_REFERENCE_DEPTH_DENSITY,wp_acc), 22.5_wp_acc, nfail)

    ! Below the cap -> no-op.
    call clear_column()
    call clear_accum()
    n = 4
    mass(1:4)    = 300.0_wp
    density(1:4) = 300.0_wp      ! 4 m total
    total_ref = column_total()
    call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                mass_base,smb_ice,runoff,t_srf,albedo,mass_split,c)
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
                                mass_base,smb_ice,runoff,t_srf,albedo,mass_split,c)

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
                                mass_base,smb_ice,runoff,t_srf,albedo,mass_split,c)

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
            call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
                                             mass_base,smb_ice,runoff,t_srf,albedo, &
                                             1,mass_max,c)
        end if

        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo,mass_split,c)

        ! Melt-like conversion of solid to liquid every third step. Mass moves
        ! between mass and mass_w, so the conserved total is unchanged.
        if (mod(istep,3) .eq. 0 .and. n .ge. 1) then
            if (mass(1) .gt. 150.0_wp) then
                mass(1)   = mass(1) - 150.0_wp
                mass_w(1) = mass_w(1) + 150.0_wp
            end if
            call merge_surface_layer(mass,mass_w,density,temperature,n, &
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
                call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
                                                 mass_base,smb_ice,runoff,t_srf,albedo, &
                                                 2,mass_max,c)
            end if
            if (n .lt. 2) then
                call split_surface_layer(mass,mass_w,density,temperature,n,2, &
                                         mass_max,mass_split)
            else
                exit
            end if
        end do

        do while (n .gt. 1 .and. mass(1) .lt. mass_min)
            call merge_surface_layer(mass,mass_w,density,temperature,n, &
                                     mass_split,mass_min,c)
        end do

        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo,mass_split,c)

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
                call split_surface_layer(mass,mass_w,density,temperature,n,15, &
                                         mass_max,mass_split)
            else
                exit
            end if
        end do

        ! --- Merge loop.
        do while (n .gt. 1 .and. mass(1) .lt. mass_min)
            call merge_surface_layer(mass,mass_w,density,temperature,n, &
                                     mass_split,mass_min,c)
        end do

        ! --- Densify a little, so the depth cap and the rho_i export can bite.
        do k = 1, n
            density(k) = min(density(k)*1.02_wp,c%rho_i*1.05_wp)
        end do

        ! --- Depth cap.
        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo,mass_split,c)

        ! --- Melt-like water production every 5th step: convert solid to
        !     liquid in place (mass leaves `mass`, appears in `mass_w`), which
        !     keeps the total budget unchanged.
        if (mod(istep,5) .eq. 0 .and. n .ge. 1) then
            if (mass(1) .gt. 40.0_wp) then
                mass(1)   = mass(1) - 40.0_wp
                mass_w(1) = mass_w(1) + 40.0_wp
            end if
            call merge_surface_layer(mass,mass_w,density,temperature,n, &
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
