program test_bessi
    ! WP8 acceptance test: the assembled BESSI column kernel.
    !
    ! The point of this test is the ORDER of operations, not the physics of
    ! any single kernel (those are covered by the WP3-WP7 tests). So the
    ! checks are built around invariants that only hold if the assembly is
    ! right: mass closure, the bare-ice early return, and the diurnal wrapper.
    !
    ! ---------------------------------------------------------------------
    ! THE MASS-CLOSURE IDENTITY
    ! ---------------------------------------------------------------------
    ! Derived from the code, not assumed. Enumerate every mass mutation in
    ! bessi_column_step_core and where the mass comes from or goes:
    !
    !   accumulation   snowfall  -> mass(1)                       [source]
    !                  rainfall  -> mass_w(1), only if mass(1) > 0
    !                  bottom export -> mass_base AND smb_ice, in equal
    !                                   amounts; its liquid -> runoff
    !   vapor flux     mass(1) or mass_w(1) += vapor_mass         [source/sink]
    !   melt           mass -> mass_w, in-column; depleted layers route
    !                  their water to runoff
    !                  shortfall -> ice_melt: runoff += it, smb_ice -= it
    !   percolation    mass_w -> runoff
    !   refreezing     mass_w -> mass, in-column
    !   bare ice       melt_mass: runoff += it, and smb_ice += net_mass_change
    !                  = vapor_mass - melt_mass
    !   densification  no mass change
    !
    ! Let S = sum over active layers of (mass + mass_w). Splitting runoff into
    ! the part sourced from the column and the part sourced from the ice
    ! underneath (bare-ice melt and ice melt):
    !
    !   S = P + vapor_snow - (runoff - melt_bare - melt_ice) - mass_base
    !
    ! and by construction
    !
    !   smb_ice = mass_base + vapor_bare - melt_bare - melt_ice
    !
    ! Substituting, mass_base cancels and the two vapor terms recombine into
    ! the single vapor_mass accumulator:
    !
    !   ---------------------------------------------------------------
    !     S + runoff + smb_ice - vapor_mass  ==  P
    !   ---------------------------------------------------------------
    !
    ! with P the cumulative precipitation mass actually accepted by the
    ! column. Note mass_base does NOT appear: it is already inside smb_ice.
    !
    ! Two upstream defects would break this identity, and the test is set up
    ! to avoid both rather than to hide them:
    !
    !   * defect 11 -- rain falling on a column with mass(1) <= 0 is silently
    !     dropped, and rain is ignored entirely on the bare-ice path. The
    !     driver therefore withholds rain on any step that begins with
    !     mass(1) <= 0, so every kilogram offered is a kilogram accepted.
    !   * defect  1 -- apply_snow_surface_vapor_mass_flux returns the
    !     unclipped vapor_mass while applying a clipped one, so when
    !     sublimation demand exceeds the surface layer the diagnostic
    !     overstates the mass removed. The closure runs therefore leave the
    !     humidity forcing off (has_relative_humidity = .FALSE., has_q_lh =
    !     .FALSE.), which makes the latent flux identically zero. A separate
    !     non-asserting probe reports the residual WITH humidity on, which is
    !     a direct measurement of defect 1.
    !
    ! Measured behaviour of the residual (gfortran -O2, wp = sp): the RELATIVE
    ! residual saturates rather than growing with run length --
    !     1 yr 6.4e-7 | 2 yr 8.2e-7 | 5 yr 9.2e-7 | 10 yr 9.5e-7
    !                 | 20 yr 9.7e-7 | 40 yr 9.8e-7
    ! so it is a bounded per-step relative bias (~8 ulp of sp, from rounding
    ! mass(1) = m_prev + m_added on a layer of order mass_max), not a random
    ! walk that eventually breaks the budget. It sits just inside the 1e-6
    ! relative threshold mandated by docs/PLAN.md section 3.1, with little
    ! headroom -- worth knowing before that threshold is tightened.

    use chion_defs
    use snow_column_utils,  only : total_snow_water_mass, surface_has_snow
    use snow_diurnal,       only : diurnal_substep_bounds
    use snow_bessi

    implicit none

    integer,  parameter :: NDAY_YEAR = 365
    real(wp), parameter :: PI_WP     = 3.14159265358979_wp

    integer :: nfail

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP8 acceptance test: snow_bessi"
    write(*,"(a)") "=========================================================="
    write(*,*)

    call test_mass_closure(nfail)
    call test_cold_dry_column(nfail)
    call test_bare_and_recover(nfail)
    call test_capacity(nfail)
    call test_bare_ice_skips_water(nfail)
    call test_diurnal(nfail)
    call test_scheme_matrix(nfail)
    call probe_defect_1(nfail)

    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP8: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP8: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

    ! =====================================================================
    ! Forcing
    ! =====================================================================

    subroutine annual_forcing(iday,c,forc)
        ! Synthetic annual cycle at 70 N: a cold accumulation season and a
        ! melt season strong enough to strip the column bare, so that a
        ! multi-year run exercises accumulation, densification, melt,
        ! percolation, refreezing, the ice-melt shortfall AND the bare-ice
        ! early return.

        implicit none

        integer,                        intent(IN)  :: iday   ! 1..365
        type(chion_const_class),        intent(IN)  :: c
        type(chion_step_forcing_class), intent(OUT) :: forc

        ! Local variables
        real(wp) :: doy, t_air

        doy = real(iday,wp)

        ! Coldest ~15 Jan, warmest ~mid July.
        t_air = 268.0_wp - 12.0_wp*cos(2.0_wp*PI_WP*(doy - 15.0_wp)/real(NDAY_YEAR,wp))

        call neutral_forcing(forc)

        forc%air_temperature = t_air
        forc%dt_days         = 1.0_wp

        if (t_air .lt. c%T0) then
            forc%snowfall_rate = 6.0e-5_wp
            forc%rainfall_rate = 0.0_wp
        else
            forc%snowfall_rate = 0.0_wp
            forc%rainfall_rate = 3.0e-5_wp
        end if

        forc%shortwave_down = max(300.0_wp*cos(2.0_wp*PI_WP*(doy - 197.0_wp) &
                                               /real(NDAY_YEAR,wp)),0.0_wp)

        forc%wind_speed         = 3.0_wp
        forc%latitude_deg       = 70.0_wp
        forc%day_of_year        = doy
        forc%solar_longitude_deg = 360.0_wp*(doy - 80.0_wp)/real(NDAY_YEAR,wp)

        return

    end subroutine annual_forcing

    subroutine neutral_forcing(forc)
        ! Everything off: no prescribed fluxes, no humidity, sea-level
        ! pressure. Callers switch on only what they need.

        implicit none

        type(chion_step_forcing_class), intent(OUT) :: forc

        forc%air_temperature = 273.15_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 0.0_wp
        forc%rainfall_rate   = 0.0_wp
        forc%shortwave_down  = 0.0_wp
        forc%wind_speed      = 0.0_wp

        forc%q_sw_net  = 0.0_wp
        forc%q_lw_down = 0.0_wp
        forc%q_sh      = 0.0_wp
        forc%q_lh      = 0.0_wp

        forc%has_q_sw_net  = .FALSE.
        forc%has_q_lw_down = .FALSE.
        forc%has_q_sh      = .FALSE.
        forc%has_q_lh      = .FALSE.

        forc%relative_humidity     = 0.0_wp
        forc%has_relative_humidity = .FALSE.

        forc%air_pressure          = 101325.0_wp
        forc%prescribed_albedo     = 0.0_wp
        forc%has_prescribed_albedo = .FALSE.

        forc%latitude_deg        = 0.0_wp
        forc%day_of_year         = 1.0_wp
        forc%solar_longitude_deg = 0.0_wp

        return

    end subroutine neutral_forcing

    ! =====================================================================
    ! Drivers
    ! =====================================================================

    subroutine run_annual_cycle(bsi,c,nyears,icol,precip)
        ! Drive one column through nyears of the synthetic annual cycle,
        ! accumulating the precipitation mass the column actually accepted.
        !
        ! Rain is withheld on any step that begins with mass(1) <= 0, because
        ! apply_accumulation would drop it (upstream defect 11) and the
        ! closure identity would then be measuring that defect rather than
        ! the assembly.

        implicit none

        type(bessi_class),       intent(INOUT) :: bsi
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(IN)    :: nyears
        integer,                 intent(IN)    :: icol
        real(wp_acc),            intent(INOUT) :: precip

        ! Local variables
        integer  :: iyr, iday
        real(wp) :: dt_seconds
        logical  :: rain_would_land
        type(chion_step_forcing_class) :: forc

        do iyr = 1, nyears
            do iday = 1, NDAY_YEAR

                call annual_forcing(iday,c,forc)

                rain_would_land = .FALSE.
                if (bsi%now%n_lay(icol) .gt. 0) then
                    if (bsi%now%mass(1,icol) .gt. 0.0_wp) rain_would_land = .TRUE.
                end if
                if (.not. rain_would_land) forc%rainfall_rate = 0.0_wp

                dt_seconds = forc%dt_days*c%seconds_per_day

                precip = precip + real(forc%snowfall_rate*dt_seconds,wp_acc) &
                                + real(forc%rainfall_rate*dt_seconds,wp_acc)

                call bessi_column_step(bsi,icol,forc,c)

            end do
        end do

        return

    end subroutine run_annual_cycle

    function column_storage(bsi,icol) result(s)
        ! Sum of solid and liquid mass held in the active layers, in dp.
        ! Deliberately NOT total_snow_water_mass: that clips negative layer
        ! masses, which would mask exactly the kind of bug this test hunts.

        implicit none

        type(bessi_class), intent(IN) :: bsi
        integer,           intent(IN) :: icol
        real(wp_acc) :: s

        ! Local variables
        integer :: k

        s = 0.0_wp_acc

        do k = 1, bsi%now%n_lay(icol)
            s = s + real(bsi%now%mass(k,icol),wp_acc) &
                  + real(bsi%now%mass_w(k,icol),wp_acc)
        end do

        return

    end function column_storage

    function closure_lhs(bsi,icol) result(lhs)
        ! S + runoff + smb_ice - vapor_mass. See the identity derived in the
        ! header. Must equal the accepted precipitation.

        implicit none

        type(bessi_class), intent(IN) :: bsi
        integer,           intent(IN) :: icol
        real(wp_acc) :: lhs

        lhs = column_storage(bsi,icol)       &
            + bsi%now%runoff(icol)           &
            + bsi%now%smb_ice(icol)          &
            - bsi%now%vapor_mass(icol)

        return

    end function closure_lhs

    ! =====================================================================
    ! Test 1 -- mass closure over a multi-year annual cycle
    ! =====================================================================

    subroutine test_mass_closure(nfail)

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        real(wp_acc)            :: precip

        write(*,"(a)") "--- 1. mass closure: S + runoff + smb_ice - vapor == precip ---"

        ! --- BESSI densification, dynamic albedo ---------------------------
        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        precip = 0.0_wp_acc
        call run_annual_cycle(bsi,c,5,1,precip)

        write(*,"(a,g16.8)")   "         precip accepted [kg m-2] = ", precip
        write(*,"(a,g16.8)")   "         storage         [kg m-2] = ", column_storage(bsi,1)
        write(*,"(a,g16.8)")   "         runoff          [kg m-2] = ", bsi%now%runoff(1)
        write(*,"(a,g16.8)")   "         smb_ice         [kg m-2] = ", bsi%now%smb_ice(1)
        write(*,"(a,g16.8)")   "         melt            [kg m-2] = ", bsi%now%melt(1)
        write(*,"(a,g16.8)")   "         refreezing      [kg m-2] = ", bsi%now%refreezing(1)
        write(*,"(a,g16.8)")   "         mass_base       [kg m-2] = ", bsi%now%mass_base(1)

        call check("5-yr cycle actually melted the column bare at least once", &
                   bsi%now%melt(1) .gt. precip, nfail)
        call check("5-yr cycle produced runoff", bsi%now%runoff(1) .gt. 0.0_wp_acc, nfail)

        call check_close("closure, bessi densification / dynamic albedo", &
                         closure_lhs(bsi,1),precip,1.0e-6_wp_acc,nfail)

        call bessi_dealloc(bsi)

        ! --- HTESSEL densification, dynamic albedo -------------------------
        ! Exercises the snapshot taken in step 7 and consumed in step 13.
        call chion_const_init(c)
        c%low_density_densification = CHION_DENSIFY_HTESSEL
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        precip = 0.0_wp_acc
        call run_annual_cycle(bsi,c,5,1,precip)

        call check_close("closure, htessel densification / dynamic albedo", &
                         closure_lhs(bsi,1),precip,1.0e-6_wp_acc,nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_mass_closure

    ! =====================================================================
    ! Test 2 -- cold dry column
    ! =====================================================================

    subroutine test_cold_dry_column(nfail)
        ! Constant snowfall at 250 K with no sunlight. Layers must build up
        ! and densify, and nothing may melt or run off.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer  :: istep
        real(wp) :: rho_top, rho_bot

        write(*,"(a)") "--- 2. cold dry column: accumulate + densify, no melt ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        call neutral_forcing(forc)
        forc%air_temperature = 250.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 6.0e-5_wp
        forc%wind_speed      = 3.0_wp

        do istep = 1, 3*NDAY_YEAR
            call bessi_column_step(bsi,1,forc,c)
        end do

        rho_top = bsi%now%density(1,1)
        rho_bot = bsi%now%density(bsi%now%n_lay(1),1)

        write(*,"(a,i0)")    "         n_lay             = ", bsi%now%n_lay(1)
        write(*,"(a,g14.6)") "         density top       = ", rho_top
        write(*,"(a,g14.6)") "         density bottom    = ", rho_bot
        write(*,"(a,g14.6)") "         t_srf             = ", bsi%now%t_srf(1)

        call check("layers accumulated", bsi%now%n_lay(1) .gt. 1, nfail)
        call check("no melt",       bsi%now%melt(1)   .eq. 0.0_wp_acc, nfail)
        call check("no runoff",     bsi%now%runoff(1) .eq. 0.0_wp_acc, nfail)
        call check("no refreezing", bsi%now%refreezing(1) .eq. 0.0_wp_acc, nfail)
        call check("no liquid water anywhere", &
                   maxval(bsi%now%mass_w(:,1)) .eq. 0.0_wp, nfail)
        call check("column densified with depth", rho_bot .gt. rho_top, nfail)
        call check("bottom density exceeds fresh-snow density", rho_bot .gt. c%rho_s, nfail)
        call check("density never exceeds ice", &
                   maxval(bsi%now%density(1:bsi%now%n_lay(1),1)) .le. c%rho_i, nfail)
        call check("surface stayed below freezing", bsi%now%t_srf(1) .lt. c%T0, nfail)
        call check("albedo relaxed toward dry-snow value", &
                   bsi%now%albedo(1) .le. c%alpha_dry .and. &
                   bsi%now%albedo(1) .ge. c%alpha_wet, nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_cold_dry_column

    ! =====================================================================
    ! Test 3 -- melt to bare, then recover
    ! =====================================================================

    subroutine test_bare_and_recover(nfail)

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer      :: istep, n_after_build, day_bare
        real(wp)     :: albedo_bare
        real(wp_acc) :: smb_before_bare, smb_after_bare, runoff_before_bare

        write(*,"(a)") "--- 3. melting column goes bare, then recovers ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        ! Phase A: build a snowpack.
        call neutral_forcing(forc)
        forc%air_temperature = 260.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 1.0e-4_wp
        forc%wind_speed      = 2.0_wp

        do istep = 1, 60
            call bessi_column_step(bsi,1,forc,c)
        end do

        n_after_build = bsi%now%n_lay(1)
        call check("phase A built a snowpack", n_after_build .gt. 0, nfail)

        ! Phase B: strip it.
        call neutral_forcing(forc)
        forc%air_temperature = 280.0_wp
        forc%dt_days         = 1.0_wp
        forc%shortwave_down  = 350.0_wp
        forc%wind_speed      = 2.0_wp

        day_bare = 0
        do istep = 1, 200
            call bessi_column_step(bsi,1,forc,c)
            if (day_bare .eq. 0) then
                if (.not. surface_has_snow(bsi%now%mass(:,1),bsi%now%n_lay(1))) day_bare = istep
            end if
        end do

        albedo_bare        = bsi%now%albedo(1)
        smb_before_bare    = bsi%now%smb_ice(1)
        runoff_before_bare = bsi%now%runoff(1)

        write(*,"(a,i0)")    "         first bare step   = ", day_bare
        write(*,"(a,g14.6)") "         albedo when bare  = ", albedo_bare
        write(*,"(a,g16.8)") "         smb_ice           = ", smb_before_bare

        call check("column went bare during phase B", day_bare .gt. 0, nfail)
        call check("column is empty after phase B", bsi%now%n_lay(1) .eq. 0, nfail)
        call check("albedo dropped to bare-ice value", &
                   abs(albedo_bare - c%alpha_ice) .lt. 1.0e-6_wp, nfail)
        call check("phase B produced runoff", runoff_before_bare .gt. 0.0_wp_acc, nfail)
        call check("bare ice is losing mass (smb_ice < 0)", &
                   smb_before_bare .lt. 0.0_wp_acc, nfail)

        ! Phase C: snow returns.
        call neutral_forcing(forc)
        forc%air_temperature = 263.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 1.0e-4_wp
        forc%wind_speed      = 2.0_wp

        call bessi_column_step(bsi,1,forc,c)

        write(*,"(a,g14.6)") "         T(1) after resnow = ", bsi%now%temperature(1,1)

        call check("one snowy step re-creates a layer", bsi%now%n_lay(1) .eq. 1, nfail)

        ! Step 3 rule: fresh snow landing on a column that STARTED bare takes
        ! the air temperature. The value cannot be checked exactly here
        ! because the energy solve runs afterwards, but the alternative is
        ! unambiguous: without the rule the layer would carry the reset value
        ! c%T0 = 273.15 K, and the solve cannot cool 8.6 kg m-2 of snow by ten
        ! degrees in one step. Anything below 268 K can only have come from
        ! the rule having fired.
        call check("fresh layer took the air temperature (step 3 rule)", &
                   bsi%now%temperature(1,1) .lt. 268.0_wp, nfail)
        call check("albedo brightened away from bare ice", &
                   bsi%now%albedo(1) .gt. c%alpha_ice, nfail)

        do istep = 1, 60
            call bessi_column_step(bsi,1,forc,c)
        end do

        smb_after_bare = bsi%now%smb_ice(1)

        call check("column recovered a multi-layer snowpack", bsi%now%n_lay(1) .gt. 1, nfail)
        call check("bare-ice ablation stopped once snow returned", &
                   abs(smb_after_bare - smb_before_bare) .lt. 1.0e-6_wp_acc*abs(smb_before_bare), &
                   nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_bare_and_recover

    ! =====================================================================
    ! Test 4 -- drive the column to Ntot capacity
    ! =====================================================================

    subroutine test_capacity(nfail)
        ! Heavy cold snowfall: the split loop fills every slot, then each
        ! further split must free one by merging the two deepest layers, and
        ! the depth cap must export from the bottom. n must never exceed Ntot.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer      :: istep, n_max
        real(wp_acc) :: precip, thickness
        integer      :: k

        write(*,"(a)") "--- 4. column driven to Ntot capacity ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        call neutral_forcing(forc)
        forc%air_temperature = 255.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 2.0e-3_wp        ! 172.8 kg m-2 d-1
        forc%wind_speed      = 4.0_wp

        precip = 0.0_wp_acc
        n_max  = 0

        do istep = 1, 400
            precip = precip + real(forc%snowfall_rate*forc%dt_days*c%seconds_per_day,wp_acc)
            call bessi_column_step(bsi,1,forc,c)
            n_max = max(n_max,bsi%now%n_lay(1))
            if (bsi%now%n_lay(1) .gt. bsi%par%Ntot) exit
        end do

        thickness = 0.0_wp_acc
        do k = 1, bsi%now%n_lay(1)
            thickness = thickness + real(bsi%now%mass(k,1),wp_acc) &
                                    /real(bsi%now%density(k,1),wp_acc)
        end do

        write(*,"(a,i0)")    "         n_lay max         = ", n_max
        write(*,"(a,g16.8)") "         mass_base         = ", bsi%now%mass_base(1)
        write(*,"(a,g16.8)") "         column thickness  = ", thickness

        call check("layer count reached Ntot", n_max .eq. bsi%par%Ntot, nfail)
        call check("layer count never exceeded Ntot", &
                   bsi%now%n_lay(1) .le. bsi%par%Ntot, nfail)
        call check("bottom export occurred (merge and/or depth cap)", &
                   bsi%now%mass_base(1) .gt. 0.0_wp_acc, nfail)
        call check("depth cap bounded the column", &
                   thickness .lt. 2.0_wp_acc*real(BESSI_REFERENCE_LAYER_COUNT,wp_acc) &
                                  *real(bsi%par%mass_split,wp_acc) &
                                  /real(BESSI_REFERENCE_DEPTH_DENSITY,wp_acc), nfail)
        call check_close("closure holds at capacity", &
                         closure_lhs(bsi,1),precip,1.0e-6_wp_acc,nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_capacity

    ! =====================================================================
    ! Test 5 -- the bare-ice early return
    ! =====================================================================

    subroutine test_bare_ice_skips_water(nfail)
        ! Two identical columns differing ONLY in the surface layer's solid
        ! mass, straddling TOL_EMPTY_LAYER. Column 1 is "bare" by
        ! surface_has_snow and must take the early return; column 2 is not.
        !
        ! The subsurface layer is wet (20 kg m-2 of liquid) and cold (263 K),
        ! so if percolation ran it would shed ~5.9 kg m-2 to runoff (pore
        ! volume 0.141 m, retention 14.1 kg m-2) and if refreezing ran it
        ! would freeze ~6.4 kg m-2 and pull the layer up to T0. Neither may
        ! happen on the bare column.
        !
        ! mass_min is lowered for this test. With the default 100 kg m-2 the
        ! merge loop inside apply_accumulation would immediately absorb the
        ! sliver of a surface layer into the layer below, and the column would
        ! no longer be bare by the time the step-5 test is reached -- which is
        ! itself worth knowing: on a default column the bare state with a wet
        ! layer underneath is only reachable through melt or sublimation, not
        ! through a thin surface layer.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc

        write(*,"(a)") "--- 5. bare-ice branch skips percolation and refreezing ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        bsi%par%mass_min = 1.0e-12_wp
        call bessi_par_validate(bsi%par)
        call bessi_alloc(bsi,2)
        call bessi_init_state(bsi,c)

        ! Column 1: surface mass below TOL_EMPTY_LAYER -> bare.
        bsi%now%n_lay(1)         = 2
        bsi%now%mass(1,1)        = 1.0e-11_wp
        bsi%now%mass(2,1)        = 100.0_wp
        bsi%now%mass_w(1,1)      = 0.0_wp
        bsi%now%mass_w(2,1)      = 20.0_wp
        bsi%now%density(1,1)     = 400.0_wp
        bsi%now%density(2,1)     = 400.0_wp
        bsi%now%temperature(1,1) = 263.0_wp
        bsi%now%temperature(2,1) = 263.0_wp

        ! Column 2: identical but with a real surface layer.
        bsi%now%n_lay(2)         = 2
        bsi%now%mass(1,2)        = 50.0_wp
        bsi%now%mass(2,2)        = 100.0_wp
        bsi%now%mass_w(1,2)      = 0.0_wp
        bsi%now%mass_w(2,2)      = 20.0_wp
        bsi%now%density(1,2)     = 400.0_wp
        bsi%now%density(2,2)     = 400.0_wp
        bsi%now%temperature(1,2) = 263.0_wp
        bsi%now%temperature(2,2) = 263.0_wp

        call check("setup: column 1 reads as bare", &
                   .not. surface_has_snow(bsi%now%mass(:,1),bsi%now%n_lay(1)), nfail)
        call check("setup: column 2 reads as snow-covered", &
                   surface_has_snow(bsi%now%mass(:,2),bsi%now%n_lay(2)), nfail)

        call neutral_forcing(forc)
        forc%air_temperature = 275.0_wp
        forc%dt_days         = 1.0_wp
        forc%shortwave_down  = 200.0_wp
        forc%wind_speed      = 2.0_wp

        call bessi_column_step(bsi,1,forc,c)
        call bessi_column_step(bsi,2,forc,c)

        write(*,"(a,g16.8)") "         bare  mass_w(2)   = ", bsi%now%mass_w(2,1)
        write(*,"(a,g16.8)") "         snowy mass_w(2)   = ", bsi%now%mass_w(2,2)

        ! The bare column: nothing may have touched the water or the heat.
        call check("bare column: liquid water untouched (no percolation)", &
                   bsi%now%mass_w(2,1) .eq. 20.0_wp, nfail)
        call check("bare column: nothing refroze", &
                   bsi%now%refreezing(1) .eq. 0.0_wp_acc, nfail)
        call check("bare column: subsurface temperature unchanged (no latent release)", &
                   bsi%now%temperature(2,1) .eq. 263.0_wp, nfail)
        call check("bare column: subsurface solid mass unchanged", &
                   bsi%now%mass(2,1) .eq. 100.0_wp, nfail)
        call check("bare column: took the bare-ice path (melt charged to runoff)", &
                   bsi%now%runoff(1) .gt. 0.0_wp_acc, nfail)
        call check("bare column: albedo forced to bare ice", &
                   abs(bsi%now%albedo(1) - c%alpha_ice) .lt. 1.0e-6_wp, nfail)

        ! The contrast case: with snow present, both DO run.
        call check("snowy column: water was processed (percolation and/or refreezing)", &
                   bsi%now%mass_w(2,2) .ne. 20.0_wp, nfail)
        call check("snowy column: refreezing occurred", &
                   bsi%now%refreezing(2) .gt. 0.0_wp_acc, nfail)
        call check("snowy column: latent release warmed the subsurface", &
                   bsi%now%temperature(2,2) .gt. 263.0_wp, nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_bare_ice_skips_water

    ! =====================================================================
    ! Test 6 -- diurnal substepping
    ! =====================================================================

    subroutine test_diurnal(nfail)
        ! 6a. A cold, sunlit, snowing column with no melt: substepping must
        !     not change the mass at all, because precipitation is a RATE and
        !     dt_days is what shrinks.
        ! 6b. A warm column: substepping resolves the midday shortwave peak,
        !     which melting rectifies, so the melt MUST increase.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi_off, bsi_on
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer      :: istep, k, n_sub
        real(wp_acc) :: s_off, s_on, melt_off, melt_on
        real(wp_acc) :: frac_sum, tol_roundoff
        real(wp)     :: h_a, h_b

        write(*,"(a)") "--- 6. diurnal substepping on vs off ---"

        ! --- 6a-i: the tiling itself ---------------------------------------
        ! Precipitation stays a RATE while dt_days is scaled by the interval
        ! fraction, so the daily precipitation total is conserved exactly if
        ! and only if the fractions sum to one. Checked directly, in dp,
        ! before any sp state is involved.
        n_sub    = 8
        frac_sum = 0.0_wp_acc
        do k = 1, n_sub
            call diurnal_substep_bounds(k,n_sub,h_a,h_b)
            frac_sum = frac_sum + (real(h_b,wp_acc) - real(h_a,wp_acc)) &
                                  /(2.0_wp_acc*real(PI_WP,wp_acc))
        end do

        call check_close("substep fractions sum to one day exactly", &
                         frac_sum,1.0_wp_acc,1.0e-12_wp_acc,nfail)

        ! --- 6a: cold, no melt -> mass must be identical -------------------
        call chion_const_init(c)

        call bessi_par_init(bsi_off%par)
        call bessi_alloc(bsi_off,1)
        call bessi_init_state(bsi_off,c)

        call bessi_par_init(bsi_on%par)
        bsi_on%par%diurnal_shortwave_substeps     = .TRUE.
        bsi_on%par%diurnal_shortwave_max_substeps = 8
        bsi_on%par%diurnal_shortwave_threshold    = 0.0_wp
        call bessi_par_validate(bsi_on%par)
        call bessi_alloc(bsi_on,1)
        call bessi_init_state(bsi_on,c)

        call neutral_forcing(forc)
        forc%air_temperature     = 266.0_wp    ! above the 265.15 K substep gate
        forc%dt_days             = 1.0_wp
        forc%snowfall_rate       = 6.0e-5_wp
        forc%shortwave_down      = 250.0_wp
        forc%wind_speed          = 2.0_wp
        forc%latitude_deg        = 70.0_wp
        forc%solar_longitude_deg = 90.0_wp     ! near northern summer solstice

        do istep = 1, 120
            call bessi_column_step(bsi_off,1,forc,c)
            call bessi_column_step(bsi_on,1,forc,c)
        end do

        s_off = column_storage(bsi_off,1)
        s_on  = column_storage(bsi_on,1)

        write(*,"(a,g16.8)") "         cold storage, off = ", s_off
        write(*,"(a,g16.8)") "         cold storage, on  = ", s_on

        ! The two runs differ only in HOW the identical daily snowfall mass is
        ! added: once, or in eight pieces. The difference that survives is
        ! pure sp store round-off on mass(1), and is bounded above by
        ! (number of additions) * eps(sp)/2 * (largest layer mass) -- the
        ! worst case in which every rounding goes the same way. Asserting
        ! against that derived bound rather than a round number keeps the
        ! check honest: it fails if anything other than round-off changes.
        tol_roundoff = 0.5_wp_acc*real(120*8,wp_acc)*real(epsilon(1.0_wp),wp_acc) &
                       *real(bsi_on%par%mass_max,wp_acc)

        write(*,"(a,es10.3)") "         sp round-off bound= ", tol_roundoff

        call check("cold case: no melt either way", &
                   bsi_off%now%melt(1) .eq. 0.0_wp_acc .and. &
                   bsi_on%now%melt(1)  .eq. 0.0_wp_acc, nfail)
        call check("daily precipitation total is conserved by substepping &
                   &(to the sp round-off bound)", &
                   abs(s_on - s_off) .le. tol_roundoff, nfail)

        call bessi_dealloc(bsi_off)
        call bessi_dealloc(bsi_on)

        ! --- 6b: warm, melting -> melt must differ -------------------------
        call chion_const_init(c)

        call bessi_par_init(bsi_off%par)
        call bessi_alloc(bsi_off,1)
        call bessi_init_state(bsi_off,c)

        call bessi_par_init(bsi_on%par)
        bsi_on%par%diurnal_shortwave_substeps     = .TRUE.
        bsi_on%par%diurnal_shortwave_max_substeps = 8
        bsi_on%par%diurnal_shortwave_threshold    = 0.0_wp
        call bessi_par_validate(bsi_on%par)
        call bessi_alloc(bsi_on,1)
        call bessi_init_state(bsi_on,c)

        ! Build the same snowpack in both, with substepping inactive
        ! (shortwave = 0 disables the substep criterion outright).
        call neutral_forcing(forc)
        forc%air_temperature = 260.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 3.0e-4_wp
        forc%wind_speed      = 2.0_wp

        do istep = 1, 60
            call bessi_column_step(bsi_off,1,forc,c)
            call bessi_column_step(bsi_on,1,forc,c)
        end do

        call check_close("identical snowpacks before the melt phase", &
                         column_storage(bsi_on,1),column_storage(bsi_off,1), &
                         1.0e-9_wp_acc,nfail)

        call neutral_forcing(forc)
        forc%air_temperature     = 272.0_wp
        forc%dt_days             = 1.0_wp
        forc%shortwave_down      = 220.0_wp
        forc%wind_speed          = 2.0_wp
        forc%latitude_deg        = 70.0_wp
        forc%solar_longitude_deg = 90.0_wp

        do istep = 1, 40
            call bessi_column_step(bsi_off,1,forc,c)
            call bessi_column_step(bsi_on,1,forc,c)
        end do

        melt_off = bsi_off%now%melt(1)
        melt_on  = bsi_on%now%melt(1)

        write(*,"(a,g16.8)") "         melt, substeps off= ", melt_off
        write(*,"(a,g16.8)") "         melt, substeps on = ", melt_on

        call check("warm case: substepping changed the melt", melt_on .ne. melt_off, nfail)
        call check("warm case: substepping increased the melt (rectification)", &
                   melt_on .gt. melt_off, nfail)

        call bessi_dealloc(bsi_off)
        call bessi_dealloc(bsi_on)

        write(*,*)

        return

    end subroutine test_diurnal

    ! =====================================================================
    ! Test 7 -- scheme matrix
    ! =====================================================================

    subroutine test_scheme_matrix(nfail)
        ! Both densification schemes x all three albedo schemes x both
        ! fresh-snow-density schemes, each over two annual cycles. Every
        ! combination must run, stay finite, and close.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        integer      :: idens, ialb, irho
        real(wp_acc) :: precip
        character(len=64) :: label
        character(len=16) :: dens_name(2), alb_name(3), rho_name(2)
        integer :: dens_flag(2), alb_flag(3), rho_flag(2)

        dens_name = ["bessi  ","htessel"]
        dens_flag = [CHION_DENSIFY_BESSI,CHION_DENSIFY_HTESSEL]

        alb_name  = ["constant  ","dynamic   ","prescribed"]
        alb_flag  = [CHION_ALBEDO_CONSTANT,CHION_ALBEDO_DYNAMIC,CHION_ALBEDO_PRESCRIBED]

        rho_name  = ["rho_const","rho_param"]
        rho_flag  = [CHION_FRESH_SNOW_DENSITY_CONSTANT,CHION_FRESH_SNOW_DENSITY_PARAMETERIZED]

        write(*,"(a)") "--- 7. every densification / albedo / fresh-snow scheme ---"

        do idens = 1, 2
        do ialb  = 1, 3
        do irho  = 1, 2

            call chion_const_init(c)
            c%low_density_densification    = dens_flag(idens)
            c%albedo_scheme                = alb_flag(ialb)
            c%fresh_snow_density_scheme    = rho_flag(irho)

            call bessi_par_init(bsi%par)
            call bessi_alloc(bsi,1)
            call bessi_init_state(bsi,c)

            precip = 0.0_wp_acc
            call run_annual_cycle_albedo(bsi,c,2,1,precip, &
                                         alb_flag(ialb) .eq. CHION_ALBEDO_PRESCRIBED)

            write(label,"(a,a,a,a,a,a)") trim(dens_name(idens)), " / ", &
                                         trim(alb_name(ialb)),   " / ", &
                                         trim(rho_name(irho)),   ""

            call check(trim(label)//": state stayed finite", &
                       all(bsi%now%mass(:,1)    .eq. bsi%now%mass(:,1))    .and. &
                       all(bsi%now%density(:,1) .eq. bsi%now%density(:,1)) .and. &
                       all(bsi%now%temperature(:,1) .eq. bsi%now%temperature(:,1)), nfail)

            call check_close(trim(label)//": closure", &
                             closure_lhs(bsi,1),precip,1.0e-6_wp_acc,nfail)

            call bessi_dealloc(bsi)

        end do
        end do
        end do

        write(*,*)

        return

    end subroutine test_scheme_matrix

    subroutine run_annual_cycle_albedo(bsi,c,nyears,icol,precip,prescribed)
        ! As run_annual_cycle, but optionally supplying a prescribed albedo so
        ! the PRESCRIBED scheme takes its intended path rather than silently
        ! behaving as dynamic (trap 9).

        implicit none

        type(bessi_class),       intent(INOUT) :: bsi
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(IN)    :: nyears
        integer,                 intent(IN)    :: icol
        real(wp_acc),            intent(INOUT) :: precip
        logical,                 intent(IN)    :: prescribed

        ! Local variables
        integer  :: iyr, iday
        real(wp) :: dt_seconds
        logical  :: rain_would_land
        type(chion_step_forcing_class) :: forc

        do iyr = 1, nyears
            do iday = 1, NDAY_YEAR

                call annual_forcing(iday,c,forc)

                if (prescribed) then
                    forc%has_prescribed_albedo = .TRUE.
                    forc%prescribed_albedo     = 0.65_wp
                end if

                rain_would_land = .FALSE.
                if (bsi%now%n_lay(icol) .gt. 0) then
                    if (bsi%now%mass(1,icol) .gt. 0.0_wp) rain_would_land = .TRUE.
                end if
                if (.not. rain_would_land) forc%rainfall_rate = 0.0_wp

                dt_seconds = forc%dt_days*c%seconds_per_day

                precip = precip + real(forc%snowfall_rate*dt_seconds,wp_acc) &
                                + real(forc%rainfall_rate*dt_seconds,wp_acc)

                call bessi_column_step(bsi,icol,forc,c)

            end do
        end do

        return

    end subroutine run_annual_cycle_albedo

    ! =====================================================================
    ! Probe -- upstream defect 1, measured rather than asserted
    ! =====================================================================

    subroutine probe_defect_1(nfail)
        ! With humidity forcing on, apply_snow_surface_vapor_mass_flux clips
        ! the mass it removes but reports the unclipped demand, so the closure
        ! identity acquires a residual exactly equal to the over-reported
        ! sublimation. This is an upstream defect, not a port bug, so it is
        ! measured and printed rather than asserted. The check that DOES run
        ! is a sign check: the residual can only ever go one way.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer      :: istep
        real(wp_acc) :: precip, residual

        write(*,"(a)") "--- probe. defect 1 (vapor-mass diagnostics not mass-closed) ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        call neutral_forcing(forc)
        forc%air_temperature       = 258.0_wp
        forc%dt_days               = 1.0_wp
        forc%snowfall_rate         = 2.0e-6_wp    ! barely any snow to sublimate
        forc%shortwave_down        = 0.0_wp
        forc%wind_speed            = 2.0_wp
        forc%relative_humidity     = 0.2_wp       ! very dry air -> sublimation
        forc%has_relative_humidity = .TRUE.

        precip = 0.0_wp_acc
        do istep = 1, 500
            precip = precip + real(forc%snowfall_rate*forc%dt_days*c%seconds_per_day,wp_acc)
            call bessi_column_step(bsi,1,forc,c)
        end do

        residual = closure_lhs(bsi,1) - precip

        write(*,"(a,g16.8)") "         precip            = ", precip
        write(*,"(a,g16.8)") "         sublimation       = ", bsi%now%sublimation(1)
        write(*,"(a,g16.8)") "         closure residual  = ", residual

        call check("residual is non-negative (over-reported sublimation only)", &
                   residual .ge. -1.0e-6_wp_acc*max(precip,1.0_wp_acc), nfail)
        call check("residual is bounded by the reported sublimation", &
                   residual .le. bsi%now%sublimation(1) &
                                 + 1.0e-6_wp_acc*max(precip,1.0_wp_acc), nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine probe_defect_1

    ! =====================================================================
    ! Check helpers (style follows tests/test_column_utils.f90)
    ! =====================================================================

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

        tol = max(abs(expected)*8.0_wp*epsilon(1.0_wp),tiny(1.0_wp))

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

        tol = max(abs(expected)*1.0e-6_wp_acc,1.0e-30_wp_acc)

        if (abs(value-expected) .le. tol) then
            write(*,"(a,a,a,g16.8)") "  ok   : ", trim(label), " = ", value
        else
            write(*,"(a,a,a,g16.8,a,g16.8)") "  FAIL : ", trim(label), &
                                             " = ", value, " expected ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_acc

    subroutine check_close(label,value,expected,reltol,nfail)
        ! Relative comparison in dp, reporting the relative residual so a
        ! marginal failure is diagnosable rather than just red.

        implicit none

        character(len=*), intent(IN)    :: label
        real(wp_acc),     intent(IN)    :: value
        real(wp_acc),     intent(IN)    :: expected
        real(wp_acc),     intent(IN)    :: reltol
        integer,          intent(INOUT) :: nfail

        ! Local variables
        real(wp_acc) :: scale, rel

        scale = max(abs(expected),1.0e-30_wp_acc)
        rel   = abs(value-expected)/scale

        if (rel .le. reltol) then
            write(*,"(a,a,a,es10.3)") "  ok   : ", trim(label), "  rel.resid = ", rel
        else
            write(*,"(a,a,a,es10.3,a,es10.3)") "  FAIL : ", trim(label), &
                          "  rel.resid = ", rel, " > ", reltol
            write(*,"(a,g20.12,a,g20.12)") "         value = ", value, &
                                           "  expected = ", expected
            nfail = nfail + 1
        end if

        return

    end subroutine check_close

end program test_bessi
