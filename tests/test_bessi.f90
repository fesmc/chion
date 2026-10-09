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
    ! Two upstream defects would break this identity; chion fixes the first
    ! and the test avoids the second rather than hiding it:
    !
    !   * defect 11 -- rain falling on a column with mass(1) <= 0 was silently
    !     dropped upstream. chion now routes it to runoff (D29), so every
    !     kilogram offered is accepted and no rain is withheld.
    !   * defect  1 -- apply_snow_surface_vapor_mass_flux returned the
    !     unclipped vapor_mass while applying a clipped one, so when
    !     sublimation demand exceeded the surface layer the diagnostic
    !     overstated the mass removed. Fixed upstream in Chion.jl 03bb445 and
    !     ported (C1): the diagnostic is the change actually applied. Test 1b
    !     asserts the identity with humidity on where nothing is clipped,
    !     test 1c where the surface layer is exhausted every step.
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
    use phys_constants, only : sec_day

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
    call test_mass_closure_humid(nfail)
    call test_mass_closure_exhausted(nfail)
    call test_cold_dry_column(nfail)
    call test_bare_and_recover(nfail)
    call test_capacity(nfail)
    call test_bare_ice_skips_water(nfail)
    call test_bare_rain_routed_once(nfail)
    call test_fresh_snow_split_on_bare(nfail)
    call test_diurnal(nfail)
    call test_scheme_matrix(nfail)
    call test_substrate_cold_content(nfail)
    call test_substrate_column(nfail)
    call test_fine_layers_column(nfail)
    call test_fine_layers_firn(nfail)

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

    subroutine run_annual_cycle(bsi,c,nyears,icol,precip,h_ice)
        ! Drive one column through nyears of the synthetic annual cycle,
        ! accumulating the precipitation mass offered. All of it is accepted:
        ! rain with no layer to hold it runs off (D29). h_ice is the host ice
        ! thickness (default 0, land: no ice substrate).

        implicit none

        type(bessi_class),       intent(INOUT) :: bsi
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(IN)    :: nyears
        integer,                 intent(IN)    :: icol
        real(wp_acc),            intent(INOUT) :: precip
        real(wp), optional,      intent(IN)    :: h_ice

        ! Local variables
        integer  :: iyr, iday
        real(wp) :: dt_seconds
        type(chion_step_forcing_class) :: forc

        do iyr = 1, nyears
            do iday = 1, NDAY_YEAR

                call annual_forcing(iday,c,forc)
                if (present(h_ice)) forc%H_ice = h_ice

                dt_seconds = forc%dt_days*real(sec_day,wp)

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
    ! Test 1b -- mass closure with humidity on (parameterized vapour mass)
    ! =====================================================================

    subroutine test_mass_closure_humid(nfail)
        ! The vapour mass of the parameterized BESSI turbulence is the
        ! humidity-gradient mass flux E*dt (Chion.jl d0146e1), applied to the
        ! surface layer and accumulated in vapor_mass. On a column that stays
        ! below T0 and keeps a surface layer far heavier than a day's
        ! sublimation, nothing is clipped (defect 1 cannot fire) and the
        ! closure identity holds to round-off.
        !
        ! The tolerance is that round-off, not the 1e-6 of test 1: the pack
        ! never melts, so every step rounds the surface layer (up to
        ! mass_max) at least twice -- snowfall and vapour -- and the storage
        ! sum and accumulators add as many again. Bound: 4 roundings of
        ! half an ulp of mass_max per step, relative to the precipitation;
        ! 2.4e-5 at sp (measured 2.7e-6), 4.4e-14 at dp (measured 2.7e-14).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)              :: bsi
        type(chion_const_class)        :: c
        type(chion_step_forcing_class) :: forc
        integer      :: iyr, iday
        real(wp)     :: dt_seconds, doy
        real(wp_acc) :: precip, reltol
        integer      :: nstep

        write(*,"(a)") "--- 1b. mass closure with humidity on (cold, no clipping) ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        precip = 0.0_wp_acc
        do iyr = 1, 3
            do iday = 1, NDAY_YEAR
                doy = real(iday,wp)
                call neutral_forcing(forc)
                forc%air_temperature = 255.0_wp - 8.0_wp*cos(2.0_wp*PI_WP*(doy - 15.0_wp) &
                                                             /real(NDAY_YEAR,wp))
                forc%dt_days         = 1.0_wp
                forc%snowfall_rate   = 3.0e-5_wp
                forc%shortwave_down  = max(200.0_wp*cos(2.0_wp*PI_WP*(doy - 197.0_wp) &
                                                        /real(NDAY_YEAR,wp)),0.0_wp)
                forc%wind_speed      = 5.0_wp
                forc%relative_humidity     = 0.5_wp
                forc%has_relative_humidity = .TRUE.

                dt_seconds = forc%dt_days*real(sec_day,wp)
                precip = precip + real(forc%snowfall_rate*dt_seconds,wp_acc)

                call bessi_column_step(bsi,1,forc,c)
            end do
        end do

        write(*,"(a,g16.8)") "         precip accepted [kg m-2] = ", precip
        write(*,"(a,g16.8)") "         vapor_mass      [kg m-2] = ", bsi%now%vapor_mass(1)
        write(*,"(a,g16.8)") "         sublimation     [kg m-2] = ", bsi%now%sublimation(1)

        call check("vapour exchange is exercised (sublimation > 1 kg m-2)", &
                   bsi%now%sublimation(1) .gt. 1.0_wp_acc, nfail)
        call check("column never melted (solid exchange only)", &
                   bsi%now%melt(1) .eq. 0.0_wp_acc, nfail)
        nstep  = 3*NDAY_YEAR
        reltol = 4.0_wp_acc*real(nstep,wp_acc)*0.5_wp_acc &
                 *real(spacing(bsi%par%mass_max),wp_acc)/precip
        call check_close("closure with humidity on (round-off bound)", &
                         closure_lhs(bsi,1),precip,reltol,nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_mass_closure_humid

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
    ! Test 12 -- bare ice over a cold substrate: cold content delays melt
    ! =====================================================================

    subroutine test_substrate_cold_content(nfail)
        ! A bare column (n = 0) over a 5-layer substrate at 263 K, under a
        ! constant prescribed surface flux Q and nothing else (no emission, no
        ! turbulence, no precipitation), so the surface balance is exactly Q.
        ! Without a substrate the bare ice is held at T0 and melts Q*dt/Lm
        ! every step. With one (Chion.jl 03bb445) the ice must first be warmed:
        !
        !   * the first steps do not melt at all -- the cold ice takes
        !     everything (Q = 20 W m-2 warms the top ~0.3 m by a few K a day);
        !   * melt only starts after the top layer has warmed;
        !   * energy closes: Lm*melt + sum_k m_k ci (T_k - 263) == Q*t, the
        !     insulated base keeping every joule in the substrate;
        !   * so the melt deficit against the no-substrate column is exactly the
        !     cold content taken up, (melt_0 - melt)*Lm == sum m ci dT.
        ! Bare ice is impermeable and solid: runoff = melt, smb_ice = -melt.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        integer,  parameter :: NSTEP = 40
        real(wp), parameter :: Q = 20.0_wp, T_ICE = 263.0_wp

        type(bessi_class)       :: bsi, bsi0
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer      :: istep, k, first_melt
        real(wp)     :: t1_at_first_melt
        real(wp_acc) :: cold_taken, total_in

        write(*,"(a)") "--- 12. bare ice over a cold substrate: cold content delays melt ---"

        call chion_const_init(c)
        c%eps_snow = 0.0_wp

        call bessi_par_init(bsi%par)
        bsi%par%ice_substrate_layers = 5
        call bessi_par_validate(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)
        bsi%now%ice_temperature = T_ICE
        bsi%now%t_srf           = T_ICE

        call bessi_par_init(bsi0%par)
        call bessi_alloc(bsi0,1)
        call bessi_init_state(bsi0,c)

        call neutral_forcing(forc)
        forc%air_temperature = 275.0_wp
        forc%has_q_sw_net    = .TRUE.
        forc%q_sw_net        = Q
        forc%has_q_lw_down   = .TRUE.
        forc%has_q_sh        = .TRUE.
        forc%has_q_lh        = .TRUE.
        forc%H_ice           = 1000.0_wp

        first_melt       = 0
        t1_at_first_melt = 0.0_wp

        do istep = 1, NSTEP
            call bessi_column_step(bsi, 1,forc,c)
            call bessi_column_step(bsi0,1,forc,c)
            if (first_melt .eq. 0) then
                if (bsi%now%melt(1) .gt. 0.0_wp_acc) then
                    first_melt       = istep
                    t1_at_first_melt = bsi%now%ice_temperature(1,1)
                end if
            end if
        end do

        cold_taken = 0.0_wp_acc
        do k = 1, 5
            cold_taken = cold_taken + real(c%rho_i*bsi%par%ice_substrate_top_thickness &
                                           *2.0_wp**(k-1),wp_acc)*real(c%ci,wp_acc) &
                                    *(real(bsi%now%ice_temperature(k,1),wp_acc) - real(T_ICE,wp_acc))
        end do
        total_in = real(Q,wp_acc)*real(NSTEP,wp_acc)*real(sec_day,wp_acc)

        write(*,"(a,i0)")    "         first melting step         = ", first_melt
        write(*,"(a,g14.6)") "         top ice T at first melt    = ", t1_at_first_melt
        write(*,"(a,g16.8)") "         melt, substrate  [kg m-2]  = ", bsi%now%melt(1)
        write(*,"(a,g16.8)") "         melt, none       [kg m-2]  = ", bsi0%now%melt(1)
        write(*,"(a,g16.8)") "         cold content taken [J m-2] = ", cold_taken

        call check("still bare (no layers)", bsi%now%n_lay(1) .eq. 0, nfail)
        call check("no substrate: melts from the first step", bsi0%now%melt(1) .gt. 0.0_wp_acc, nfail)
        call check("substrate: no melt in the first step (cold ice takes Q)", first_melt .gt. 1, nfail)
        call check("substrate: top ice warmed to near T0 before melting", &
                   t1_at_first_melt .gt. T_ICE + 5.0_wp, nfail)
        call check("substrate: melts eventually", first_melt .gt. 0 .and. first_melt .le. NSTEP, nfail)
        call check("substrate stays at or below T0", &
                   all(bsi%now%ice_temperature(:,1) .le. c%T0), nfail)
        call check_close("no substrate: melt == Q t/Lm", bsi0%now%melt(1), &
                         total_in/real(c%Lm,wp_acc), 1.0e-6_wp_acc, nfail)
        call check_close("substrate: Lm*melt + cold content taken == Q t", &
                         real(c%Lm,wp_acc)*bsi%now%melt(1) + cold_taken, total_in, &
                         1.0e-5_wp_acc, nfail)
        call check_close("melt deficit == cold content taken", &
                         real(c%Lm,wp_acc)*(bsi0%now%melt(1) - bsi%now%melt(1)), cold_taken, &
                         1.0e-3_wp_acc, nfail)
        call check_close("substrate: runoff == melt", bsi%now%runoff(1), bsi%now%melt(1), &
                         1.0e-12_wp_acc, nfail)
        call check_close("substrate: smb_ice == -melt", bsi%now%smb_ice(1), -bsi%now%melt(1), &
                         1.0e-12_wp_acc, nfail)

        call bessi_dealloc(bsi)
        call bessi_dealloc(bsi0)

        write(*,*)

        return

    end subroutine test_substrate_cold_content

    ! =====================================================================
    ! Test 13 -- the substrate in a full column; land; reset
    ! =====================================================================

    subroutine test_substrate_column(nfail)
        ! (a) Mass closure over a 5-yr annual cycle with the substrate on: the
        !     substrate carries no mass, so the identity is unchanged, while
        !     the run itself must differ from the no-substrate one.
        ! (b) A land column (H_ice = 0) with the substrate configured is
        !     bit-identical to a no-substrate run, and its substrate untouched
        !     (D34).
        ! (c) A column reset restores the substrate to temperature_init (D34).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi, bsi0
        type(chion_const_class) :: c
        real(wp_acc)            :: precip, precip0
        real(wp)                :: tmin

        write(*,"(a)") "--- 13. ice substrate in a full column; land; reset ---"

        call chion_const_init(c)

        call bessi_par_init(bsi%par)
        bsi%par%ice_substrate_layers = 5
        call bessi_par_validate(bsi%par)
        call bessi_alloc(bsi,2)
        call bessi_init_state(bsi,c)

        call bessi_par_init(bsi0%par)
        call bessi_alloc(bsi0,1)
        call bessi_init_state(bsi0,c)

        ! (a) column 1 on ice
        precip = 0.0_wp_acc
        call run_annual_cycle(bsi,c,5,1,precip,h_ice=1000.0_wp)
        precip0 = 0.0_wp_acc
        call run_annual_cycle(bsi0,c,5,1,precip0)

        tmin = minval(bsi%now%ice_temperature(:,1))
        write(*,"(a,g16.8,a,g16.8)") "         melt with / without substrate = ", &
                                     bsi%now%melt(1), " / ", bsi0%now%melt(1)
        write(*,"(a,g14.6)")         "         coldest substrate layer [K]  = ", tmin

        call check("substrate cooled below T0 over the cycle", tmin .lt. c%T0 - 1.0_wp, nfail)
        call check("substrate changes the melt", bsi%now%melt(1) .ne. bsi0%now%melt(1), nfail)
        call check_close("closure with the substrate on", closure_lhs(bsi,1),precip, &
                         1.0e-6_wp_acc,nfail)

        ! (b) column 2 on land, same configuration
        precip = 0.0_wp_acc
        call run_annual_cycle(bsi,c,5,2,precip,h_ice=0.0_wp)

        call check("land: identical to no substrate (n_lay, mass, temperature)", &
                   bsi%now%n_lay(2) .eq. bsi0%now%n_lay(1) .and. &
                   all(bsi%now%mass(:,2) .eq. bsi0%now%mass(:,1)) .and. &
                   all(bsi%now%temperature(:,2) .eq. bsi0%now%temperature(:,1)), nfail)
        call check("land: identical to no substrate (melt, runoff, smb_ice, t_srf)", &
                   bsi%now%melt(2) .eq. bsi0%now%melt(1) .and. &
                   bsi%now%runoff(2) .eq. bsi0%now%runoff(1) .and. &
                   bsi%now%smb_ice(2) .eq. bsi0%now%smb_ice(1) .and. &
                   bsi%now%t_srf(2) .eq. bsi0%now%t_srf(1), nfail)
        call check("land: substrate untouched", &
                   all(bsi%now%ice_temperature(:,2) .eq. bsi%par%temperature_init), nfail)

        ! (c) reset
        call bessi_reset_columns(bsi,c,[1])
        call check("reset: substrate back to temperature_init", &
                   all(bsi%now%ice_temperature(:,1) .eq. bsi%par%temperature_init), nfail)

        call bessi_dealloc(bsi)
        call bessi_dealloc(bsi0)

        write(*,*)

        return

    end subroutine test_substrate_column

    subroutine test_fine_layers_column(nfail)
        ! Fine near-surface layers (Chion.jl 03bb445, plan C4) over a 5-yr
        ! annual cycle that melts the column bare every summer: mass closure
        ! (the remesh moves mass, never creates or exports it), the limited
        ! layers back at their target thickness after every step that ends
        ! with snow and a layer below them (both remesh halves: snowfall caps
        ! down, melt fills up), and a run different from the default column.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi, bsi0
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        real(wp_acc) :: precip, precip0
        real(wp)     :: h(NEAR_SURFACE_LAYERS), dz
        integer      :: iyr, iday, k, n, n_max, n_checked, n_off
        logical      :: off

        write(*,"(a)") "--- 14. fine near-surface layers in a full column ---"

        h = [0.02_wp, 0.05_wp, 0.10_wp, 0.30_wp]

        call chion_const_init(c)

        call bessi_par_init(bsi%par)
        bsi%par%near_surface_layer_max_thicknesses = h
        call bessi_par_validate(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        call bessi_par_init(bsi0%par)
        call bessi_alloc(bsi0,1)
        call bessi_init_state(bsi0,c)

        precip    = 0.0_wp_acc
        n_max     = 0
        n_checked = 0
        n_off     = 0

        do iyr = 1, 5
            do iday = 1, NDAY_YEAR

                call annual_forcing(iday,c,forc)
                precip = precip + real(forc%snowfall_rate*forc%dt_days*real(sec_day,wp),wp_acc) &
                                + real(forc%rainfall_rate*forc%dt_days*real(sec_day,wp),wp_acc)
                call bessi_column_step(bsi,1,forc,c)

                n     = bsi%now%n_lay(1)
                n_max = max(n_max,n)
                if (.not. surface_has_snow(bsi%now%mass(:,1),n)) cycle

                off = .FALSE.
                do k = 1, min(NEAR_SURFACE_LAYERS,n-1)
                    dz  = bsi%now%mass(k,1)/bsi%now%density(k,1)
                    off = off .or. abs(dz - h(k)) .gt. 1.0e-5_wp*h(k)
                end do
                if (n .gt. 1) n_checked = n_checked + 1
                if (off) n_off = n_off + 1

            end do
        end do

        precip0 = 0.0_wp_acc
        call run_annual_cycle(bsi0,c,5,1,precip0)

        write(*,"(a,i0,a,i0)")       "         steps checked / off target   = ", n_checked, " / ", n_off
        write(*,"(a,i0)")            "         max layer count              = ", n_max
        write(*,"(a,g16.8,a,g16.8)") "         melt with / without          = ", &
                                     bsi%now%melt(1), " / ", bsi0%now%melt(1)

        call check("fine layers: geometry checked on many steps", n_checked .gt. 300, nfail)
        call check("fine layers: limited layers at target after every step", n_off .eq. 0, nfail)
        call check("fine layers: the cycle still melts the column bare", &
                   bsi%now%melt(1) .gt. precip, nfail)
        call check("fine layers change the result", bsi%now%melt(1) .ne. bsi0%now%melt(1), nfail)
        call check_close("closure with fine layers", closure_lhs(bsi,1),precip, &
                         1.0e-6_wp_acc,nfail)

        call bessi_dealloc(bsi)
        call bessi_dealloc(bsi0)

        write(*,*)

        return

    end subroutine test_fine_layers_column

    subroutine test_fine_layers_firn(nfail)
        ! A firn column under sustained cold accumulation with fine
        ! near-surface layers, 4 years at 8.6 kg m-2 d-1: Ntot capacity,
        ! bottom merges and the depth cap all act. chion (C4b, D32) splits
        ! and merges the first layer below the fine ones by mass, so the
        ! column keeps several layers there, each interior one (5 .. n-1;
        ! the bottom one collects bottom merges, as in BESSI without fine
        ! layers) within [mass_min, mass_max] after every step. Chion.jl,
        ! and legacy_chion builds, keep everything below the fine layers in
        ! layer 5. Mass closure in both.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        real(wp_acc) :: precip
        integer      :: istep, k, n, n_max, n_min_late, n_out
        real(wp)     :: m5_max

        write(*,"(a)") "--- 15. fine layers over a firn column under sustained accumulation ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        bsi%par%near_surface_layer_max_thicknesses = [0.02_wp, 0.05_wp, 0.10_wp, 0.30_wp]
        call bessi_par_validate(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)

        call neutral_forcing(forc)
        forc%air_temperature = 250.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 1.0e-4_wp          ! 8.64 kg m-2 d-1
        forc%shortwave_down  = 50.0_wp
        forc%wind_speed      = 4.0_wp

        precip     = 0.0_wp_acc
        n_max      = 0
        n_min_late = huge(1)
        n_out      = 0
        m5_max     = 0.0_wp

        do istep = 1, 4*NDAY_YEAR

            precip = precip + real(forc%snowfall_rate*forc%dt_days*real(sec_day,wp),wp_acc)
            call bessi_column_step(bsi,1,forc,c)

            n     = bsi%now%n_lay(1)
            n_max = max(n_max,n)
            if (istep .gt. NDAY_YEAR) n_min_late = min(n_min_late,n)
            if (n .ge. 5) m5_max = max(m5_max,bsi%now%mass(5,1))

            do k = 5, n-1
                if (bsi%now%mass(k,1) .lt. bsi%par%mass_min .or. &
                    bsi%now%mass(k,1) .gt. bsi%par%mass_max) n_out = n_out + 1
            end do

        end do

        write(*,"(a,i0,a,i0)")  "         layer count: max / min after year 1 = ", n_max, " / ", n_min_late
        write(*,"(a,g14.6)")    "         largest layer-5 mass                = ", m5_max
        write(*,"(a,i0)")       "         interior sub-fine layers out of bounds = ", n_out
        write(*,"(a,g16.8)")    "         mass_base                           = ", bsi%now%mass_base(1)

        if (NEAR_SURFACE_SPLIT_MERGE_BELOW) then
            call check("C4b: several layers below the fine ones after year 1 (n >= 8)", &
                       n_min_late .ge. 8, nfail)
            call check("C4b: column reaches Ntot and never exceeds it", &
                       n_max .eq. bsi%par%Ntot, nfail)
            call check("C4b: interior sub-fine layers within [mass_min, mass_max]", &
                       n_out .eq. 0, nfail)
        else
            call check("legacy: everything below the fine layers in layer 5 (n <= 5)", &
                       n_max .le. 5, nfail)
            call check("legacy: layer 5 outgrows mass_max", m5_max .gt. bsi%par%mass_max, nfail)
        end if
        call check("bottom export occurred (depth cap / bottom merge)", &
                   bsi%now%mass_base(1) .gt. 0.0_wp_acc, nfail)
        call check_close("closure, firn column with fine layers", closure_lhs(bsi,1),precip, &
                         1.0e-6_wp_acc,nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_fine_layers_firn

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
            precip = precip + real(forc%snowfall_rate*forc%dt_days*real(sec_day,wp),wp_acc)
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
                   thickness .lt. 2.0_wp_acc*real(BESSI_REFERENCE_SNOW_DEPTH_M,wp_acc), nfail)
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

    subroutine test_bare_rain_routed_once(nfail)
        ! Rain on a bare column runs off exactly once (D29). Two bare columns
        ! with only rain falling:
        !
        !   1. n = 0: nothing can hold the rain, so it must reach runoff.
        !   2. a surface layer below TOL_EMPTY_LAYER: bare by surface_has_snow,
        !      but mass(1) > 0, so the rain lands in mass_w(1) and must NOT
        !      also reach runoff. Chion.jl dev_nils 8fff530 counts it twice.
        !
        ! The bare-ice branch adds its own ice melt to runoff, so the rain's
        ! share is runoff - melt.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        real(wp_acc) :: rain_mass

        integer :: n_ice

        do n_ice = 0, 5, 5

        write(*,"(a,i0,a)") "--- 5b. rain on a bare column is routed exactly once (ice substrate ", &
                            n_ice, ") ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        bsi%par%mass_min = 1.0e-12_wp
        bsi%par%ice_substrate_layers = n_ice
        call bessi_par_validate(bsi%par)
        call bessi_alloc(bsi,2)
        call bessi_init_state(bsi,c)

        ! Column 1: no layers.
        bsi%now%n_lay(1) = 0

        ! Column 2: a sliver of a surface layer over a dry layer.
        bsi%now%n_lay(2)         = 2
        bsi%now%mass(1,2)        = 1.0e-11_wp
        bsi%now%mass(2,2)        = 100.0_wp
        bsi%now%mass_w(:,2)      = 0.0_wp
        bsi%now%density(1:2,2)   = 400.0_wp
        bsi%now%temperature(1:2,2) = 263.0_wp

        call neutral_forcing(forc)
        forc%air_temperature = 275.0_wp
        forc%dt_days         = 1.0_wp
        forc%rainfall_rate   = 2.0e-5_wp
        forc%shortwave_down  = 200.0_wp
        forc%wind_speed      = 2.0_wp
        forc%H_ice           = 1000.0_wp

        rain_mass = real(forc%rainfall_rate,wp_acc)*real(forc%dt_days,wp_acc) &
                   *real(sec_day,wp_acc)

        call bessi_column_step(bsi,1,forc,c)
        call bessi_column_step(bsi,2,forc,c)

        call check("n = 0: still bare", bsi%now%n_lay(1) .eq. 0, nfail)
        call check_close("n = 0: rain reaches runoff once", &
                         bsi%now%runoff(1) - bsi%now%melt(1), &
                         rain_mass, 1.0e-6_wp_acc, nfail)

        call check("sliver: still bare", &
                   .not. surface_has_snow(bsi%now%mass(:,2),bsi%now%n_lay(2)), nfail)
        call check_close("sliver: rain held in mass_w(1)", &
                         real(bsi%now%mass_w(1,2),wp_acc), rain_mass, 1.0e-6_wp_acc, nfail)
        call check("sliver: rain not also in runoff", &
                   abs(bsi%now%runoff(2) - bsi%now%melt(2)) .le. 1.0e-9_wp_acc, nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        end do

        return

    end subroutine test_bare_rain_routed_once

    ! =====================================================================
    ! Test 5c -- fresh snow on a bare column, split in the same step
    ! =====================================================================

    subroutine test_fresh_snow_split_on_bare(nfail)
        ! Chion.jl 03bb445: snow falling on a bare column takes the air
        ! temperature in EVERY new layer. A snowfall above mass_max splits in
        ! accumulation; before the fix the split copied the reset slot
        ! temperature (temperature_init = 273 K) into layer 2.
        !
        ! Cold air, no sunlight: nothing can warm the new layers above T_a
        ! within the step beyond the small conductive/longwave exchange.

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer :: n

        write(*,"(a)") "--- 5c. fresh snow split on a bare column takes T_a in every layer ---"

        call chion_const_init(c)
        call bessi_par_init(bsi%par)
        call bessi_alloc(bsi,1)
        call bessi_init_state(bsi,c)
        bsi%now%n_lay(1) = 0

        call neutral_forcing(forc)
        forc%air_temperature = 250.0_wp
        forc%dt_days         = 1.0_wp
        forc%snowfall_rate   = 1.5_wp*bsi%par%mass_max/real(sec_day,wp)
        forc%wind_speed      = 2.0_wp

        call bessi_column_step(bsi,1,forc,c)

        n = bsi%now%n_lay(1)
        call check("snowfall above mass_max split into >= 2 layers", n .ge. 2, nfail)
        call check("every new layer is near the air temperature, none at temperature_init", &
                   all(abs(bsi%now%temperature(1:n,1) - forc%air_temperature) .lt. 2.0_wp), nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_fresh_snow_split_on_bare

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

        ! --- 6c: an unsplit day keeps its forcing (D39) --------------------
        ! Polar night (80 N, solar longitude 270) with shortwave in the
        ! forcing: the criterion does not split the day. chion steps it with
        ! the forcing as given, i.e. exactly as with substepping off;
        ! legacy_chion (Chion.jl) averages it over [-pi, pi], which zeroes the
        ! shortwave, i.e. exactly as substepping off with no shortwave.
        call bessi_par_init(bsi_off%par)
        call bessi_alloc(bsi_off,1)
        call bessi_init_state(bsi_off,c)

        call bessi_par_init(bsi_on%par)
        bsi_on%par%diurnal_shortwave_substeps     = .TRUE.
        bsi_on%par%diurnal_shortwave_max_substeps = 8
        bsi_on%par%diurnal_temperature_cycle      = .TRUE.
        bsi_on%par%diurnal_temperature_amplitude  = 1.0_wp
        call bessi_par_validate(bsi_on%par)
        call bessi_alloc(bsi_on,1)
        call bessi_init_state(bsi_on,c)

        call neutral_forcing(forc)
        forc%air_temperature = 260.0_wp
        forc%snowfall_rate   = 3.0e-4_wp
        forc%wind_speed      = 2.0_wp
        do istep = 1, 30
            call bessi_column_step(bsi_off,1,forc,c)
            call bessi_column_step(bsi_on,1,forc,c)
        end do

        forc%air_temperature     = 266.0_wp
        forc%snowfall_rate       = 0.0_wp
        forc%shortwave_down      = 100.0_wp
        forc%latitude_deg        = 80.0_wp
        forc%solar_longitude_deg = 270.0_wp
        forc%day_of_year         = 356.0_wp
        do istep = 1, 10
            call bessi_column_step(bsi_on,1,forc,c)
            if (DIURNAL_SINGLE_INTERVAL_AVERAGED) forc%shortwave_down = 0.0_wp
            call bessi_column_step(bsi_off,1,forc,c)
            forc%shortwave_down = 100.0_wp
        end do

        if (DIURNAL_SINGLE_INTERVAL_AVERAGED) then
            call check("unsplit polar-night day: shortwave zeroed (legacy_chion, Chion.jl)", &
                       bsi_on%now%t_srf(1) .eq. bsi_off%now%t_srf(1) .and. &
                       all(bsi_on%now%temperature(:,1) .eq. bsi_off%now%temperature(:,1)), nfail)
        else
            call check("unsplit polar-night day: forcing kept (D39), = substeps off", &
                       bsi_on%now%t_srf(1) .eq. bsi_off%now%t_srf(1) .and. &
                       all(bsi_on%now%temperature(:,1) .eq. bsi_off%now%temperature(:,1)), nfail)
        end if

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

                dt_seconds = forc%dt_days*real(sec_day,wp)

                precip = precip + real(forc%snowfall_rate*dt_seconds,wp_acc) &
                                + real(forc%rainfall_rate*dt_seconds,wp_acc)

                call bessi_column_step(bsi,icol,forc,c)

            end do
        end do

        return

    end subroutine run_annual_cycle_albedo

    ! =====================================================================
    ! Test 1c -- mass closure when sublimation exhausts the surface layer
    ! =====================================================================

    subroutine test_mass_closure_exhausted(nfail)
        ! Very dry air and barely any snowfall: every step's sublimation
        ! demand exceeds the surface layer, so the solid exchange is clipped
        ! at zero mass. Upstream defect 1 reported the unclipped demand and
        ! left a closure residual equal to the over-reported sublimation;
        ! since Chion.jl 03bb445 (ported in C1) the diagnostic is the change
        ! actually applied, and the identity holds to round-off.
        !
        ! Tolerance: 4 roundings of half an ulp of the precipitation total
        ! per step, as in test 1b (the layer masses here are far smaller).

        implicit none

        integer, intent(INOUT) :: nfail

        ! Local variables
        type(bessi_class)       :: bsi
        type(chion_const_class) :: c
        type(chion_step_forcing_class) :: forc
        integer      :: istep
        real(wp_acc) :: precip, demand, reltol

        write(*,"(a)") "--- 1c. mass closure, sublimation exhausts the surface layer ---"

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
            precip = precip + real(forc%snowfall_rate*forc%dt_days*real(sec_day,wp),wp_acc)
            call bessi_column_step(bsi,1,forc,c)
        end do

        ! The unclipped demand of one step on the initial T, for the coverage
        ! check below: sublimation must have been limited by the mass available.
        demand = real(forc%snowfall_rate*forc%dt_days*real(sec_day,wp),wp_acc)

        write(*,"(a,g16.8)") "         precip            = ", precip
        write(*,"(a,g16.8)") "         sublimation       = ", bsi%now%sublimation(1)
        write(*,"(a,g16.8)") "         closure residual  = ", closure_lhs(bsi,1) - precip

        call check("the column is stripped bare (sublimation clipped)", &
                   bsi%now%n_lay(1) .eq. 0 .or. bsi%now%mass(1,1) .lt. demand, nfail)
        call check("sublimation removed (nearly) all the snowfall", &
                   bsi%now%sublimation(1) .gt. 0.9_wp_acc*precip, nfail)
        reltol = 4.0_wp_acc*500.0_wp_acc*0.5_wp_acc*real(spacing(real(precip,wp)),wp_acc)/precip
        call check_close("closure with an exhausted surface layer (defect 1 fixed)", &
                         closure_lhs(bsi,1),precip,reltol,nfail)

        call bessi_dealloc(bsi)

        write(*,*)

        return

    end subroutine test_mass_closure_exhausted

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
