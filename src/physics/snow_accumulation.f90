module snow_accumulation
    ! Snowfall and rainfall accumulation, fresh-snow density, and the
    ! layer-structure enforcement that follows mass addition.
    !
    ! Port of Chion.jl/src/processes/accumulation.jl.
    !
    ! CALLING CONVENTION: contiguous column slices plus the active layer count
    ! n, per docs/porting_notes.md D8. n is intent(INOUT) because the split,
    ! merge and depth-cap routines change it.
    !
    ! DEPENDENCY ON WP4: everything in the "layer-structure enforcement" section
    ! is delegated to snow_layers. All five calls are confined to the bottom of
    ! apply_accumulation, so a change to those signatures touches nothing else
    ! in this module. Two of them differ from the Julia argument lists: WP4
    ! dropped Ntot from merge_surface_layer, and Ntot and dt_seconds from
    ! enforce_snow_depth_cap, because neither body uses them.
    !
    ! PRESERVED QUIRKS:
    !   * Rain alone never creates a snow layer, and rain is only added to
    !     mass_w(1) when mass(1) > 0 STRICTLY (not > TOL_TINY). A column that
    !     is bare and receives only rain routes nothing here: the rain is
    !     simply dropped by this routine, and the bare-ice branch of the column
    !     step handles it. Preserved as-is; flagged upstream.
    !   * The split loop's out-of-slots branch is asymmetric: with Ntot <= 2 it
    !     calls free_slot_for_surface_split, otherwise merge_bottom_layer. The
    !     Ntot <= 2 case cannot merge a bottom layer without destroying the only
    !     other layer, hence the special case.
    !   * The density mix is volume-weighted, not mass-weighted:
    !         rho_new = m_new / (m_prev/rho_prev + m_added/rho_fresh)
    !     so it conserves layer volume, not mean density.
    !   * The fresh-snow density clamp lower bound 50 kg m-3 is a bare literal
    !     upstream (accumulation.jl:18,25); named here.

    use chion_defs, only : wp, wp_acc, chion_const_class, &
                           CHION_FRESH_SNOW_DENSITY_CONSTANT, &
                           CHION_FRESH_SNOW_DENSITY_PARAMETERIZED, &
                           io_unit_err
    use snow_albedo, only : albedo_refresh_from_snowfall

    use snow_layers, only : split_surface_layer, merge_surface_layer, &
                            merge_bottom_layer, free_slot_for_surface_split, &
                            enforce_snow_depth_cap

    implicit none

    private

    ! Fresh-snow density clamp bounds. The upper bound is c%rho_i.
    real(wp), parameter, public :: FRESH_SNOW_DENSITY_MIN = 50.0_wp   ! [kg m-3]

    public :: fresh_snow_density
    public :: apply_accumulation

contains

    function fresh_snow_density(c,air_temperature,wind_speed) result(rho_fresh)
        ! Chion.jl/src/processes/accumulation.jl:12-26.
        !
        !   constant       : clamp(rho_s, 50, rho_i)
        !   parameterized  : clamp(a + b*(T - T0) + c*sqrt(max(V,0)), 50, rho_i)
        !
        ! Both schemes clamp, so the constant scheme is not simply "rho_s" --
        ! a badly configured rho_s is silently pulled into range.

        implicit none

        type(chion_const_class), intent(IN) :: c
        real(wp),                intent(IN) :: air_temperature   ! [K]
        real(wp),                intent(IN) :: wind_speed        ! [m s-1]
        real(wp) :: rho_fresh                                    ! [kg m-3]

        ! Local variables
        real(wp) :: wind

        select case(c%fresh_snow_density_scheme)

            case(CHION_FRESH_SNOW_DENSITY_CONSTANT)

                rho_fresh = c%rho_s

            case(CHION_FRESH_SNOW_DENSITY_PARAMETERIZED)

                wind = max(wind_speed,0.0_wp)

                rho_fresh = c%rho_s_a &
                          + c%rho_s_b*(air_temperature - c%T0) &
                          + c%rho_s_c*sqrt(wind)

            case DEFAULT

                ! Chion.jl uses "if constant ... else parameterized", so an
                ! unrecognized flag silently selects the parameterized form.
                ! chion refuses instead, mirroring docs/porting_notes.md D5.
                ! Unreachable in practice: the flag is validated at parse time
                ! by chion_fresh_snow_density_scheme_flag.
                write(io_unit_err,*) "fresh_snow_density:: Error: &
                                     &fresh_snow_density_scheme not recognized."
                write(io_unit_err,*) "fresh_snow_density_scheme should be one of: &
                                     &[CHION_FRESH_SNOW_DENSITY_CONSTANT, &
                                     &CHION_FRESH_SNOW_DENSITY_PARAMETERIZED]"
                write(io_unit_err,*) "fresh_snow_density_scheme = ", c%fresh_snow_density_scheme
                stop "Program stopped."

        end select

        rho_fresh = min(max(rho_fresh,FRESH_SNOW_DENSITY_MIN),c%rho_i)

        return

    end function fresh_snow_density

    subroutine apply_accumulation(mass,mass_w,density,temperature,n, &
                                  mass_base,smb_ice,runoff,t_srf,albedo, &
                                  c,Ntot,mass_max,mass_split,mass_min, &
                                  snowfall_rate,rainfall_rate,dt_seconds, &
                                  air_temperature,wind_speed)
        ! Chion.jl/src/processes/accumulation.jl:36-152
        ! (_apply_accumulation_resolved!). Step 1 of column_step_core!.
        !
        ! Flow:
        !   0. empty column: create a surface layer only if snowfall > 0,
        !      otherwise set albedo = alpha_ice and return.
        !   1. snowfall: add mass to layer 1, mix its density by volume, and
        !      brighten the albedo.
        !   2. rainfall: add to mass_w(1) if mass(1) > 0.
        !   3. split loop  while mass(1) > mass_max
        !   4. merge loop  while n > 1 and mass(1) < mass_min
        !   5. depth cap.
        !
        ! NOTE the albedo call here is the snowfall REFRESH, not the aging
        ! update. Aging happens later in the step, exactly once per call
        ! (docs/PLAN.md section 5, item 5).

        implicit none

        real(wp),                intent(INOUT) :: mass(:)         ! (Ntot) [kg m-2] solid
        real(wp),                intent(INOUT) :: mass_w(:)       ! (Ntot) [kg m-2] liquid
        real(wp),                intent(INOUT) :: density(:)      ! (Ntot) [kg m-3]
        real(wp),                intent(INOUT) :: temperature(:)  ! (Ntot) [K]
        integer,                 intent(INOUT) :: n               ! active layers
        real(wp_acc),            intent(INOUT) :: mass_base       ! [kg m-2] cumulative
        real(wp_acc),            intent(INOUT) :: smb_ice         ! [kg m-2] cumulative
        real(wp_acc),            intent(INOUT) :: runoff          ! [kg m-2] cumulative
        real(wp),                intent(INOUT) :: t_srf           ! [K] surface temperature
        real(wp),                intent(INOUT) :: albedo          ! [1]
        type(chion_const_class), intent(IN)    :: c
        integer,                 intent(IN)    :: Ntot            ! layer capacity
        real(wp),                intent(IN)    :: mass_max        ! [kg m-2] split trigger
        real(wp),                intent(IN)    :: mass_split      ! [kg m-2] mass left below
        real(wp),                intent(IN)    :: mass_min        ! [kg m-2] merge trigger
        real(wp),                intent(IN)    :: snowfall_rate   ! [kg m-2 s-1]
        real(wp),                intent(IN)    :: rainfall_rate   ! [kg m-2 s-1]
        real(wp),                intent(IN)    :: dt_seconds      ! [s]
        real(wp),                intent(IN)    :: air_temperature ! [K]
        real(wp),                intent(IN)    :: wind_speed      ! [m s-1]

        ! Local variables
        real(wp) :: m_prev, m_added, m_new
        real(wp) :: rho_fresh, rho_prev, rho_new

        ! --- Step 0: empty column -------------------------------------------
        if (n .eq. 0) then
            if (snowfall_rate .gt. 0.0_wp) then
                n = 1
            else
                albedo = c%alpha_ice
                return
            end if
        end if

        ! --- Step 1: snowfall -----------------------------------------------
        if (snowfall_rate .gt. 0.0_wp) then

            m_prev  = mass(1)
            m_added = snowfall_rate*dt_seconds
            m_new   = m_prev + m_added

            rho_fresh = fresh_snow_density(c,air_temperature,wind_speed)

            ! A newly created (or reset) layer has density 0; fall back to the
            ! fresh-snow density so the volume mix below is well posed.
            if (density(1) .gt. 0.0_wp) then
                rho_prev = density(1)
            else
                rho_prev = rho_fresh
            end if

            if (m_new .gt. 0.0_wp) then
                ! Volume-weighted mix: conserves m_prev/rho_prev + m_added/rho_fresh.
                rho_new    = m_new/(m_prev/rho_prev + m_added/rho_fresh)
                density(1) = rho_new
            end if

            mass(1) = m_new

            call albedo_refresh_from_snowfall(albedo,c,m_added)

        end if

        ! --- Step 2: rainfall -----------------------------------------------
        ! Strict mass(1) > 0, and only into the surface layer.
        if (mass(1) .gt. 0.0_wp .and. rainfall_rate .gt. 0.0_wp) then
            mass_w(1) = mass_w(1) + rainfall_rate*dt_seconds
        end if

        ! --- Step 3: split loop ---------------------------------------------
        do while (n .gt. 0 .and. mass(1) .gt. mass_max)

            if (n .eq. Ntot) then

                if (Ntot .le. 2) then
                    call free_slot_for_surface_split(mass,mass_w,density,temperature,n, &
                                                     mass_base,smb_ice,runoff,t_srf,albedo, &
                                                     Ntot,mass_max,c)
                else
                    call merge_bottom_layer(mass,mass_w,density,temperature,n, &
                                            mass_base,smb_ice,c)
                end if

                ! Re-test after freeing a slot: either the column emptied or the
                ! surface is already back under mass_max.
                if (n .eq. 0) exit
                if (mass(1) .le. mass_max) exit

            end if

            call split_surface_layer(mass,mass_w,density,temperature,n, &
                                     Ntot,mass_max,mass_split)

        end do

        ! --- Step 4: merge loop ---------------------------------------------
        do while (n .gt. 1 .and. mass(1) .lt. mass_min)
            call merge_surface_layer(mass,mass_w,density,temperature,n, &
                                     mass_split,mass_min,c)
        end do

        ! --- Step 5: depth cap ----------------------------------------------
        call enforce_snow_depth_cap(mass,mass_w,density,temperature,n, &
                                    mass_base,smb_ice,runoff,t_srf,albedo, &
                                    mass_split,c)

        return

    end subroutine apply_accumulation

end module snow_accumulation
