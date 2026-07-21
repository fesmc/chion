module snow_diagnostics
    ! Column and domain diagnostics.
    !
    ! Port of Chion.jl/src/diagnostics.jl:
    !   _column_summary            (diagnostics.jl:123-143)  -> column_summary
    !   _summarize_domain_state_kernel! (:165-211)           -> summarize_domain_state
    !   _summarize_year_state_kernel!   (:220-240)           -> summarize_year_state
    !   _state_dict / get_state / print_state (:17-90)       -> get_state / print_state
    !
    ! CALLING CONVENTION: the single-column routines take contiguous column
    ! slices, e.g. mass(:,icol), plus the active layer count n. See
    ! docs/porting_notes.md D8.
    !
    ! =====================================================================
    ! WARNING -- THREE MUTUALLY INCONSISTENT REDUCTIONS COEXIST
    !
    ! chion now contains three different "sum over a column" behaviours, all
    ! of which are present in Chion.jl and all of which are ported verbatim:
    !
    !   1. snow_column_utils::bulk_snow_density
    !      Sums BOTH the mass numerator and the thickness denominator only
    !      over layers passing (mass > 0 .and. density > TOL_TINY).
    !
    !   2. snow_diagnostics::column_summary  (this module)
    !      Sums the thickness denominator only over qualifying layers, but
    !      the mass numerator UNCONDITIONALLY. The two therefore disagree
    !      whenever a layer has mass but near-zero density.
    !
    !   3. snow_column_utils::total_snow_water_mass
    !      Clips negative layer values at zero; column_summary's wet_mass and
    !      liquid_water do NOT.
    !
    ! These are not accidents of the port. They gate/report different things
    ! and the acceptance test tests/test_diagnostics.f90 asserts the
    ! disagreement explicitly. A future "cleanup" that unifies them must
    ! fail that test. See docs/PLAN.md section 5, item 1.
    ! =====================================================================
    !
    ! LAYERS BEYOND n ARE NEVER TOUCHED. Stale values above the active layer
    ! count are ignored, not zeroed, exactly as in Julia's 1:n loop.
    !
    ! PRECISION: every sum is accumulated in real(wp_acc). These are
    ! reductions over up to Ntot layers feeding a division; see
    ! docs/PLAN.md section 3.1.

    use chion_defs, only : wp, wp_acc, TOL_TINY, io_unit_err, chion_const_class

    implicit none

    private

    ! Maximum number of layers a single-column snapshot can hold. Snapshots
    ! are a diagnostic convenience, not part of the kernel path, so a fixed
    ! bound keeps the type assignable and free of allocation.
    integer, parameter, public :: SNAPSHOT_MAX_LAYERS = 64

    type snow_column_snapshot_class
        ! Fortran analogue of Chion.jl's _state_dict (diagnostics.jl:17-50).
        ! Profile entries 1:n are valid; entries above n are set to zero.

        integer  :: n                            ! active layer count

        real(wp) :: mass(SNAPSHOT_MAX_LAYERS)     ! [kg m-2] solid, per layer
        real(wp) :: mass_w(SNAPSHOT_MAX_LAYERS)   ! [kg m-2] liquid, per layer
        real(wp) :: density(SNAPSHOT_MAX_LAYERS)  ! [kg m-3] per layer
        real(wp) :: thickness(SNAPSHOT_MAX_LAYERS)! [m] per layer

        real(wp) :: total_mass                    ! [kg m-2] unconditional sum
        real(wp) :: total_liquid_water            ! [kg m-2] unconditional sum
        real(wp) :: total_wet_mass                ! [kg m-2] solid + liquid
        real(wp) :: total_thickness               ! [m] guarded sum
        real(wp) :: bulk_density                  ! [kg m-3] see column_summary

        real(wp) :: surface_temperature           ! [K] T(1), or T0 when n == 0
        real(wp) :: albedo                        ! [1]
        real(wp) :: smb_ice                       ! [kg m-2] cumulative
    end type snow_column_snapshot_class

    public :: snow_column_snapshot_class

    public :: column_summary
    public :: summarize_domain_state
    public :: summarize_year_state
    public :: get_state
    public :: print_state

contains

    pure subroutine column_summary(mass,mass_w,density,n, &
                                   thickness,wet_mass,bulk_density,liquid_water)
        ! Chion.jl/src/diagnostics.jl:123-143, argument order following the
        ! Julia signature _column_summary(N, mass, mass_w, density, idx).
        !
        ! The four inclusion criteria are deliberately NOT the same:
        !
        !   thickness    sums mass(k)/density(k) ONLY where
        !                mass(k) > 0 .and. density(k) > TOL_TINY.
        !   wet_mass     sums mass(k) + mass_w(k) UNCONDITIONALLY.
        !                Negative layer values are NOT clipped -- unlike
        !                snow_column_utils::total_snow_water_mass, which does
        !                clip. Both behaviours are ported; see the module head.
        !   liquid_water sums mass_w(k) unconditionally, also unclipped.
        !   bulk_density = solid_mass_total / thickness_total, where
        !                solid_mass_total is the UNCONDITIONAL mass sum --
        !                it therefore includes layers that were excluded from
        !                the thickness sum. This is what makes it disagree
        !                with snow_column_utils::bulk_snow_density, whose
        !                numerator is restricted to the same qualifying
        !                layers as its denominator. Guarded to 0 when
        !                thickness_total <= TOL_TINY.
        !
        ! Note bulk_density has NO "total_mass <= 0" guard, again unlike
        ! bulk_snow_density. A column whose masses sum to a negative value
        ! but whose qualifying layers give a positive thickness reports a
        ! negative bulk density. Preserved as-is.
        !
        ! Layers k > n are never read.

        implicit none

        real(wp), intent(IN)  :: mass(:)        ! (Ntot) [kg m-2] solid
        real(wp), intent(IN)  :: mass_w(:)      ! (Ntot) [kg m-2] liquid
        real(wp), intent(IN)  :: density(:)     ! (Ntot) [kg m-3]
        integer,  intent(IN)  :: n              ! active layer count
        real(wp), intent(OUT) :: thickness      ! [m]
        real(wp), intent(OUT) :: wet_mass       ! [kg m-2]
        real(wp), intent(OUT) :: bulk_density   ! [kg m-3]
        real(wp), intent(OUT) :: liquid_water   ! [kg m-2]

        ! Local variables
        integer      :: k
        real(wp_acc) :: solid, liquid, rho
        real(wp_acc) :: thickness_acc, wet_mass_acc, solid_mass_acc, liquid_acc

        thickness_acc  = 0.0_wp_acc
        wet_mass_acc   = 0.0_wp_acc
        solid_mass_acc = 0.0_wp_acc
        liquid_acc     = 0.0_wp_acc

        do k = 1, n

            solid  = real(mass(k),   wp_acc)
            liquid = real(mass_w(k), wp_acc)
            rho    = real(density(k),wp_acc)

            solid_mass_acc = solid_mass_acc + solid
            wet_mass_acc   = wet_mass_acc   + solid + liquid
            liquid_acc     = liquid_acc     + liquid

            ! Julia: `if solid > zero(solid) && rho > EPS_TINY`. Both halves
            ! are cheap and side-effect free, so a Fortran .and. is safe here
            ! -- no short-circuit dependence (docs/porting_notes.md D10).
            if (solid .gt. 0.0_wp_acc .and. rho .gt. TOL_TINY) then
                thickness_acc = thickness_acc + solid/rho
            end if

        end do

        if (thickness_acc .gt. TOL_TINY) then
            bulk_density = real(solid_mass_acc/thickness_acc,wp)
        else
            bulk_density = 0.0_wp
        end if

        thickness    = real(thickness_acc,wp)
        wet_mass     = real(wet_mass_acc, wp)
        liquid_water = real(liquid_acc,   wp)

        return

    end subroutine column_summary

    subroutine summarize_domain_state(mass,mass_w,density,n_lay, &
                                      thickness,wet_mass,bulk_density,liquid_water)
        ! Chion.jl/src/diagnostics.jl:165-211, restricted to the four fields
        ! that are actually computed. The kernel's other nine outputs
        ! (base_mass, smb_ice, runoff, melt, refreezing, vapor_mass,
        ! sublimation, latent_heat_flux_sum, albedo) are verbatim copies of
        ! state arrays; in Fortran they are plain assignments and belong to
        ! the caller (WP11/WP14), not here.

        implicit none

        real(wp), intent(IN)  :: mass(:,:)      ! (Ntot,ncol) [kg m-2]
        real(wp), intent(IN)  :: mass_w(:,:)    ! (Ntot,ncol) [kg m-2]
        real(wp), intent(IN)  :: density(:,:)   ! (Ntot,ncol) [kg m-3]
        integer,  intent(IN)  :: n_lay(:)       ! (ncol)
        real(wp), intent(OUT) :: thickness(:)   ! (ncol) [m]
        real(wp), intent(OUT) :: wet_mass(:)    ! (ncol) [kg m-2]
        real(wp), intent(OUT) :: bulk_density(:)! (ncol) [kg m-3]
        real(wp), intent(OUT) :: liquid_water(:)! (ncol) [kg m-2]

        ! Local variables
        integer :: i, ncol

        ncol = size(n_lay)

        call check_domain_shapes("summarize_domain_state",mass,mass_w,density,ncol)

        if (size(thickness)    .ne. ncol .or. size(wet_mass)     .ne. ncol .or. &
            size(bulk_density) .ne. ncol .or. size(liquid_water) .ne. ncol) then
            write(io_unit_err,*) "summarize_domain_state:: Error: output arrays must have length ncol."
            write(io_unit_err,*) "ncol = ", ncol
            write(io_unit_err,*) "sizes = ", size(thickness), size(wet_mass), &
                                             size(bulk_density), size(liquid_water)
            stop "Program stopped."
        end if

        do i = 1, ncol
            call column_summary(mass(:,i),mass_w(:,i),density(:,i),n_lay(i), &
                                thickness(i),wet_mass(i),bulk_density(i),liquid_water(i))
        end do

        return

    end subroutine summarize_domain_state

    subroutine summarize_year_state(mass,mass_w,density,n_lay,mass_base, &
                                    thickness,wet_mass,bulk_density,base_mass)
        ! Chion.jl/src/diagnostics.jl:220-240. The equilibrium-year subset:
        ! the same column_summary call, discarding liquid_water, plus a
        ! pass-through of the cumulative basal mass.
        !
        ! mass_base is wp_acc (a cumulative accumulator, docs/PLAN.md
        ! section 3.1); the reported base_mass is wp, matching the rest of
        ! the summary and the eventual NetCDF output.

        implicit none

        real(wp),     intent(IN)  :: mass(:,:)      ! (Ntot,ncol) [kg m-2]
        real(wp),     intent(IN)  :: mass_w(:,:)    ! (Ntot,ncol) [kg m-2]
        real(wp),     intent(IN)  :: density(:,:)   ! (Ntot,ncol) [kg m-3]
        integer,      intent(IN)  :: n_lay(:)       ! (ncol)
        real(wp_acc), intent(IN)  :: mass_base(:)   ! (ncol) [kg m-2] cumulative
        real(wp),     intent(OUT) :: thickness(:)   ! (ncol) [m]
        real(wp),     intent(OUT) :: wet_mass(:)    ! (ncol) [kg m-2]
        real(wp),     intent(OUT) :: bulk_density(:)! (ncol) [kg m-3]
        real(wp),     intent(OUT) :: base_mass(:)   ! (ncol) [kg m-2]

        ! Local variables
        integer  :: i, ncol
        real(wp) :: liquid_water_discard

        ncol = size(n_lay)

        call check_domain_shapes("summarize_year_state",mass,mass_w,density,ncol)

        if (size(mass_base)    .ne. ncol .or. size(thickness) .ne. ncol .or. &
            size(wet_mass)     .ne. ncol .or. size(base_mass) .ne. ncol .or. &
            size(bulk_density) .ne. ncol) then
            write(io_unit_err,*) "summarize_year_state:: Error: array lengths must equal ncol."
            write(io_unit_err,*) "ncol = ", ncol
            stop "Program stopped."
        end if

        do i = 1, ncol
            call column_summary(mass(:,i),mass_w(:,i),density(:,i),n_lay(i), &
                                thickness(i),wet_mass(i),bulk_density(i),liquid_water_discard)
            base_mass(i) = real(mass_base(i),wp)
        end do

        return

    end subroutine summarize_year_state

    subroutine get_state(snap,mass,mass_w,density,temperature,n,albedo,smb_ice,c)
        ! Chion.jl/src/diagnostics.jl:17-71 (_state_dict / get_state), for a
        ! single column.
        !
        ! HAZARD IN THE JULIA ORIGINAL, deliberately NOT reproduced:
        ! _state_dict computes the per-layer thickness as the unguarded
        ! elementwise division
        !     thickness = active_solid_mass ./ active_density
        ! (diagnostics.jl:32) and then total_thickness = sum(thickness)
        ! (:45). An active layer with zero density therefore yields Inf, and
        ! a zero-mass zero-density layer yields NaN, both of which poison the
        ! total. The kernel path (_column_summary) guards the same division
        ! and is the safe one. chion uses the guarded form here, so
        ! total_thickness and bulk_density agree with column_summary and a
        ! -ffpe-trap=zero,invalid build does not trap on a diagnostic call.
        ! Reported upstream; see docs/porting_notes.md.
        !
        ! Note total_mass, total_liquid_water and total_wet_mass remain the
        ! unconditional (unclipped) sums, matching both _state_dict and
        ! column_summary.

        implicit none

        type(snow_column_snapshot_class), intent(OUT) :: snap
        real(wp),                intent(IN) :: mass(:)        ! (Ntot) [kg m-2]
        real(wp),                intent(IN) :: mass_w(:)      ! (Ntot) [kg m-2]
        real(wp),                intent(IN) :: density(:)     ! (Ntot) [kg m-3]
        real(wp),                intent(IN) :: temperature(:) ! (Ntot) [K]
        integer,                 intent(IN) :: n              ! active layers
        real(wp),                intent(IN) :: albedo         ! [1]
        real(wp_acc),            intent(IN) :: smb_ice        ! [kg m-2] cumulative
        type(chion_const_class), intent(IN) :: c

        ! Local variables
        integer      :: k
        real(wp)     :: thickness_total, wet_mass_total, bulk_density_now, liquid_total
        real(wp_acc) :: solid_acc, rho_acc, solid_mass_acc

        if (n .lt. 0 .or. n .gt. min(size(mass),SNAPSHOT_MAX_LAYERS)) then
            write(io_unit_err,*) "get_state:: Error: active layer count out of range."
            write(io_unit_err,*) "n, size(mass), SNAPSHOT_MAX_LAYERS = ", &
                                 n, size(mass), SNAPSHOT_MAX_LAYERS
            stop "Program stopped."
        end if

        snap%n = n

        snap%mass      = 0.0_wp
        snap%mass_w    = 0.0_wp
        snap%density   = 0.0_wp
        snap%thickness = 0.0_wp

        ! Profiles, plus the per-layer thickness guarded exactly as
        ! column_summary guards its sum, so profile and total are consistent.
        solid_mass_acc = 0.0_wp_acc

        do k = 1, n

            snap%mass(k)    = mass(k)
            snap%mass_w(k)  = mass_w(k)
            snap%density(k) = density(k)

            solid_acc = real(mass(k),   wp_acc)
            rho_acc   = real(density(k),wp_acc)

            solid_mass_acc = solid_mass_acc + solid_acc

            if (solid_acc .gt. 0.0_wp_acc .and. rho_acc .gt. TOL_TINY) then
                snap%thickness(k) = real(solid_acc/rho_acc,wp)
            else
                snap%thickness(k) = 0.0_wp
            end if

        end do

        call column_summary(mass,mass_w,density,n, &
                            thickness_total,wet_mass_total,bulk_density_now,liquid_total)

        snap%total_thickness    = thickness_total
        snap%total_wet_mass     = wet_mass_total
        snap%bulk_density       = bulk_density_now
        snap%total_liquid_water = liquid_total
        snap%total_mass         = real(solid_mass_acc,wp)

        if (n .eq. 0) then
            snap%surface_temperature = c%T0
        else
            snap%surface_temperature = temperature(1)
        end if

        snap%albedo  = albedo
        snap%smb_ice = real(smb_ice,wp)

        return

    end subroutine get_state

    subroutine print_state(snap,icol)
        ! Chion.jl/src/diagnostics.jl:79-90 (print_state), plus the layer
        ! profile, which Julia leaves to the caller to pull out of the dict.

        implicit none

        type(snow_column_snapshot_class), intent(IN) :: snap
        integer, optional,                intent(IN) :: icol

        ! Local variables
        integer :: k, idx

        idx = 1
        if (present(icol)) idx = icol

        write(*,"(a)") repeat("=",60)
        write(*,"(a)") "Snowpack Column State"
        write(*,"(a)") repeat("=",60)
        write(*,"(a,i0)")          "Column index:    ", idx
        write(*,"(a,i0)")          "Active layers:   ", snap%n
        write(*,"(a,f14.2,a)")     "Total mass:      ", snap%total_mass,      " kg m-2"
        write(*,"(a,f14.2,a)")     "Liquid water:    ", snap%total_liquid_water," kg m-2"
        write(*,"(a,f14.2,a)")     "Total wet mass:  ", snap%total_wet_mass,  " kg m-2"
        write(*,"(a,f14.3,a)")     "Total thickness: ", snap%total_thickness, " m"
        write(*,"(a,f14.3,a)")     "Bulk density:    ", snap%bulk_density,    " kg m-3"
        write(*,"(a,f14.3,a)")     "Surface temp:    ", snap%surface_temperature," K"
        write(*,"(a,f14.3)")       "Surface albedo:  ", snap%albedo
        write(*,"(a,f14.2,a)")     "smb_ice:         ", snap%smb_ice,         " kg m-2"

        if (snap%n .gt. 0) then
            write(*,"(a)") ""
            write(*,"(a4,4a14)") "k", "mass", "mass_w", "density", "thickness"
            do k = 1, snap%n
                write(*,"(i4,4g14.6)") k, snap%mass(k), snap%mass_w(k), &
                                          snap%density(k), snap%thickness(k)
            end do
        end if

        write(*,"(a)") ""

        return

    end subroutine print_state

    subroutine check_domain_shapes(routine,mass,mass_w,density,ncol)
        ! Shared shape validation for the two domain reductions.

        implicit none

        character(len=*), intent(IN) :: routine
        real(wp),         intent(IN) :: mass(:,:)
        real(wp),         intent(IN) :: mass_w(:,:)
        real(wp),         intent(IN) :: density(:,:)
        integer,          intent(IN) :: ncol

        if (size(mass,2) .ne. ncol .or. size(mass_w,2) .ne. ncol .or. &
            size(density,2) .ne. ncol) then
            write(io_unit_err,*) trim(routine)//":: Error: layer arrays must be (Ntot,ncol)."
            write(io_unit_err,*) "ncol = ", ncol
            write(io_unit_err,*) "shape(mass), shape(mass_w), shape(density) = ", &
                                 shape(mass), shape(mass_w), shape(density)
            stop "Program stopped."
        end if

        if (size(mass,1) .ne. size(mass_w,1) .or. size(mass,1) .ne. size(density,1)) then
            write(io_unit_err,*) trim(routine)//":: Error: layer arrays must share Ntot."
            write(io_unit_err,*) "Ntot = ", size(mass,1), size(mass_w,1), size(density,1)
            stop "Program stopped."
        end if

        return

    end subroutine check_domain_shapes

end module snow_diagnostics
