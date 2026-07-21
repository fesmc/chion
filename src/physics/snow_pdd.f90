module snow_pdd
    ! Bulk positive-degree-day surface mass balance.
    !
    ! Port of Chion.jl/src/processes/pdd.jl (branch main), cross-checked
    ! against smbpal/src/smb_pdd.f90 (calc_ablation_pdd, calc_temp_effective),
    ! which is the PDD scheme actually in production use in yelmox.
    !
    ! WARNING. Chion.jl's PDDModel is known NOT to be fully working, so this
    ! module is a port of a defective reference. The physics below reproduces
    ! Chion.jl deliberately, defect for defect, so that the two models can be
    ! compared field by field in WP16. Every known defect is catalogued, with
    ! its physical consequence and a recommended upstream fix, in
    ! docs/pdd_defects.md. Read that file before changing anything here.
    !
    ! This is a BULK model: no layers, no temperature profile, no density, no
    ! energy balance, no albedo. State is four per-column scalars.
    !
    ! Units. Degree-day factors are [kg m-2 K-1 day-1], PDD sums are [K day],
    ! all masses are [kg m-2] (equivalently mm w.e.).
    !
    ! CALLING CONVENTION follows docs/porting_notes.md D8: the column kernel
    ! takes the four state scalars for one column, not (array, index).

    use chion_defs, only : wp, wp_acc, io_unit_err, &
                           chion_const_class, chion_forcing_class, &
                           chion_step_forcing_class

    implicit none

    private

    ! === PDD flavour =========================================================
    !
    ! Chion.jl selects between these two implicitly, on the length of the
    ! timestep: _is_monthly_pdd_step (pdd.jl:13) fires when
    ! 27 <= dt_days <= 32. chion replaces that with an explicit parameter.
    ! A change of physics must never be triggered by a timestep length; see
    ! docs/pdd_defects.md D4. The default reproduces Julia for daily forcing.

    integer, parameter, public :: CHION_PDD_SIMPLE = 1   ! max(T-T0,0)*dt
    integer, parameter, public :: CHION_PDD_PISM   = 2   ! Calov-Greve integral

    ! === Numerical constants =================================================
    !
    ! Standard normal pdf normalisation, 1/sqrt(2*pi), and 1/sqrt(2) for the
    ! erfc argument. Both dp: the expectation integral is evaluated entirely in
    ! wp_acc, matching smbpal, which explicitly promotes to real(8) around erfc
    ! "to avoid potential underflow errors at very low temperatures"
    ! (smb_pdd.f90:91).

    real(wp_acc), parameter :: PDD_INV_SQRT_2PI = 0.398942280401432678_wp_acc
    real(wp_acc), parameter :: PDD_INV_SQRT_2   = 0.707106781186547524_wp_acc

    ! === Parameters ==========================================================

    type pdd_par_class
        ! Mirrors Chion.jl PDDModel (src/models.jl:89-95), plus the explicit
        ! method flag that replaces the implicit monthly trigger.

        real(wp) :: ddf_snow             ! [kg m-2 K-1 d-1] snow degree-day factor
        real(wp) :: ddf_ice              ! [kg m-2 K-1 d-1] ice degree-day factor
        real(wp) :: refreezing_fraction  ! [1] fraction of SNOW melt that refreezes
        real(wp) :: temperature_sigma    ! [K] std dev of unresolved T variability

        integer  :: pdd_method           ! CHION_PDD_SIMPLE | CHION_PDD_PISM
    end type pdd_par_class

    ! === State ===============================================================

    type pdd_state_class
        ! Chion.jl PDDState (src/state.jl:180-185): four (ncol) vectors.
        !
        ! snowpack_swe is a reservoir, so wp. The other three are cumulative
        ! and never reset, so they MUST be wp_acc -- see docs/PLAN.md 3.1.
        !
        ! NOTE snowpack_swe has no cap, no aging and no densification, so it
        ! grows without bound in any column that never melts (defect D3). If
        ! that defect is not fixed upstream, this field will eventually reach
        ! magnitudes where sp cannot resolve a daily snowfall increment, and it
        ! must then be promoted to wp_acc as well.

        integer :: ncol

        real(wp),     allocatable :: snowpack_swe(:)   ! [kg m-2] snow reservoir
        real(wp_acc), allocatable :: smb_ice(:)        ! [kg m-2] cumulative
        real(wp_acc), allocatable :: runoff(:)         ! [kg m-2] cumulative
        real(wp_acc), allocatable :: pdd_sum(:)        ! [K d]    cumulative
    end type pdd_state_class

    type pdd_class
        type(pdd_par_class)   :: par
        type(pdd_state_class) :: now
    end type pdd_class

    public :: pdd_par_class
    public :: pdd_state_class
    public :: pdd_class

    public :: pdd_par_init
    public :: pdd_par_validate
    public :: pdd_method_flag

    public :: pdd_alloc
    public :: pdd_dealloc
    public :: pdd_init_state
    public :: pdd_reset_column

    public :: pdd_normal_cdf
    public :: pdd_expected_positive_temperature
    public :: pdd_degree_days
    public :: pdd_step_mass

    public :: pdd_column_apply
    public :: pdd_column_step
    public :: pdd_step

contains

    ! =====================================================================
    ! Parameters
    ! =====================================================================

    subroutine pdd_par_init(par)
        ! Chion.jl PDDModel keyword defaults (src/models.jl:97-107).
        ! WP13 will override these from the namelist.

        implicit none

        type(pdd_par_class), intent(OUT) :: par

        par%ddf_snow            = 3.0_wp
        par%ddf_ice             = 8.0_wp
        par%refreezing_fraction = 0.6_wp
        par%temperature_sigma   = 5.0_wp

        ! Default reproduces Chion.jl for sub-monthly forcing, which is the
        ! only regime Chion.jl's daily path exercises.
        par%pdd_method          = CHION_PDD_SIMPLE

        call pdd_par_validate(par)

        return

    end subroutine pdd_par_init

    subroutine pdd_par_validate(par)
        ! The four validations from Chion.jl's PDDModel constructor
        ! (src/models.jl:104-107), promoted from Julia `error()` to a stop.

        implicit none

        type(pdd_par_class), intent(IN) :: par

        if (par%ddf_snow .le. 0.0_wp) then
            write(io_unit_err,*) "pdd_par_validate:: Error: ddf_snow must be positive."
            write(io_unit_err,*) "ddf_snow = ", par%ddf_snow
            stop "Program stopped."
        end if

        if (par%ddf_ice .lt. 0.0_wp) then
            write(io_unit_err,*) "pdd_par_validate:: Error: ddf_ice must be non-negative."
            write(io_unit_err,*) "ddf_ice = ", par%ddf_ice
            stop "Program stopped."
        end if

        if (par%refreezing_fraction .lt. 0.0_wp .or. par%refreezing_fraction .gt. 1.0_wp) then
            write(io_unit_err,*) "pdd_par_validate:: Error: refreezing_fraction must be in [0,1]."
            write(io_unit_err,*) "refreezing_fraction = ", par%refreezing_fraction
            stop "Program stopped."
        end if

        if (par%temperature_sigma .le. 0.0_wp) then
            write(io_unit_err,*) "pdd_par_validate:: Error: temperature_sigma must be positive."
            write(io_unit_err,*) "temperature_sigma = ", par%temperature_sigma
            stop "Program stopped."
        end if

        if (par%pdd_method .ne. CHION_PDD_SIMPLE .and. &
            par%pdd_method .ne. CHION_PDD_PISM) then
            write(io_unit_err,*) "pdd_par_validate:: Error: pdd_method not recognized."
            write(io_unit_err,*) "pdd_method should be one of: ['simple','pism']"
            write(io_unit_err,*) "pdd_method = ", par%pdd_method
            stop "Program stopped."
        end if

        return

    end subroutine pdd_par_validate

    function pdd_method_flag(name) result(flag)
        ! Map a namelist string onto a PDD flavour flag. Follows the pattern of
        ! chion_albedo_scheme_flag in chion_defs.

        implicit none

        character(len=*), intent(IN) :: name
        integer :: flag

        select case(trim(adjustl(name)))
            case("simple")
                flag = CHION_PDD_SIMPLE
            case("pism","calov_greve","calov-greve")
                flag = CHION_PDD_PISM
            case DEFAULT
                write(io_unit_err,*) "pdd_method_flag:: Error: pdd method not recognized."
                write(io_unit_err,*) "pdd_method should be one of: ['simple','pism']"
                write(io_unit_err,*) "pdd_method = ", trim(name)
                stop "Program stopped."
        end select

        return

    end function pdd_method_flag

    ! =====================================================================
    ! Allocation and state
    ! =====================================================================

    subroutine pdd_alloc(pdd,ncol)

        implicit none

        type(pdd_class), intent(INOUT) :: pdd
        integer,         intent(IN)    :: ncol

        call pdd_dealloc(pdd)

        if (ncol .le. 0) then
            write(io_unit_err,*) "pdd_alloc:: Error: ncol must be positive."
            write(io_unit_err,*) "ncol = ", ncol
            stop "Program stopped."
        end if

        pdd%now%ncol = ncol

        allocate(pdd%now%snowpack_swe(ncol))
        allocate(pdd%now%smb_ice(ncol))
        allocate(pdd%now%runoff(ncol))
        allocate(pdd%now%pdd_sum(ncol))

        call pdd_init_state(pdd)

        return

    end subroutine pdd_alloc

    subroutine pdd_dealloc(pdd)

        implicit none

        type(pdd_class), intent(INOUT) :: pdd

        if (allocated(pdd%now%snowpack_swe)) deallocate(pdd%now%snowpack_swe)
        if (allocated(pdd%now%smb_ice))      deallocate(pdd%now%smb_ice)
        if (allocated(pdd%now%runoff))       deallocate(pdd%now%runoff)
        if (allocated(pdd%now%pdd_sum))      deallocate(pdd%now%pdd_sum)

        pdd%now%ncol = 0

        return

    end subroutine pdd_dealloc

    subroutine pdd_init_state(pdd)
        ! Cold start. Chion.jl PDDState(model) zeroes all four fields
        ! (src/state.jl:187-195). Note there is no equivalent of BESSI's
        ! density_init / temperature_init: a PDD column starts bare.

        implicit none

        type(pdd_class), intent(INOUT) :: pdd

        pdd%now%snowpack_swe = 0.0_wp
        pdd%now%smb_ice      = 0.0_wp_acc
        pdd%now%runoff       = 0.0_wp_acc
        pdd%now%pdd_sum      = 0.0_wp_acc

        return

    end subroutine pdd_init_state

    subroutine pdd_reset_column(pdd,icol)
        ! Reset one column, for use when the active mask deactivates it
        ! (mirrors Chion.jl _reset_model_columns!, src/runtime.jl:178-190).

        implicit none

        type(pdd_class), intent(INOUT) :: pdd
        integer,         intent(IN)    :: icol

        pdd%now%snowpack_swe(icol) = 0.0_wp
        pdd%now%smb_ice(icol)      = 0.0_wp_acc
        pdd%now%runoff(icol)       = 0.0_wp_acc
        pdd%now%pdd_sum(icol)      = 0.0_wp_acc

        return

    end subroutine pdd_reset_column

    ! =====================================================================
    ! Degree-day integrals
    ! =====================================================================

    elemental function pdd_normal_cdf(z) result(cdf)
        ! Standard normal CDF.
        !
        ! DEVIATION (allowed by docs/PLAN.md 4.1): Chion.jl uses the
        ! Abramowitz-Stegun rational-polynomial approximation
        ! (_normal_cdf, pdd.jl:16-21), whose absolute error is ~7.5e-8.
        ! chion uses the exact identity Phi(z) = erfc(-z/sqrt(2))/2, which is
        ! what smbpal does (smb_pdd.f90:94). The two agree to ~1e-7, i.e. below
        ! sp round-off, and erfc is correct in the far tails where the
        ! polynomial is not. See docs/pdd_defects.md, deviations section.
        !
        ! Evaluated in wp_acc for the reason smbpal states explicitly: erfc
        ! underflows for strongly negative temperatures.

        implicit none

        real(wp_acc), intent(IN) :: z
        real(wp_acc) :: cdf

        cdf = 0.5_wp_acc*erfc(-z*PDD_INV_SQRT_2)

        return

    end function pdd_normal_cdf

    elemental function pdd_expected_positive_temperature(temp_c,sigma) result(teff)
        ! Expectation of max(T,0) for T ~ N(temp_c, sigma^2):
        !
        !     E[T+] = sigma*phi(z) + Tbar*Phi(z),   z = Tbar/sigma
        !
        ! Calov-Greve semi-analytical PDD integral. Chion.jl
        ! _expected_positive_temperature (pdd.jl:23-26); smbpal
        ! calc_temp_effective (smb_pdd.f90:56-111). The two are the same
        ! formula; only the Phi implementation differs (see pdd_normal_cdf).
        !
        ! Limits, all exercised by tests/test_pdd.f90:
        !     Tbar -> +inf   ->  Tbar
        !     Tbar  =  0     ->  sigma/sqrt(2*pi)
        !     Tbar -> -inf   ->  0
        !
        ! temp_c is in CELSIUS, not Kelvin.

        implicit none

        real(wp_acc), intent(IN) :: temp_c    ! [degC] mean air temperature
        real(wp_acc), intent(IN) :: sigma     ! [K]    std dev of variability
        real(wp_acc) :: teff                  ! [degC] effective positive temp

        ! Local variables
        real(wp_acc) :: z, pdf, cdf

        z   = temp_c/sigma
        pdf = PDD_INV_SQRT_2PI*exp(-0.5_wp_acc*z*z)
        cdf = pdd_normal_cdf(z)

        teff = sigma*pdf + temp_c*cdf

        return

    end function pdd_expected_positive_temperature

    function pdd_degree_days(air_temperature,dt_days,par,c) result(pdd)
        ! Positive degree days accumulated over one step [K d].
        !
        ! simple : max(T - T0, 0) * dt_days                (Chion.jl pdd.jl:8)
        ! pism   : dt_days * E[T+]                         (Chion.jl pdd.jl:28)
        !
        ! DEVIATION: Chion.jl chooses between these on dt_days alone
        ! (27 <= dt <= 32 -> pism). chion uses par%pdd_method. Defect D4.
        !
        ! DEVIATION: T0 comes from the constants struct, not the hard-coded
        ! 273.15 of pdd.jl:9 and :29.
        !
        ! Returned in wp_acc so that pdd_sum closes tightly and so that the
        ! erfc branch is not truncated to sp before being multiplied by dt.

        implicit none

        real(wp),                intent(IN) :: air_temperature   ! [K]
        real(wp),                intent(IN) :: dt_days           ! [d]
        type(pdd_par_class),     intent(IN) :: par
        type(chion_const_class), intent(IN) :: c
        real(wp_acc) :: pdd                                      ! [K d]

        ! Local variables
        real(wp_acc) :: temp_c, dt

        temp_c = real(air_temperature,wp_acc) - real(c%T0,wp_acc)
        dt     = real(dt_days,wp_acc)

        select case(par%pdd_method)

            case(CHION_PDD_SIMPLE)
                pdd = max(temp_c,0.0_wp_acc)*dt

            case(CHION_PDD_PISM)
                pdd = dt*pdd_expected_positive_temperature(temp_c, &
                                        real(par%temperature_sigma,wp_acc))

            case DEFAULT
                ! Unreachable if pdd_par_validate has been called. Kept so a
                ! corrupted flag fails loudly rather than silently selecting
                ! one of the two flavours (contrast Chion.jl's if/else style,
                ! see docs/porting_notes.md D5).
                write(io_unit_err,*) "pdd_degree_days:: Error: pdd_method not recognized."
                write(io_unit_err,*) "pdd_method = ", par%pdd_method
                stop "Program stopped."

        end select

        return

    end function pdd_degree_days

    elemental function pdd_step_mass(rate,dt_days,c) result(mass)
        ! Convert a precipitation rate [kg m-2 s-1] to a per-step mass
        ! [kg m-2], clipping negative rates. Chion.jl _pdd_step_mass
        ! (pdd.jl:11), except that seconds_per_day comes from the constants
        ! struct rather than the hard-coded 86400.0.
        !
        ! Returned in wp_acc: this feeds the cumulative accumulators directly
        ! and is the reference term of the mass-balance identity.

        implicit none

        real(wp),                intent(IN) :: rate      ! [kg m-2 s-1]
        real(wp),                intent(IN) :: dt_days   ! [d]
        type(chion_const_class), intent(IN) :: c
        real(wp_acc) :: mass                             ! [kg m-2]

        mass = max(real(rate,wp_acc),0.0_wp_acc) &
             * real(dt_days,wp_acc)*real(c%seconds_per_day,wp_acc)

        return

    end function pdd_step_mass

    ! =====================================================================
    ! The core
    ! =====================================================================

    subroutine pdd_column_apply(snowpack_swe,smb_ice,runoff,pdd_sum, &
                                snowfall,rainfall,pdd,par)
        ! The six-line PDD core, one column, one step.
        ! Chion.jl _pdd_apply_column_pdd! (pdd.jl:32-58).
        !
        ! ORDER MATTERS and is reproduced exactly:
        !
        !   pdd_sum        += pdd
        !   available_snow  = snowpack_swe + snowfall    ! snowfall BEFORE melt
        !   snow_melt       = min(available_snow, ddf_snow*pdd)
        !   remaining_pdd   = max(pdd - snow_melt/ddf_snow, 0)
        !   ice_melt        = ddf_ice*remaining_pdd
        !   refrozen        = refreezing_fraction*snow_melt
        !   snowpack_swe    = available_snow - snow_melt + refrozen
        !   smb_ice        += snowfall - snow_melt + refrozen - ice_melt
        !   runoff         += rainfall + snow_melt - refrozen + ice_melt
        !
        ! Snowfall is added to the reservoir BEFORE melt is applied, so snow
        ! falling within a step can be melted within the same step. Rainfall
        ! contributes to runoff ONLY -- it never enters the snowpack, never
        ! refreezes and never releases latent heat. Ice melt is drawn only
        ! from the degree days left over after the snow reservoir is exhausted,
        ! which is the same construction as smbpal's
        !     melt_ice = max(0,(melt_snow-acc)*mm_ice/mm_snow)
        ! (smb_pdd.f90:43), rewritten in terms of remaining degree days.
        !
        ! DEFECTS reproduced here deliberately, see docs/pdd_defects.md:
        !   D1  refrozen has no cold-content and no capacity limit; smbpal
        !       caps it with refrz = min(melt_snow, acc*f_refrz_max).
        !   D2  smb_ice is credited with the change in snowpack_swe as well as
        !       with the flux to the ice, so it is a whole-column mass change,
        !       not the "net mass forcing to the ice sheet" it is documented
        !       to be. This is why the three-reservoir identity
        !       dsmb + drunoff == snowfall + rainfall - dsnowpack does NOT
        !       hold; what holds instead is dsmb + drunoff == snowfall +
        !       rainfall. Both are asserted in tests/test_pdd.f90.
        !   D6  refrozen re-enters snowpack_swe as ordinary snow, so it can be
        !       melted again at ddf_snow and refrozen again next step.
        !
        ! All arithmetic in wp_acc. snowpack_swe is wp on entry and exit;
        ! everything between is dp so that the closure identity holds to
        ! round-off rather than to sp resolution.

        implicit none

        real(wp),            intent(INOUT) :: snowpack_swe   ! [kg m-2]
        real(wp_acc),        intent(INOUT) :: smb_ice        ! [kg m-2] cumulative
        real(wp_acc),        intent(INOUT) :: runoff         ! [kg m-2] cumulative
        real(wp_acc),        intent(INOUT) :: pdd_sum        ! [K d]    cumulative
        real(wp_acc),        intent(IN)    :: snowfall       ! [kg m-2] this step
        real(wp_acc),        intent(IN)    :: rainfall       ! [kg m-2] this step
        real(wp_acc),        intent(IN)    :: pdd            ! [K d]    this step
        type(pdd_par_class), intent(IN)    :: par

        ! Local variables
        real(wp_acc) :: ddf_snow, ddf_ice, f_refrz
        real(wp_acc) :: available_snow, snow_melt, remaining_pdd
        real(wp_acc) :: ice_melt, refrozen

        ddf_snow = real(par%ddf_snow,wp_acc)
        ddf_ice  = real(par%ddf_ice,wp_acc)
        f_refrz  = real(par%refreezing_fraction,wp_acc)

        pdd_sum = pdd_sum + pdd

        available_snow = real(snowpack_swe,wp_acc) + snowfall
        snow_melt      = min(available_snow,ddf_snow*pdd)

        ! ddf_snow > 0 is enforced by pdd_par_validate. The guard is kept
        ! because Chion.jl has one in its kernel path (pdd.jl:50) and none in
        ! its two vectorised paths (pdd.jl:229, :271) -- an inconsistency
        ! between three copies of the same six lines (defect D7).
        if (ddf_snow .gt. 0.0_wp_acc) then
            remaining_pdd = max(pdd - snow_melt/ddf_snow,0.0_wp_acc)
        else
            remaining_pdd = 0.0_wp_acc
        end if

        ice_melt = ddf_ice*remaining_pdd
        refrozen = f_refrz*snow_melt

        snowpack_swe = real(available_snow - snow_melt + refrozen,wp)
        smb_ice      = smb_ice + (snowfall - snow_melt + refrozen - ice_melt)
        runoff       = runoff  + (rainfall + snow_melt - refrozen + ice_melt)

        return

    end subroutine pdd_column_apply

    subroutine pdd_column_step(snowpack_swe,smb_ice,runoff,pdd_sum,forc,par,c)
        ! One column, one step, from the standard per-column forcing contract.
        ! Chion.jl _pdd_step_column! (pdd.jl:60-90), generalised to both PDD
        ! flavours.
        !
        ! PDD uses only three of the 22 chion_step_forcing_class fields
        ! (air_temperature, snowfall_rate, rainfall_rate) plus dt_days. It
        ! takes the full type anyway, so that the model kernels remain
        ! interchangeable behind the WP11 dispatcher.

        implicit none

        real(wp),                        intent(INOUT) :: snowpack_swe
        real(wp_acc),                    intent(INOUT) :: smb_ice
        real(wp_acc),                    intent(INOUT) :: runoff
        real(wp_acc),                    intent(INOUT) :: pdd_sum
        type(chion_step_forcing_class),  intent(IN)    :: forc
        type(pdd_par_class),             intent(IN)    :: par
        type(chion_const_class),         intent(IN)    :: c

        ! Local variables
        real(wp_acc) :: snowfall, rainfall, pdd

        snowfall = pdd_step_mass(forc%snowfall_rate,forc%dt_days,c)
        rainfall = pdd_step_mass(forc%rainfall_rate,forc%dt_days,c)
        pdd      = pdd_degree_days(forc%air_temperature,forc%dt_days,par,c)

        call pdd_column_apply(snowpack_swe,smb_ice,runoff,pdd_sum, &
                              snowfall,rainfall,pdd,par)

        return

    end subroutine pdd_column_step

    subroutine pdd_step(pdd,forc,dt_days,active_idx,c)
        ! Advance all ACTIVE columns by one step.
        !
        ! DEVIATION: Chion.jl threads active_indices down to the PDD stepper
        ! (src/runtime.jl:290-300) but no PDD method has a parameter to receive
        ! it -- every PDD path loops over all columns, so deactivated columns
        ! keep accumulating. chion honours the mask. Defect D8.

        implicit none

        type(pdd_class),           intent(INOUT) :: pdd
        type(chion_forcing_class), intent(IN)    :: forc
        real(wp),                  intent(IN)    :: dt_days
        integer,                   intent(IN)    :: active_idx(:)
        type(chion_const_class),   intent(IN)    :: c

        ! Local variables
        integer      :: i, icol
        real(wp_acc) :: snowfall, rainfall, dpdd

        if (forc%ncol .ne. pdd%now%ncol) then
            write(io_unit_err,*) "pdd_step:: Error: forcing column count must match state."
            write(io_unit_err,*) "forc%ncol, pdd%now%ncol = ", forc%ncol, pdd%now%ncol
            stop "Program stopped."
        end if

        !$omp parallel do default(shared) private(i,icol,snowfall,rainfall,dpdd)
        do i = 1, size(active_idx)

            icol = active_idx(i)

            snowfall = pdd_step_mass(forc%snowfall_rate(icol),dt_days,c)
            rainfall = pdd_step_mass(forc%rainfall_rate(icol),dt_days,c)
            dpdd     = pdd_degree_days(forc%air_temperature(icol),dt_days,pdd%par,c)

            call pdd_column_apply(pdd%now%snowpack_swe(icol), &
                                  pdd%now%smb_ice(icol),      &
                                  pdd%now%runoff(icol),       &
                                  pdd%now%pdd_sum(icol),      &
                                  snowfall,rainfall,dpdd,pdd%par)

        end do
        !$omp end parallel do

        return

    end subroutine pdd_step

end module snow_pdd
