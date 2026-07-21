module snow_itm
    ! Insolation-Temperature-Melt (ITM) snowpack model.
    !
    ! =====================================================================
    ! SOURCE OF TRUTH IS *NOT* Chion.jl. ITM does not exist there --
    ! build_model(:itm,...) deliberately errors (Chion.jl
    ! test/test_case_api.jl:538). This module is a port of
    !     ~/models/smbpal/src/smb_itm.f90
    ! (equivalently ~/models/yelmox/libs/smbpal/smb_itm.f90), which is the
    ! scheme actually in production use in yelmox. See docs/PLAN.md WP12.
    !
    ! The physics is UNCHANGED: calc_itm, calc_atmos_transmissivity,
    ! calc_albedo_planet, itm_c_lat, calc_albedo_surface, and the budget
    ! sequence of calc_snowpack_budget_step are line-for-line equivalents of
    ! the smbpal originals. Acceptance for this WP is numerical equivalence
    ! with smbpal (tests/test_itm.f90), because chion is to replace smbpal
    ! in yelmox and that equivalence is the evidence the migration is safe.
    ! =====================================================================
    !
    ! ---------------------------------------------------------------------
    ! INTERFACE CHANGE FROM smbpal: chion DOES NOT OWN AN INSOLATION MODULE
    !
    ! smbpal computes top-of-atmosphere insolation internally, calling
    !     now%S = calc_insol_day(day, lats, insol_time, fldr=par%insol_fldr)
    ! (smbpal/src/smbpal.f90:414) from the vendored `insolation` package,
    ! which also needs an orbital-parameter data directory on disk.
    !
    ! chion has no such dependency. Insolation arrives as FORCING:
    !     S = forc%q_sw_net        when forc%has_q_sw_net
    !     S = forc%shortwave_down  otherwise
    ! Both are [W m-2] in chion_step_forcing_class, the same units and the
    ! same meaning (downward shortwave at the top of the atmosphere, since
    ! ITM applies its own atmospheric transmissivity `atrans` to it).
    !
    ! *** MIGRATION REQUIREMENT FOR yelmox (WP19 phase A) ***
    ! The host must now supply insolation itself, either from the vendored
    ! yelmox/libs/insol/ package (the like-for-like choice, reproducing
    ! smbpal exactly) or from the climate forcing. If the host writes
    ! shortwave_down but leaves has_q_sw_net .FALSE., ITM behaves as smbpal
    ! did. If the host prescribes q_sw_net, note that ITM still multiplies
    ! by atrans*(1-alb_s) -- q_sw_net is used as a TOA-like driver here, not
    ! as an already-net surface flux. Prescribing a true net surface flux
    ! would double-count the albedo. This is inherited from the ITM
    ! formulation, not introduced here.
    ! ---------------------------------------------------------------------
    !
    ! ADDITIONAL HOST-SUPPLIED INPUTS. chion_step_forcing_class was designed
    ! for BESSI and does not carry ITM's three geometry/vegetation inputs:
    ! surface elevation z_srf, ice thickness H_ice, and the annual
    ! positive-degree-day total PDDs. They are passed to itm_step as
    ! explicit arguments rather than being bolted onto the shared forcing
    ! type, which would change a contract shared with BESSI and PDD.
    ! PDDs in smbpal is a WHOLE-YEAR total recomputed once per year
    ! (smbpal.f90:378-385), not a per-step quantity; the host must supply it
    ! with the same meaning.
    !
    ! UNITS. smbpal works in [mm w.e. d-1]; chion forcing is [kg m-2 s-1].
    ! The conversion (x SEC_DAY) is done inside itm_step, so all ITM state
    ! remains in smbpal's units and is directly comparable.

    use chion_defs, only : wp, wp_acc, io_unit_err, chion_step_forcing_class
    use nml,        only : nml_read

    implicit none

    private

    ! Local physical constants, values taken verbatim from
    ! smbpal/src/smb_itm.f90:12-14. NOTE L_M = 3.35e5 differs from chion's
    ! chion_const_class%Lm = 3.34e5 (Chion.jl's value); the smbpal value is
    ! kept here because ITM's melt rate is calibrated against it and the
    ! acceptance test is equivalence with smbpal. Do not "reconcile" these
    ! without re-tuning itm_c/itm_t.
    real(wp), parameter, public :: ITM_SEC_DAY = 86400.0_wp   ! [s d-1]
    real(wp), parameter, public :: ITM_RHO_W   = 1.0e3_wp     ! [kg m-3] pure water
    real(wp), parameter, public :: ITM_L_M     = 3.35e5_wp    ! [J kg-1] latent heat of melting

    type itm_par_class
        ! smbpal itm_par_class (smb_itm.f90:16-27), same names and order, plus
        ! firn_fac (see below).

        real(wp) :: trans_a               ! [1] atmospheric transmissivity, intercept
        real(wp) :: trans_b               ! [m-1/2] transmissivity elevation slope
        real(wp) :: trans_c               ! [1] unused by the active transmissivity form

        real(wp) :: itm_c                 ! [W m-2] ITM offset
        real(wp) :: itm_t                 ! [W m-2 K-1] ITM temperature coefficient
        real(wp) :: itm_b                 ! [W m-2 deg-1] latitude slope of itm_c
        real(wp) :: itm_lat0              ! [deg N] reference latitude; |lat0| >= 90 disables

        real(wp) :: H_snow_max            ! [mm w.e.] snowpack cap
        real(wp) :: Pmaxfrac              ! [1] refreezing factor

        real(wp) :: H_snow_crit_desert    ! [mm w.e.] critical snow depth, PDDs <= 100
        real(wp) :: H_snow_crit_forest    ! [mm w.e.] critical snow depth, PDDs >= 1000
        real(wp) :: melt_crit             ! [mm d-1] melt rate at which snow albedo goes wet

        real(wp) :: alb_ocean             ! [1]
        real(wp) :: alb_land              ! [1]
        real(wp) :: alb_forest            ! [1]
        real(wp) :: alb_ice               ! [1]
        real(wp) :: alb_snow_dry          ! [1]
        real(wp) :: alb_snow_wet          ! [1]

        ! DEVIATION: firn_fac lives in smbpal_param_class, not itm_par_class,
        ! because smbpal applies calc_temp_surf once per year outside the ITM
        ! module (smbpal.f90:457). chion computes tsrf inside itm_step, so
        ! the parameter has to travel with the ITM parameters.
        real(wp) :: firn_fac              ! [K (mm w.e.)-1] firn warming per unit net refreezing
    end type itm_par_class

    type itm_state_class
        ! Per-column scalars. ncol columns, packed list (docs/PLAN.md section 1).
        !
        ! PRECISION (docs/PLAN.md section 3.1). smbpal's
        ! calc_snowpack_budget_step produces *rates* for the step: melt,
        ! runoff, refrz, smb, smbi and melt_net are all divided by dt on the
        ! way out (smb_itm.f90:196-202) and are overwritten every call.
        ! Nothing inside the routine accumulates; smbpal's accumulation lives
        ! in smbpal_average, outside. So none of the eight named state fields
        ! is cumulative, and all eight are wp -- promoting them to wp_acc
        ! would buy nothing and would change the values chion reports.
        !
        ! H_snow is prognostic but bounded above by H_snow_max (default
        ! 5000 mm w.e.) and never sums small increments onto a large total
        ! without also being drawn down, so it stays wp too. It is also the
        ! quantity the smbpal equivalence test compares, and smbpal itself
        ! carries it in sp.
        !
        ! The wp_acc block below is the genuinely cumulative set. These have
        ! no smbpal counterpart -- they exist because chion's other models
        ! (BESSI, PDD) expose cumulative smb_ice/runoff/melt/refreezing and
        ! WP11/WP14 need one shape for all three. They are summed every step
        ! and never reset, which is exactly the case docs/porting_notes.md D1
        ! shows to be unsafe in sp (80 kg m-2 drift over 100 yr), hence
        ! wp_acc. They are pure additions and do not affect the smbpal
        ! comparison.

        integer :: ncol

        ! --- The eight state fields, smbpal semantics, per-step values -----
        real(wp), allocatable :: H_snow(:)    ! [mm w.e.] snowpack thickness (prognostic)
        real(wp), allocatable :: alb_s(:)     ! [1] surface albedo
        real(wp), allocatable :: smb(:)       ! [mm w.e. d-1] total surface mass balance
        real(wp), allocatable :: smbi(:)      ! [mm w.e. d-1] mass balance seen by the ice sheet
        real(wp), allocatable :: melt(:)      ! [mm w.e. d-1] total melt (snow + ice)
        real(wp), allocatable :: runoff(:)    ! [mm w.e. d-1] runoff
        real(wp), allocatable :: refrz(:)     ! [mm w.e. d-1] refreezing
        real(wp), allocatable :: tsrf(:)      ! [K] surface temperature

        ! Net melt, smbpal's melt_net. Not in the eight, but it is an output
        ! of calc_snowpack_budget_step and is what tsrf is built from, so it
        ! is retained rather than discarded as a local.
        real(wp), allocatable :: melt_net(:)  ! [mm w.e. d-1] refrz - melt

        ! --- Cumulative accumulators (chion addition; must be wp_acc) ------
        real(wp_acc), allocatable :: smb_cum(:)     ! [mm w.e.]
        real(wp_acc), allocatable :: smbi_cum(:)    ! [mm w.e.]
        real(wp_acc), allocatable :: melt_cum(:)    ! [mm w.e.]
        real(wp_acc), allocatable :: runoff_cum(:)  ! [mm w.e.]
        real(wp_acc), allocatable :: refrz_cum(:)   ! [mm w.e.]
    end type itm_state_class

    type itm_class
        type(itm_par_class)   :: par
        type(itm_state_class) :: now
    end type itm_class

    public :: itm_par_class
    public :: itm_state_class
    public :: itm_class

    public :: itm_par_load
    public :: itm_alloc
    public :: itm_dealloc
    public :: itm_init_state
    public :: itm_step

    ! Physics, exposed so the acceptance test and any host diagnostic can
    ! call them directly.
    public :: calc_itm
    public :: calc_atmos_transmissivity
    public :: calc_albedo_planet
    public :: calc_albedo_surface
    public :: itm_c_lat
    public :: calc_temp_surf

contains

    subroutine itm_par_load(par,filename,init,group)
        ! smbpal itm_par_load (smb_itm.f90:35-77), one nml_read per
        ! parameter, plus firn_fac.

        implicit none

        type(itm_par_class),        intent(INOUT) :: par
        character(len=*),           intent(IN)    :: filename
        logical,          optional, intent(IN)    :: init
        character(len=*), optional, intent(IN)    :: group

        ! Local variables
        logical            :: init_pars
        character(len=56)  :: nml_group

        nml_group = "itm"
        if (present(group)) nml_group = trim(group)

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        call nml_read(filename,nml_group,"trans_a",           par%trans_a,           init=init_pars)
        call nml_read(filename,nml_group,"trans_b",           par%trans_b,           init=init_pars)
        call nml_read(filename,nml_group,"trans_c",           par%trans_c,           init=init_pars)
        call nml_read(filename,nml_group,"itm_c",             par%itm_c,             init=init_pars)
        call nml_read(filename,nml_group,"itm_t",             par%itm_t,             init=init_pars)
        call nml_read(filename,nml_group,"itm_b",             par%itm_b,             init=init_pars)
        call nml_read(filename,nml_group,"itm_lat0",          par%itm_lat0,          init=init_pars)
        call nml_read(filename,nml_group,"H_snow_max",        par%H_snow_max,        init=init_pars)
        call nml_read(filename,nml_group,"Pmaxfrac",          par%Pmaxfrac,          init=init_pars)
        call nml_read(filename,nml_group,"H_snow_crit_desert",par%H_snow_crit_desert,init=init_pars)
        call nml_read(filename,nml_group,"H_snow_crit_forest",par%H_snow_crit_forest,init=init_pars)
        call nml_read(filename,nml_group,"melt_crit",         par%melt_crit,         init=init_pars)
        call nml_read(filename,nml_group,"alb_ocean",         par%alb_ocean,         init=init_pars)
        call nml_read(filename,nml_group,"alb_land",          par%alb_land,          init=init_pars)
        call nml_read(filename,nml_group,"alb_forest",        par%alb_forest,        init=init_pars)
        call nml_read(filename,nml_group,"alb_ice",           par%alb_ice,           init=init_pars)
        call nml_read(filename,nml_group,"alb_snow_dry",      par%alb_snow_dry,      init=init_pars)
        call nml_read(filename,nml_group,"alb_snow_wet",      par%alb_snow_wet,      init=init_pars)
        call nml_read(filename,nml_group,"firn_fac",          par%firn_fac,          init=init_pars)

        return

    end subroutine itm_par_load

    subroutine itm_alloc(itm,ncol)
        ! Allocate the per-column state. No physics, no initial values --
        ! those come from itm_init_state, following the yelmo convention.

        implicit none

        type(itm_class), intent(INOUT) :: itm
        integer,         intent(IN)    :: ncol

        call itm_dealloc(itm)

        if (ncol .le. 0) then
            write(io_unit_err,*) "itm_alloc:: Error: ncol must be positive."
            write(io_unit_err,*) "ncol = ", ncol
            stop "Program stopped."
        end if

        itm%now%ncol = ncol

        allocate(itm%now%H_snow(ncol))
        allocate(itm%now%alb_s(ncol))
        allocate(itm%now%smb(ncol))
        allocate(itm%now%smbi(ncol))
        allocate(itm%now%melt(ncol))
        allocate(itm%now%runoff(ncol))
        allocate(itm%now%refrz(ncol))
        allocate(itm%now%tsrf(ncol))
        allocate(itm%now%melt_net(ncol))

        allocate(itm%now%smb_cum(ncol))
        allocate(itm%now%smbi_cum(ncol))
        allocate(itm%now%melt_cum(ncol))
        allocate(itm%now%runoff_cum(ncol))
        allocate(itm%now%refrz_cum(ncol))

        return

    end subroutine itm_alloc

    subroutine itm_dealloc(itm)

        implicit none

        type(itm_class), intent(INOUT) :: itm

        if (allocated(itm%now%H_snow))     deallocate(itm%now%H_snow)
        if (allocated(itm%now%alb_s))      deallocate(itm%now%alb_s)
        if (allocated(itm%now%smb))        deallocate(itm%now%smb)
        if (allocated(itm%now%smbi))       deallocate(itm%now%smbi)
        if (allocated(itm%now%melt))       deallocate(itm%now%melt)
        if (allocated(itm%now%runoff))     deallocate(itm%now%runoff)
        if (allocated(itm%now%refrz))      deallocate(itm%now%refrz)
        if (allocated(itm%now%tsrf))       deallocate(itm%now%tsrf)
        if (allocated(itm%now%melt_net))   deallocate(itm%now%melt_net)
        if (allocated(itm%now%smb_cum))    deallocate(itm%now%smb_cum)
        if (allocated(itm%now%smbi_cum))   deallocate(itm%now%smbi_cum)
        if (allocated(itm%now%melt_cum))   deallocate(itm%now%melt_cum)
        if (allocated(itm%now%runoff_cum)) deallocate(itm%now%runoff_cum)
        if (allocated(itm%now%refrz_cum))  deallocate(itm%now%refrz_cum)

        itm%now%ncol = 0

        return

    end subroutine itm_dealloc

    subroutine itm_init_state(itm,H_snow)
        ! Cold start. smbpal seeds the snowpack at its maximum thickness
        ! (smbpal.f90:117: smb%now%H_snow = smb%par%itm%H_snow_max), which is
        ! preserved as the default so a chion ITM run starting from scratch
        ! matches an smbpal run starting from scratch. Pass H_snow to
        ! override, e.g. from a restart.

        implicit none

        type(itm_class),    intent(INOUT) :: itm
        real(wp), optional, intent(IN)    :: H_snow(:)

        if (.not. allocated(itm%now%H_snow)) then
            write(io_unit_err,*) "itm_init_state:: Error: state not allocated. Call itm_alloc first."
            stop "Program stopped."
        end if

        if (present(H_snow)) then
            if (size(H_snow) .ne. itm%now%ncol) then
                write(io_unit_err,*) "itm_init_state:: Error: H_snow must have length ncol."
                write(io_unit_err,*) "ncol, size(H_snow) = ", itm%now%ncol, size(H_snow)
                stop "Program stopped."
            end if
            itm%now%H_snow = H_snow
        else
            itm%now%H_snow = itm%par%H_snow_max
        end if

        itm%now%alb_s    = itm%par%alb_snow_dry
        itm%now%smb      = 0.0_wp
        itm%now%smbi     = 0.0_wp
        itm%now%melt     = 0.0_wp
        itm%now%runoff   = 0.0_wp
        itm%now%refrz    = 0.0_wp
        itm%now%melt_net = 0.0_wp
        itm%now%tsrf     = 273.15_wp

        itm%now%smb_cum    = 0.0_wp_acc
        itm%now%smbi_cum   = 0.0_wp_acc
        itm%now%melt_cum   = 0.0_wp_acc
        itm%now%runoff_cum = 0.0_wp_acc
        itm%now%refrz_cum  = 0.0_wp_acc

        return

    end subroutine itm_init_state

    subroutine itm_step(itm,icol,forc,z_srf,H_ice,PDDs)
        ! Advance one column by one step.
        !
        ! Port of smbpal calc_snowpack_budget_step (smb_itm.f90:80-206), in
        ! the same order and with the same expressions. smbpal's routine is
        ! `elemental` over (nx,ny); chion calls it once per column from the
        ! dispatcher's loop, so it is a plain subroutine here.
        !
        ! Note on smbpal's definitions, restated (smb_itm.f90:85-87):
        !     SMB    = sf   + rf   - runoff     [kg m-2 d-1] == [mm d-1]
        !     runoff = rain + melt - refrz      [kg m-2 d-1] == [mm d-1]
        !
        ! smbpal requires dt >= 1 day (smb_itm.f90:92). That is not enforced
        ! here either -- see the "To avoid numerical issues with dt>1" clamp
        ! at the H_snow update, which is the mechanism smbpal relies on.

        implicit none

        type(itm_class), intent(INOUT) :: itm
        integer,         intent(IN)    :: icol
        type(chion_step_forcing_class), intent(IN) :: forc
        real(wp),        intent(IN)    :: z_srf   ! [m] surface elevation
        real(wp),        intent(IN)    :: H_ice   ! [m] ice thickness
        real(wp),        intent(IN)    :: PDDs    ! [K d] annual positive degree days

        ! Local variables
        real(wp) :: dt, lat, t2m, S
        real(wp) :: pr, sf, rf
        real(wp) :: H_snow, alb_s
        real(wp) :: smb, smbi, melt, runoff, refrz, melt_net
        real(wp) :: itm_c_now, melt_pot, atrans, rfac
        real(wp) :: melted_snow, melted_ice, snow_to_ice
        real(wp) :: refrz_rain, refrz_snow

        if (icol .lt. 1 .or. icol .gt. itm%now%ncol) then
            write(io_unit_err,*) "itm_step:: Error: column index out of range."
            write(io_unit_err,*) "icol, ncol = ", icol, itm%now%ncol
            stop "Program stopped."
        end if

        ! === Unpack the forcing =========================================
        dt  = forc%dt_days
        lat = forc%latitude_deg
        t2m = forc%air_temperature

        ! Insolation. chion has no insolation module; the host supplies it.
        ! See the module header for the yelmox migration consequence.
        if (forc%has_q_sw_net) then
            S = forc%q_sw_net
        else
            S = forc%shortwave_down
        end if

        ! [kg m-2 s-1] -> [mm w.e. d-1], smbpal's working units.
        sf = forc%snowfall_rate*ITM_SEC_DAY
        rf = forc%rainfall_rate*ITM_SEC_DAY
        pr = sf + rf
        ! smbpal receives (pr,sf) and forms rf = pr - sf (smb_itm.f90:107).
        ! chion receives snowfall and rainfall separately and forms pr
        ! instead; algebraically identical, and it avoids a cancellation.

        H_snow = itm%now%H_snow(icol)

        ! === Preliminary surface albedo =================================
        ! No `melt` argument yet, so calc_albedo_surface assumes high melt
        ! (melt_crit + 1), which is what avoids an albedo jump at melt onset.
        alb_s = calc_albedo_surface(itm%par,z_srf,H_ice,H_snow,PDDs)

        ! === Accumulation ===============================================
        H_snow = H_snow + sf*dt

        ! === Potential melt from the ITM scheme =========================
        ! NOTE trans_c is read and stored but unused: smbpal's three-argument
        ! transmissivity call is commented out (smb_itm.f90:116-117) in
        ! favour of the two-argument sqrt-elevation form. Preserved.
        atrans = calc_atmos_transmissivity(z_srf,itm%par%trans_a,itm%par%trans_b)

        if (abs(itm%par%itm_lat0) .lt. 90.0_wp) then
            itm_c_now = itm_c_lat(itm%par%itm_c,itm%par%itm_b,itm%par%itm_lat0,lat)
        else
            itm_c_now = itm%par%itm_c
        end if

        ! calc_itm takes air temperature in CELSIUS. smbpal hard-codes the
        ! 273.15 conversion at the call site (smb_itm.f90:123); preserved
        ! rather than routed through chion_const_class%T0, so that a host
        ! that retunes T0 cannot silently retune ITM's melt.
        melt_pot = calc_itm(S,t2m-273.15_wp,alb_s,atrans,itm_c_now,itm%par%itm_t)

        ! === Partition melt between snow and ice ========================
        if (melt_pot*dt .gt. H_snow) then
            ! All snow melted; the remaining energy melts ice.
            melted_snow = H_snow
            melted_ice  = melt_pot*dt - H_snow
        else
            ! Snow melt uses all the energy; none left for ice.
            melted_snow = melt_pot*dt
            melted_ice  = 0.0_wp
        end if

        melt = melted_snow + melted_ice

        H_snow = H_snow - melted_snow

        ! smbpal: "To avoid numerical issues with dt>1"
        H_snow = max(H_snow,0.0_wp)

        ! === Albedo again, now knowing the actual melt rate =============
        alb_s = calc_albedo_surface(itm%par,z_srf,H_ice,H_snow,PDDs,melt=melt/dt)

        ! === Refreezing =================================================
        ! Refreezing fraction increases linearly to 1 as H_snow reaches 1 m
        ! w.e. NOTE the literal 1e-3 guard on pr and the literal 1e3 scale:
        ! both are smbpal's, carried over verbatim.
        rfac = itm%par%Pmaxfrac*sf/max(1.0e-3_wp,pr)
        rfac = rfac + min(1.0_wp,H_snow/1.0e3_wp)*(1.0_wp - rfac)

        ! Capacity is the snowpack thickness. Rain claims it first.
        refrz_rain = min(rf*dt*rfac,H_snow)
        refrz_snow = min(melted_snow*rfac,H_snow-refrz_rain)
        refrz      = refrz_snow + refrz_rain

        ! Refrozen mass leaves the snowpack as ice; excess above the cap
        ! also becomes ice.
        snow_to_ice = refrz
        H_snow      = H_snow - refrz

        if (H_snow .gt. itm%par%H_snow_max) then
            snow_to_ice = snow_to_ice + (H_snow - itm%par%H_snow_max)
            H_snow      = itm%par%H_snow_max
        end if

        ! === Budgets ====================================================
        runoff = (melted_snow - refrz_snow) + (rf*dt - refrz_rain) + melted_ice

        smb = (sf + rf)*dt - runoff

        ! smbi = -melted_ice for negative mass balance,
        !      =     new_ice for positive mass balance.
        smbi = snow_to_ice + refrz - melted_ice

        ! Net melt, for the surface-temperature adjustment.
        if (H_ice .gt. 0.0_wp) then
            melt_net = refrz - melt
        else
            ! refrz is zero here, but smbpal keeps the term for consistency.
            melt_net = refrz - melted_snow
        end if

        ! === Back to daily rates [mm w.e. d-1] ==========================
        melt     = melt/dt
        melt_net = melt_net/dt
        runoff   = runoff/dt
        refrz    = refrz/dt
        smb      = smb/dt
        smbi     = smbi/dt

        ! === Store ======================================================
        itm%now%H_snow(icol)   = H_snow
        itm%now%alb_s(icol)    = alb_s
        itm%now%smb(icol)      = smb
        itm%now%smbi(icol)     = smbi
        itm%now%melt(icol)     = melt
        itm%now%runoff(icol)   = runoff
        itm%now%refrz(icol)    = refrz
        itm%now%melt_net(icol) = melt_net

        ! Surface temperature. smbpal applies calc_temp_surf once per year to
        ! the ANNUAL MEAN t2m and melt_net (smbpal.f90:457); chion applies it
        ! per step to the step values. Over a full year of equal-length steps
        ! the two agree only where the min(T0,...) cap is inactive -- the cap
        ! makes it a nonlinear function, so a per-step mean is not the same
        ! as a function of the mean. Recorded as a deviation; the host can
        ! recover smbpal's exact behaviour by averaging tsrf's inputs itself.
        itm%now%tsrf(icol) = calc_temp_surf(t2m,H_ice,melt_net,itm%par%firn_fac)

        ! Cumulative accumulators, in wp_acc. Rates x dt, so these are the
        ! integrated quantities in [mm w.e.].
        itm%now%smb_cum(icol)    = itm%now%smb_cum(icol)    + real(smb,   wp_acc)*real(dt,wp_acc)
        itm%now%smbi_cum(icol)   = itm%now%smbi_cum(icol)   + real(smbi,  wp_acc)*real(dt,wp_acc)
        itm%now%melt_cum(icol)   = itm%now%melt_cum(icol)   + real(melt,  wp_acc)*real(dt,wp_acc)
        itm%now%runoff_cum(icol) = itm%now%runoff_cum(icol) + real(runoff,wp_acc)*real(dt,wp_acc)
        itm%now%refrz_cum(icol)  = itm%now%refrz_cum(icol)  + real(refrz, wp_acc)*real(dt,wp_acc)

        return

    end subroutine itm_step

    ! =====================================================================
    ! Physics -- unchanged from smbpal
    ! =====================================================================

    pure function calc_albedo_surface(par,z_srf,H_ice,H_snow,PDDs,melt) result(alb)
        ! smbpal calc_albedo_surface (smb_itm.f90:213-275).
        !
        ! Critical snow depth follows vegetation type as diagnosed from PDDs:
        !     PDDs <= 100     desert, H_snow_crit_desert (10 mm)
        !     100..1000       tundra, linear in PDDs
        !     PDDs >= 1000    forest, H_snow_crit_forest (100 mm)
        !
        ! `melt` is optional exactly as in smbpal: when absent the routine
        ! assumes melt above melt_crit, i.e. wet snow. That default is what
        ! keeps the first (pre-melt) albedo call from producing a jump in
        ! albedo at the onset of the melt season.

        implicit none

        type(itm_par_class), intent(IN) :: par
        real(wp),            intent(IN) :: z_srf    ! [m]
        real(wp),            intent(IN) :: H_ice    ! [m]
        real(wp),            intent(IN) :: H_snow   ! [mm w.e.]
        real(wp),            intent(IN) :: PDDs     ! [K d]
        real(wp), optional,  intent(IN) :: melt     ! [mm w.e. d-1]
        real(wp) :: alb

        ! Local variables
        real(wp) :: H_snow_crit, depth, as_snow, alb_bg, melt_now

        if (PDDs .le. 100.0_wp) then
            H_snow_crit = par%H_snow_crit_desert
        else if (PDDs .le. 1000.0_wp) then
            H_snow_crit = par%H_snow_crit_desert &
                        + (par%H_snow_crit_forest - par%H_snow_crit_desert) &
                          *(PDDs - 100.0_wp)/(1000.0_wp - 100.0_wp)
        else
            H_snow_crit = par%H_snow_crit_forest
        end if

        depth = min(H_snow/H_snow_crit,1.0_wp)

        melt_now = par%melt_crit + 1.0_wp
        if (present(melt)) melt_now = melt

        ! Background (snow-free) albedo from the ground type.
        ! NOTE the three branches are ordered, and the second tests
        ! H_ice .eq. 0.0 exactly -- an ice thickness of 1e-6 m selects the
        ! ice branch. smbpal's test, preserved.
        if (z_srf .le. 0.0_wp) then
            alb_bg = par%alb_ocean
        else if (z_srf .gt. 0.0_wp .and. H_ice .eq. 0.0_wp) then
            alb_bg = par%alb_land  *(1000.0_wp - min(PDDs,1000.0_wp))/1000.0_wp &
                   + par%alb_forest*(           min(PDDs,1000.0_wp))/1000.0_wp
        else
            alb_bg = par%alb_ice
        end if

        as_snow = par%alb_snow_dry
        if (melt_now .gt. par%melt_crit) as_snow = par%alb_snow_wet

        alb = alb_bg + depth*(as_snow - alb_bg)

        return

    end function calc_albedo_surface

    pure function calc_albedo_planet(alb_s,a,b) result(alb_p)
        ! smbpal calc_albedo_planet (smb_itm.f90:277-288). Not used by the
        ! budget step; retained because it is part of the ITM scheme's public
        ! surface and hosts may want it for radiation diagnostics.

        implicit none

        real(wp), intent(IN) :: alb_s     ! [1] surface albedo
        real(wp), intent(IN) :: a         ! [1]
        real(wp), intent(IN) :: b         ! [1]
        real(wp) :: alb_p                 ! [1] planetary albedo

        alb_p = a + b*alb_s

        return

    end function calc_albedo_planet

    pure function calc_atmos_transmissivity(z_srf,a,b) result(at)
        ! smbpal calc_atmos_transmissivity (smb_itm.f90:290-302):
        !     at = a + b*max(z_srf,0)**0.5
        !
        ! DEVIATION (round-off only): smbpal writes `max(z_srf,0.d0)**0.5`,
        ! mixing an sp variable with a dp literal, so gfortran evaluates the
        ! max and the power in dp and rounds the result back to sp. chion
        ! keeps the whole expression in wp. The difference is at sp round-off
        ! and is measured by the equivalence test.

        implicit none

        real(wp), intent(IN) :: z_srf     ! [m]
        real(wp), intent(IN) :: a         ! [1]
        real(wp), intent(IN) :: b         ! [m-1/2]
        real(wp) :: at                    ! [1]

        at = a + b*max(z_srf,0.0_wp)**0.5_wp

        return

    end function calc_atmos_transmissivity

    pure function calc_itm(S,t2m,alb_s,atrans,c,t) result(melt)
        ! smbpal calc_itm (smb_itm.f90:311-327). Insolation-temperature-melt:
        ! potential melt from absorbed shortwave plus a linear temperature
        ! term, clipped at zero.
        !
        ! t2m is in DEGREES CELSIUS.
        !
        ! DEVIATION (round-off only): smbpal's `(atrans*(1.d0 - alb_s)*S ...)`
        ! and `max(melt,0.d0)*sec_day*1d3` promote to dp via the literals.
        ! chion evaluates in wp throughout.

        implicit none

        real(wp), intent(IN) :: S         ! [W m-2] insolation
        real(wp), intent(IN) :: t2m       ! [degC] air temperature
        real(wp), intent(IN) :: alb_s     ! [1] surface albedo
        real(wp), intent(IN) :: atrans    ! [1] atmospheric transmissivity
        real(wp), intent(IN) :: c         ! [W m-2] ITM offset
        real(wp), intent(IN) :: t         ! [W m-2 K-1] ITM temperature coefficient
        real(wp) :: melt                  ! [mm w.e. d-1]

        ! Potential melt [m s-1]
        melt = (atrans*(1.0_wp - alb_s)*S + c + t*t2m)/(ITM_RHO_W*ITM_L_M)

        ! [m s-1] -> [mm d-1], positive melt only
        melt = max(melt,0.0_wp)*ITM_SEC_DAY*1.0e3_wp

        return

    end function calc_itm

    pure function itm_c_lat(c,b,lat0,lat) result(c2D)
        ! smbpal itm_c_lat (smb_itm.f90:330-343). Linear latitude dependence
        ! of the ITM offset.

        implicit none

        real(wp), intent(IN) :: c         ! [W m-2] offset at lat0
        real(wp), intent(IN) :: b         ! [W m-2 deg-1]
        real(wp), intent(IN) :: lat0      ! [deg N]
        real(wp), intent(IN) :: lat       ! [deg N]
        real(wp) :: c2D                   ! [W m-2]

        c2D = c + b*(lat - lat0)

        return

    end function itm_c_lat

    pure function calc_temp_surf(tann,H_ice,melt_net,fac) result(ts)
        ! smbpal calc_temp_surf (smbpal.f90:644-661). Lives in smbpal.f90,
        ! not smb_itm.f90, but tsrf is part of chion's ITM state and yelmox
        ! consumes it (yelmox/libs/yelmox_domain.f90:1167), so it is ported
        ! here.
        !
        ! NOTE the 273.15 cap is hard-coded in smbpal; preserved.

        implicit none

        real(wp), intent(IN) :: tann      ! [K] air temperature
        real(wp), intent(IN) :: H_ice     ! [m]
        real(wp), intent(IN) :: melt_net  ! [mm w.e. d-1] refrz - melt
        real(wp), intent(IN) :: fac       ! [K (mm w.e.)-1]
        real(wp) :: ts                    ! [K]

        if (H_ice .gt. 0.0_wp) then
            ! Positive melt_net (net refreezing) warms the firn.
            ts = tann + fac*max(0.0_wp,melt_net)
            ! Cap at the freezing point on the ice sheet.
            ts = min(273.15_wp,ts)
        else
            ts = tann
        end if

        return

    end function calc_temp_surf

end module snow_itm
