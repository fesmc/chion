program test_defs
    ! WP1 acceptance test: chion_defs compiles standalone, the defaults match
    ! Chion.jl, the scheme-name mappings work, and the forcing/grid containers
    ! allocate and pack correctly.

    use chion_defs

    implicit none

    type(chion_const_class)        :: c
    type(chion_forcing_class)      :: forc
    type(chion_grid_class)         :: grd
    type(chion_step_forcing_class) :: fc

    integer, parameter :: ncol = 6
    logical :: active(ncol)
    integer :: nfail

    nfail = 0

    write(*,"(a)") "=========================================================="
    write(*,"(a)") " chion WP1 acceptance test"
    write(*,"(a)") "=========================================================="
    write(*,*)

    ! --- Precision -------------------------------------------------------
    ! wp = sp for state and interfaces (yelmo-compatible); wp_acc = dp for
    ! cumulative accumulators. See the precision policy in chion_defs.
    call check("wp is single precision",  wp .eq. kind(1.0),   nfail)
    call check("wp == sp",                wp .eq. sp,          nfail)
    call check("wp_acc is double",        wp_acc .eq. kind(1.d0), nfail)
    call check("wp_acc /= wp",            wp_acc .ne. wp,      nfail)

    ! --- Tolerances ------------------------------------------------------
    call check("TOL_TINY        = 1e-12", abs(TOL_TINY        - 1.0e-12_wp_acc) .lt. 1.0e-24_wp_acc, nfail)
    call check("TOL_EMPTY_LAYER = 1e-10", abs(TOL_EMPTY_LAYER - 1.0e-10_wp_acc) .lt. 1.0e-22_wp_acc, nfail)
    call check("tolerances are distinct", TOL_TINY .ne. TOL_EMPTY_LAYER, nfail)
    call check("tolerances are dp",       kind(TOL_TINY) .eq. wp_acc, nfail)

    ! The reason TOL_TINY must be compared against dp-computed quantities:
    ! in sp arithmetic it vanishes against anything of order 1 or larger.
    call check("TOL_TINY vanishes in sp arithmetic at O(1)", &
               1.0_wp + real(TOL_TINY,wp) .eq. 1.0_wp, nfail)
    call check("TOL_TINY survives in dp arithmetic at O(1)", &
               1.0_wp_acc + TOL_TINY .ne. 1.0_wp_acc, nfail)

    ! --- Constants defaults, against Chion.jl/src/constants.jl:173-200 ---
    call chion_const_init(c)

    call check_val("rho_s",           c%rho_s,          315.0_wp,       nfail)
    call check_val("rho_i",           c%rho_i,          917.0_wp,       nfail)
    call check_val("rho_w",           c%rho_w,          1000.0_wp,      nfail)
    call check_val("rho_s_a",         c%rho_s_a,        109.0_wp,       nfail)
    call check_val("rho_s_b",         c%rho_s_b,        6.0_wp,         nfail)
    call check_val("rho_s_c",         c%rho_s_c,        26.0_wp,        nfail)
    call check_val("Ki",              c%Ki,             2.1_wp,         nfail)
    call check_val("ci",              c%ci,             2110.0_wp,      nfail)
    call check_val("cw",              c%cw,             4181.0_wp,      nfail)
    call check_val("Lm",              c%Lm,             334000.0_wp,    nfail)
    call check_val("Lv",              c%Lv,             2.501e6_wp,     nfail)
    call check_val("cp_air",          c%cp_air,         1003.0_wp,      nfail)
    call check_val("D_sh",            c%D_sh,           10.0_wp,        nfail)
    call check_val("alpha_dry",       c%alpha_dry,      0.81_wp,        nfail)
    call check_val("alpha_wet",       c%alpha_wet,      0.70_wp,        nfail)
    call check_val("alpha_ice",       c%alpha_ice,      0.30_wp,        nfail)
    call check_val("max_lwc_albedo",  c%max_lwc_albedo, 0.10_wp,        nfail)
    call check_val("eps_air",         c%eps_air,        0.80_wp,        nfail)
    call check_val("eps_snow",        c%eps_snow,       0.98_wp,        nfail)
    call check_val("sigma_sb",        c%sigma_sb,       5.670373e-8_wp, nfail)
    call check_val("T0",              c%T0,             273.15_wp,      nfail)
    call check_val("seconds_per_day", c%seconds_per_day,86400.0_wp,     nfail)

    call check("default albedo_scheme = dynamic", &
               c%albedo_scheme .eq. CHION_ALBEDO_DYNAMIC, nfail)
    call check("default fresh_snow_density_scheme = constant", &
               c%fresh_snow_density_scheme .eq. CHION_FRESH_SNOW_DENSITY_CONSTANT, nfail)
    call check("default low_density_densification = bessi", &
               c%low_density_densification .eq. CHION_DENSIFY_BESSI, nfail)

    ! --- Scheme name mapping, including the Chion.jl aliases -------------
    call check("albedo 'constant'",   chion_albedo_scheme_flag("constant")   .eq. CHION_ALBEDO_CONSTANT,   nfail)
    call check("albedo 'dynamic'",    chion_albedo_scheme_flag("dynamic")    .eq. CHION_ALBEDO_DYNAMIC,    nfail)
    call check("albedo 'prescribed'", chion_albedo_scheme_flag("prescribed") .eq. CHION_ALBEDO_PRESCRIBED, nfail)
    call check("albedo alias 'bessi'  -> constant", &
               chion_albedo_scheme_flag("bessi")  .eq. CHION_ALBEDO_CONSTANT, nfail)
    call check("albedo alias 'legacy' -> constant", &
               chion_albedo_scheme_flag("legacy") .eq. CHION_ALBEDO_CONSTANT, nfail)

    call check("fresh snow 'constant'", &
               chion_fresh_snow_density_scheme_flag("constant") .eq. CHION_FRESH_SNOW_DENSITY_CONSTANT, nfail)
    call check("fresh snow alias 'htessel' -> parameterized", &
               chion_fresh_snow_density_scheme_flag("htessel") .eq. CHION_FRESH_SNOW_DENSITY_PARAMETERIZED, nfail)

    call check("densify 'bessi'",   chion_densify_scheme_flag("bessi")   .eq. CHION_DENSIFY_BESSI,   nfail)
    call check("densify 'htessel'", chion_densify_scheme_flag("htessel") .eq. CHION_DENSIFY_HTESSEL, nfail)

    ! --- Enum validation --------------------------------------------------
    call chion_check_enum("chion","model","bessi","bessi|pdd|itm")
    call chion_check_enum("chion","model","itm",  "bessi|pdd|itm")
    call chion_check_enum("chion","model","pdd",  "bessi|pdd|itm")
    write(*,"(a)") "  ok   : chion_check_enum accepts all three model names"

    ! --- Forcing container ------------------------------------------------
    call chion_forcing_alloc(forc,ncol)

    call check("forcing ncol",              forc%ncol .eq. ncol,                       nfail)
    call check("forcing arrays allocated",  allocated(forc%air_temperature),           nfail)
    call check("forcing size",              size(forc%snowfall_rate) .eq. ncol,        nfail)
    call check("no prescribed fluxes",      .not. any(forc%has_q_sw_net),              nfail)
    call check("default pressure sea level", &
               all(abs(forc%air_pressure - DEF_SEA_LEVEL_AIR_PRESSURE) .lt. 1.0e-9_wp), nfail)

    ! --- Grid, no spatial coordinates ------------------------------------
    call chion_grid_init(grd,ncol)

    call check("grid ncol",         grd%ncol .eq. ncol,      nfail)
    call check("all active",        grd%n_active .eq. ncol,  nfail)
    call check("no spatial coords", .not. grd%has_spatial,   nfail)

    ! --- Active-mask packing ---------------------------------------------
    active    = .FALSE.
    active(2) = .TRUE.
    active(5) = .TRUE.
    active(6) = .TRUE.
    call chion_grid_set_active(grd,active)

    call check("n_active = 3",   grd%n_active .eq. 3,       nfail)
    call check("active_idx(1)=2", grd%active_idx(1) .eq. 2, nfail)
    call check("active_idx(2)=5", grd%active_idx(2) .eq. 5, nfail)
    call check("active_idx(3)=6", grd%active_idx(3) .eq. 6, nfail)

    active = .FALSE.
    call chion_grid_set_active(grd,active)
    call check("n_active = 0 when all inactive", grd%n_active .eq. 0, nfail)

    ! --- Step forcing type is usable --------------------------------------
    fc%air_temperature = 263.15_wp
    fc%dt_days         = 1.0_wp
    fc%snowfall_rate   = 1.0e-5_wp
    fc%has_q_sw_net    = .FALSE.
    call check("step forcing assignable", &
               abs(fc%air_temperature - 263.15_wp) .lt. 1.0e-12_wp .and. .not. fc%has_q_sw_net, nfail)

    ! --- Cleanup ----------------------------------------------------------
    call chion_forcing_dealloc(forc)
    call chion_grid_dealloc(grd)
    call check("forcing deallocated", .not. allocated(forc%air_temperature), nfail)
    call check("grid deallocated",    .not. allocated(grd%active),           nfail)

    write(*,*)
    call chion_const_print(c)
    write(*,*)
    write(*,"(a)") "=========================================================="
    if (nfail .eq. 0) then
        write(*,"(a)") " WP1: ALL CHECKS PASSED"
        write(*,"(a)") "=========================================================="
    else
        write(*,"(a,i0,a)") " WP1: ", nfail, " CHECK(S) FAILED"
        write(*,"(a)") "=========================================================="
        stop 1
    end if

contains

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

        ! Relative tolerance at sp resolution. A dp-style 1e-14 would fail
        ! every check now that wp = sp.
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

end program test_defs
