! chion_physics.f90 -- the physics modules as ONE translation unit.
!
! The modules stay one file each (snow_*.f90) for editing; this file only
! includes them, in dependency order, so that the compiler sees all of them at
! once and inlines across modules (small helpers such as surface_has_snow,
! safe_positive or the turbulence called per layer or per substep), as -ipo
! would, but into an ordinary object: libchion.a needs no link-time
! optimization from its hosts. The module names and .mod files are unchanged.
! config/Makefile_chion.mk compiles this file into chion_physics.o; when
! adding a physics module, include it here after the modules it uses.

include "snow_column_utils.f90"
include "snow_layers.f90"
include "snow_vapor.f90"
include "snow_seb_semix.f90"
include "snow_turbulence.f90"
include "snow_diurnal.f90"
include "snow_surface_fluxes.f90"
include "snow_energy.f90"
include "snow_percolation.f90"
include "snow_refreezing.f90"
include "snow_melt.f90"
include "snow_albedo.f90"
include "snow_albedo_semix.f90"
include "snow_densify.f90"
include "snow_accumulation.f90"
include "snow_diagnostics.f90"
include "snow_bessi.f90"
include "snow_pdd.f90"
include "snow_itm.f90"
