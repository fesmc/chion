###############################################
##
## Rules for individual chion files
##
###############################################
#
# One explicit rule per object, with hand-written dependencies. There is no
# auto-dependency generation and no wildcard, matching yelmo. When adding a
# source file, add its rule here AND add it to the appropriate object list at
# the bottom.
#
# Build order matters: chion_defs must come first, since every other module
# uses it. Object lists are ordered so that `ar` receives dependencies before
# dependents.

## chion base ##################################

# chion_defs is the ONLY preprocessed source. The .F90 extension makes both
# gfortran and ifort preprocess it without a -cpp/-fpp flag; the Makefile
# passes -DCHION_DP when precision=dp. See docs/porting_notes.md D19.
$(objdir)/chion_defs.o: $(srcdir)/chion_defs.F90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# Monthly -> daily forcing interface. Model-neutral, depends only on chion_defs.
$(objdir)/chion_forcing_monthly.o: $(srcdir)/chion_forcing_monthly.f90 \
						  	$(objdir)/chion_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

## chion physics ###############################
#
# Single-column kernels. Each operates on contiguous column slices (see
# docs/porting_notes.md D8) and depends only on chion_defs plus, where noted,
# the layer utilities.

$(objdir)/snow_column_utils.o: $(physdir)/snow_column_utils.f90 \
						  	$(objdir)/chion_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_layers.o: $(physdir)/snow_layers.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# snow_surface_fluxes uses snow_layers for the depleted-surface removal and
# surface-merge loops inside apply_snow_surface_vapor_mass_flux.
$(objdir)/snow_surface_fluxes.o: $(physdir)/snow_surface_fluxes.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o \
						  	$(objdir)/snow_layers.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_energy.o: $(physdir)/snow_energy.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o \
						  	$(objdir)/snow_surface_fluxes.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_percolation.o: $(physdir)/snow_percolation.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_refreezing.o: $(physdir)/snow_refreezing.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_melt.o: $(physdir)/snow_melt.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o \
						  	$(objdir)/snow_layers.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_albedo.o: $(physdir)/snow_albedo.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_albedo_semix.o: $(physdir)/snow_albedo_semix.f90 \
						  	$(objdir)/chion_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_densify.o: $(physdir)/snow_densify.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_diurnal.o: $(physdir)/snow_diurnal.f90 \
						  	$(objdir)/chion_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_accumulation.o: $(physdir)/snow_accumulation.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o \
						  	$(objdir)/snow_layers.o $(objdir)/snow_albedo.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_diagnostics.o: $(physdir)/snow_diagnostics.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_pdd.o: $(physdir)/snow_pdd.f90 \
						  	$(objdir)/chion_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_bessi.o: $(physdir)/snow_bessi.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o \
						  	$(objdir)/snow_layers.o $(objdir)/snow_albedo.o \
						  	$(objdir)/snow_albedo_semix.o \
						  	$(objdir)/snow_accumulation.o $(objdir)/snow_densify.o \
						  	$(objdir)/snow_diurnal.o $(objdir)/snow_surface_fluxes.o \
						  	$(objdir)/snow_energy.o $(objdir)/snow_melt.o \
						  	$(objdir)/snow_percolation.o $(objdir)/snow_refreezing.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/snow_itm.o: $(physdir)/snow_itm.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

## chion core ##################################
#
# Dispatcher, public API and facade. These sit above the physics modules and
# must be archived after them.

$(objdir)/chion_model.o: $(srcdir)/chion_model.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_bessi.o \
						  	$(objdir)/snow_pdd.o $(objdir)/snow_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/chion_api.o: $(srcdir)/chion_api.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/chion_model.o \
						  	$(objdir)/snow_bessi.o $(objdir)/snow_pdd.o $(objdir)/snow_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/chion_io.o: $(srcdir)/chion_io.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/chion_api.o \
						  	$(objdir)/snow_diagnostics.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

$(objdir)/chion.o: $(srcdir)/chion.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/chion_model.o $(objdir)/chion_api.o \
						  	$(objdir)/snow_bessi.o $(objdir)/snow_pdd.o $(objdir)/snow_itm.o \
						  	$(objdir)/snow_diagnostics.o $(objdir)/chion_forcing_monthly.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

## driver-layer modules ########################
#
# NOT part of libchion.a. These sit ABOVE the library, in the host/driver
# layer, and are compiled straight into a driver (see the `grid` target).
# Keeping them out of the library preserves chion's contract that the host
# supplies forcing and insolation -- the library has no insolation module and
# no knowledge of where raw datasets live.
#
# insolation is vendored under libs/insol (top-of-atmosphere daily insolation,
# Laskar LA2004 tables in libs/insol/input). It uses fesm-utils' interp1D for
# the spline, so it needs INC_FESMUTILS like everything else.

$(objdir)/insolation.o: $(libsdir)/insol/insolation.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

# chion_domain assembles standardized forcing from raw ice_data. It sits on
# libchion.a (chion_defs), fesm-utils' coords/ncio, and the vendored insolation.
$(objdir)/chion_domain.o: $(libsdir)/domains/chion_domain.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

###############################################
##
## Object lists
##
###############################################

chion_base =    $(objdir)/chion_defs.o \
                $(objdir)/chion_forcing_monthly.o

chion_physics = $(objdir)/snow_column_utils.o \
				$(objdir)/snow_layers.o \
				$(objdir)/snow_surface_fluxes.o \
				$(objdir)/snow_energy.o \
				$(objdir)/snow_percolation.o \
				$(objdir)/snow_refreezing.o \
				$(objdir)/snow_melt.o \
				$(objdir)/snow_albedo.o \
				$(objdir)/snow_albedo_semix.o \
				$(objdir)/snow_densify.o \
				$(objdir)/snow_diurnal.o \
				$(objdir)/snow_accumulation.o \
				$(objdir)/snow_diagnostics.o \
				$(objdir)/snow_bessi.o \
				$(objdir)/snow_pdd.o \
				$(objdir)/snow_itm.o

chion_core =    $(objdir)/chion_model.o \
                $(objdir)/chion_api.o \
                $(objdir)/chion_io.o \
                $(objdir)/chion.o

# Driver-layer objects, linked into drivers that need domain loading and
# insolation (currently chion_grid.x). Grows with the domain modules (WP: see
# chion_domain, added alongside the domain loaders).
driver_objs =   $(objdir)/insolation.o \
                $(objdir)/chion_domain.o
