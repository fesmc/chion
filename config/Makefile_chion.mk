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

$(objdir)/chion_defs.o: $(srcdir)/chion_defs.f90
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

$(objdir)/snow_surface_fluxes.o: $(physdir)/snow_surface_fluxes.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
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

$(objdir)/snow_itm.o: $(physdir)/snow_itm.f90 \
						  	$(objdir)/chion_defs.o $(objdir)/snow_column_utils.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

###############################################
##
## Object lists
##
###############################################

chion_base =    $(objdir)/chion_defs.o

chion_physics = $(objdir)/snow_column_utils.o \
				$(objdir)/snow_layers.o \
				$(objdir)/snow_surface_fluxes.o \
				$(objdir)/snow_energy.o \
				$(objdir)/snow_percolation.o \
				$(objdir)/snow_refreezing.o \
				$(objdir)/snow_melt.o \
				$(objdir)/snow_albedo.o \
				$(objdir)/snow_densify.o \
				$(objdir)/snow_diurnal.o \
				$(objdir)/snow_accumulation.o \
				$(objdir)/snow_diagnostics.o \
				$(objdir)/snow_pdd.o \
				$(objdir)/snow_itm.o
