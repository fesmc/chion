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
# uses it.

## chion base ##################################

$(objdir)/chion_defs.o: $(srcdir)/chion_defs.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

## chion physics ###############################
#
# Single-column kernels. Each operates on plain arrays and depends only on
# chion_defs plus, where noted, the layer utilities. Populated by WP3-WP12.

$(objdir)/snow_column_utils.o: $(physdir)/snow_column_utils.f90 \
						  	$(objdir)/chion_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c -o $@ $<

###############################################
##
## Object lists
##
###############################################

chion_base =    $(objdir)/chion_defs.o

chion_physics = $(objdir)/snow_column_utils.o
