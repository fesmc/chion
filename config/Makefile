.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir  = src
physdir = src/physics
objdir  = libchion/include
bindir  = libchion/bin
testdir = tests

# Command-line options at make call
debug  ?= 0
openmp ?= 0

## COMPILER CONFIGURATION ##
# Assembled by `configme` from a compiler fragment, a machine fragment and an
# auto-detected netCDF block. The fragment defines the compiler (FC, FFLAGS,
# FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and any
# host-specific overrides (e.g. LFLAGS_EXTRA).
#
# Note the compiler fragment's FFLAGS_BASE references $(objdir), which is why
# the PATH options above must be defined before this point.

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC and LFLAGS_EXTRA set there. Defines the fesm-utils
# paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

###############################################
##
## List of chion rules and source files
##
###############################################

include config/Makefile_chion.mk

###############################################
##
## Compilation of complete programs
##
###############################################

check-hash:
	./check_githash > git_chion.txt

# Static library compilation
chion-static: $(chion_base) $(chion_physics)
	ar rc $(objdir)/libchion.a $(chion_base) $(chion_physics)
	ranlib $(objdir)/libchion.a
	./check_githash > git_chion.txt
	@echo " "
	@echo "    $(objdir)/libchion.a is ready."
	@echo " "

# Test and driver programs that use chion-static
defs : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_defs.x $(testdir)/test_defs.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo " "
		@echo "    test_defs.x is ready."
		@echo " "

utils : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_column_utils.x $(testdir)/test_column_utils.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo " "
		@echo "    test_column_utils.x is ready."
		@echo " "

layers : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_layers.x $(testdir)/test_layers.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_layers.x is ready."

energy : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_energy.x $(testdir)/test_energy.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_energy.x is ready."

water : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_water.x $(testdir)/test_water.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_water.x is ready."

surface : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_surface.x $(testdir)/test_surface.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_surface.x is ready."

pdd : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_pdd.x $(testdir)/test_pdd.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_pdd.x is ready."

diagnostics : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_diagnostics.x $(testdir)/test_diagnostics.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_diagnostics.x is ready."

itm : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_itm.x $(testdir)/test_itm.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_itm.x is ready."

# All acceptance tests
tests : defs utils layers energy water surface pdd diagnostics itm

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make chion-static : compiles libchion.a, the chion static library."
	@echo " make defs         : compiles test_defs.x, the chion_defs acceptance test."
	@echo " make utils        : compiles test_column_utils.x, the snow_column_utils test."
	@echo " make tests        : compiles all acceptance tests."
	@echo " make clean        : cleans object files"
	@echo ""
	@echo " Options: debug=1 (bounds/traps), debug=2 (profiling), openmp=1"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
	rm -f ./git_chion.txt
