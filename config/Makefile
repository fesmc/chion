.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug        ?= 0
openmp       ?= 1
precision    ?= sp
legacy_chion ?= 0

# legacy_chion=1 reverts chion's deliberate physics corrections to Chion.jl's
# values, so validation/ can prove the port is faithful without every upstream
# bug chion fixes turning into a WP16 gate failure. NOT a production setting:
# it selects physics believed to be wrong. See src/chion_defs.F90.

# Working precision of the chion state and public interfaces (docs/PLAN.md
# section 3.1). `sp` is the production setting and the yelmo/yelmox interface
# kind; `dp` exists so that chion can be compared against Chion.jl, which is
# Float64 throughout, without sp round-off masquerading as a porting error --
# and so that the sp-vs-dp difference is a measured number rather than an
# assumption. See docs/porting_notes.md D19.
#
# The nine cumulative accumulators are wp_acc = dp in BOTH builds; this switch
# moves only `wp`. In a dp build wp and wp_acc are the same kind, which is
# harmless: chion declares no generic interfaces, so nothing becomes ambiguous.
ifeq ($(precision), dp)
	psuffix = -dp
else ifeq ($(precision), sp)
	psuffix =
else
$(error precision must be 'sp' or 'dp', got '$(precision)')
endif

ifeq ($(legacy_chion), 1)
	psuffix := $(psuffix)-legacy
endif

# PATH options
#
# The sp build keeps the unsuffixed paths, so the artifact locations named by
# .configme/manifest.toml and .runme/info.json are unchanged and the dp build
# is purely additive. The two builds therefore coexist without `make clean`.
#
# NOTE (pre-existing): openmp=0 and openmp=1 still share these directories, so
# switching openmp does need a clean. Only fesm-utils' include dir is swapped.
srcdir  = src
physdir = src/physics
libsdir = libs
objdir  = libchion/include$(psuffix)
bindir  = libchion/bin$(psuffix)
testdir = tests

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

# chion overrides the compiler fragment's DFLAGS_DEBUG to drop `underflow` from
# the FPE trap set (the fragment, and yelmo, trap it).
#
# Gradual underflow to zero is correct and expected behaviour throughout this
# model: exp() of a large negative argument appears in the HTESSEL thermal
# metamorphism term, the PDD normal-CDF tails, the albedo decay law and the
# densification Arrhenius factors. Trapping underflow turns correctly-rounded
# results into SIGILL, and made three of the ten acceptance tests abort under
# debug=1 while passing under -O2. invalid, zero and overflow are still trapped,
# which is where real bugs show up.
#
# Remove this override to restore yelmo's flag set.
DFLAGS_DEBUG = -w -g -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

# Select the working precision. src/chion_defs.F90 is the only preprocessed
# source; its .F90 extension makes both gfortran and ifort preprocess it
# automatically, so no -cpp/-fpp flag is needed and no configme compiler
# fragment has to know about this switch.
ifeq ($(precision), dp)
	DFLAGS += -DCHION_DP
endif

ifeq ($(legacy_chion), 1)
	DFLAGS += -DCHION_LEGACY
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
chion-static: $(chion_base) $(chion_physics) $(chion_core)
	ar rc $(objdir)/libchion.a $(chion_base) $(chion_physics) $(chion_core)
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

wp7 : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_wp7.x $(testdir)/test_wp7.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_wp7.x is ready."

bessi : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_bessi.x $(testdir)/test_bessi.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_bessi.x is ready."

api : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_api.x $(testdir)/test_api.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_api.x is ready."

io : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_io.x $(testdir)/test_io.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_io.x is ready."

monthly : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_monthly.x $(testdir)/test_monthly.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_monthly.x is ready."

domain : chion-static $(driver_objs)
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_domain.x $(testdir)/test_domain.f90 \
			$(driver_objs) -L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_domain.x is ready."

# Driver programs (note: not test_*.x -- these are the runme executables)
column : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/chion_column.x $(testdir)/chion_column.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    chion_column.x is ready."

grid : chion-static $(driver_objs)
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/chion_grid.x $(testdir)/chion_grid.f90 \
			$(driver_objs) -L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    chion_grid.x is ready."

itm : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_itm.x $(testdir)/test_itm.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_itm.x is ready."

# All acceptance tests
tests : defs utils layers energy water surface wp7 bessi pdd diagnostics itm api io monthly

# Driver executables
drivers : column grid

all : tests drivers

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
	@echo " Options: debug=1 (bounds/traps), debug=2 (profiling), openmp=0 (serial; omp is default)"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
	rm -f ./git_chion.txt
