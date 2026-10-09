.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Command-line options at make call
debug        ?= 0
openmp       ?= 1
precision    ?= sp
legacy_chion ?= 0
fpsafe       ?= 0

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

# fpsafe=1 builds with value-safe optimization (-O2, no fast-math) instead of
# the machine fragment's DFLAGS_NODEBUG, in its own directories, so it never
# overwrites the production objects. The acceptance tests and validation/ use
# it: they check results to one ulp, which value-unsafe optimization (ifx
# -Ofast on Levante) does not preserve, and ifx -Ofast miscompiles
# tests/test_layers.f90. Production builds (yelmox, chion_grid runs) keep the
# machine's flags.
ifeq ($(fpsafe), 1)
	psuffix := $(psuffix)-fpsafe
else ifneq ($(fpsafe), 0)
$(error fpsafe must be 0 or 1, got '$(fpsafe)')
endif

# PATH options
#
# The sp build keeps the unsuffixed paths, so the artifact locations named by
# .configme/manifest.toml and .runme/info.json are unchanged and the dp build
# is purely additive. The two builds therefore coexist without `make clean`.
#
# openmp=0 and openmp=1 share these directories (as in yelmo; yelmox links
# libchion/include for both). The objects instead depend on a stamp naming the
# OpenMP setting they were built with ($(ompstamp), below), so switching
# openmp rebuilds them: `configme install` (serial, then omp) leaves an OpenMP
# library, and a host's `make openmp=0|1` always gets the variant it asked for.
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

# Value-safe optimization per compiler (fpsafe=1): the compiler fragments' own
# portable defaults, which the machine fragments override with -Ofast.
DFLAGS_FPSAFE_gfortran = -O2
DFLAGS_FPSAFE_ifx      = -O2 -fp-model precise
DFLAGS_FPSAFE_ifort    = -O2 -fp-model precise

ifeq ($(fpsafe), 1)
    DFLAGS_NODEBUG := $(DFLAGS_FPSAFE_$(COMPILER))
    ifeq ($(DFLAGS_NODEBUG),)
        $(error fpsafe=1: no value-safe flags for COMPILER='$(COMPILER)')
    endif
endif

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

# Each build flavour (precision, legacy_chion, fpsafe) has its own object and
# binary directories, created on demand. Order-only prerequisites: a directory's
# mtime changes whenever a file in it does, so it must not trigger rebuilds.
$(objdir) $(bindir):
	mkdir -p $@

$(chion_base) $(chion_physics) $(chion_core) $(driver_objs): | $(objdir)

# The OpenMP setting the objects in $(objdir) were compiled with. A missing
# stamp (first build, or openmp switched) is created newer than every object,
# so all of them are rebuilt with the current flags; the other setting's stamp
# is removed so that switching back rebuilds again.
ompstamp = $(objdir)/.openmp$(openmp)

$(ompstamp): | $(objdir)
	rm -f $(objdir)/.openmp*
	touch $@

$(chion_base) $(chion_physics) $(chion_core) $(driver_objs): $(ompstamp)

# chion objects embed fesm-utils types: rebuild them when libfesmutils.a
# changes, or objects keep a stale type layout (fesmc/FastHydrology#10).
$(chion_base) $(chion_physics) $(chion_core): $(FESMUTILSLIBDIR)/libfesmutils.a

# Static library compilation. The archive is a file target, so its mtime only
# changes when an object does (hosts such as yelmox depend on it to rebuild
# their own objects).
chion-static: $(objdir)/libchion.a

$(objdir)/libchion.a: $(chion_base) $(chion_physics) $(chion_core)
	ar rc $@ $(chion_base) $(chion_physics) $(chion_core)
	ranlib $@
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

seb : chion-static
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -o $(bindir)/test_seb.x $(testdir)/test_seb.f90 \
			-L${CURDIR}/$(objdir) -lchion $(LFLAGS)
		@echo "    test_seb.x is ready."

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
test_targets = defs utils layers energy water surface seb wp7 bessi pdd diagnostics itm api io monthly
tests : $(test_targets)

# Driver executables
driver_targets = column grid
drivers : $(driver_targets)

$(test_targets) $(driver_targets) domain : | $(bindir)

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
	@echo " make drivers      : compiles chion_column.x and chion_grid.x."
	@echo " make all          : tests + drivers."
	@echo " make clean        : cleans object files"
	@echo ""
	@echo " Options: debug=1 (bounds/traps), debug=2 (profiling), openmp=0 (serial; omp is default)"
	@echo "          precision=dp (libchion/*-dp), legacy_chion=1 (*-legacy; validation only)"
	@echo "          fpsafe=1 (*-fpsafe): value-safe -O2 instead of the machine's -Ofast."
	@echo "          Use it for the acceptance tests and validation/:"
	@echo "            make all fpsafe=1 [precision=dp [legacy_chion=1]]"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so $(objdir)/.openmp*
	rm -rf *.x.dSYM
	rm -f ./git_chion.txt
