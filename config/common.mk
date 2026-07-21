# Shared build configuration for chion (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). This file references
# variables those provide: FFLAGS / FFLAGS_OPENMP (compiler) and LIB_NC
# (machine or auto-detected netCDF).
#
# chion has exactly one Fortran dependency, fesm-utils, from which it uses
# precision, nml, ncio and variable_io. It needs no linear solver, no FFTW and
# no PETSc -- the only implicit solve in the model is a tridiagonal conduction
# system, handled in-tree by the Thomas algorithm.

# Dependency paths (serial build by default).
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

# OpenMP build (make openmp=1): swap the serial fesm-utils build for the
# OpenMP variant and append the compiler's OpenMP flag (FFLAGS_OPENMP, set in
# the compiler fragment). chion parallelizes over snowpack columns, which are
# fully independent, so this is the only parallelism in the package.
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    FFLAGS += $(FFLAGS_OPENMP)
endif

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LFLAGS_EXTRA)
