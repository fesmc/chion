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

# Dependency paths. OpenMP is the DEFAULT build (openmp ?= 1); chion
# parallelizes over snowpack columns, which are fully independent, so results
# are identical to serial. Runtime thread count is OMP_NUM_THREADS. Build serial
# with `make openmp=0`.
FESMUTILSROOT = fesm-utils
ifeq ($(openmp), 0)
    FESMUTILSLIBDIR = ${FESMUTILSROOT}/include-serial
else
    FESMUTILSLIBDIR = ${FESMUTILSROOT}/include-omp

    FFLAGS += $(FFLAGS_OPENMP)
endif
INC_FESMUTILS = -I${FESMUTILSLIBDIR}
LIB_FESMUTILS = -L${FESMUTILSLIBDIR} -lfesmutils

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LFLAGS_EXTRA)
