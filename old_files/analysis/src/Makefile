# -*- makefile -*-
.PHONY: all clean

F2PY = f2py 
#F2PY_FLAGS = --quiet
F2PY_FLAGS +=-DF2PY_REPORT_ON_ARRAY_COPY=1 
F2PY_FLAGS +=--noarch --f90flags='${ARCHFLAGS} -fno-range-check' --fcompiler=gnu95

GOFR_S=gofr.f90 gofr_sub.f90 gofr_mod.f90
GOFR_M=gofr

MSTE_S=mste.f90 mste_sub.f90 mste_mod.f90
MSTE_M=mste

MINK_S=minkowski.f90 minkowski_sub.f90 minkowski_mod.f90
MINK_M=minkowski

LIND_S=lindemann.f90 lindemann_sub.f90 lindemann_mod.f90
LIND_M=lindemann

all: ${MINK_M} ${MSTE_M} ${GOFR_M} ${LIND_M}

${MINK_M}: ${MINK_S}
	${F2PY} ${F2PY_FLAGS} -c ${MINK_S} -m ${MINK_M}

${GOFR_M}: ${GOFR_S}
	${F2PY} ${F2PY_FLAGS} -c ${GOFR_S} -m ${GOFR_M}

${MSTE_M}: ${MSTE_S}
	${F2PY} ${F2PY_FLAGS} -c ${MSTE_S} -m ${MSTE_M}

${LIND_M}: ${LIND_S}
	${F2PY} ${F2PY_FLAGS} -c ${LIND_S} -m ${LIND_M}

clean:
	${RM} ${MINK_M}.so ${MSTE_M}.so ${GOFR_M}.so ${LIND_M}.so
