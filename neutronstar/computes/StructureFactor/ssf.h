#ifndef SSF_H
#define SSF_H

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#include "lebedev.h"
#ifdef __cplusplus
extern "C" {
#endif
  void ssf(double *x, int *type, int natoms, double size, int npoints,
           int naver, int nrep, int *pairs, int npairs, double *k, double *gr);
#ifdef __cplusplus
}
#endif

#endif
