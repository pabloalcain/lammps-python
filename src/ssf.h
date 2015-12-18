#ifndef SSF_H
#define SSF_H

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif
  void ssf(double *x, int *type, int natoms, double size, int npoints, int nrep, double *k, double *gr);
#ifdef __cplusplus
}
#endif
#endif
