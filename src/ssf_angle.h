#ifndef SSF_ANGLE_H
#define SSF_ANGLE_H

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif
  void ssf_angle(double *x, int *type, int natoms, double size, int npoints, int naver, int nrep, double *k, double *gr);
#ifdef __cplusplus
}
void random_rotation(double u1, double u2, double u3, double *i, double *vector);


#endif
#endif
