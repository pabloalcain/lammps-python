#ifndef SSF_H
#define SSF_H

#include "math.h"
#include "stdio.h"
#ifdef __cplusplus
extern "C" {
#endif
  void ssf(double *x, int *type, int natoms, double size, int npoints, double *k, double *gr);
#ifdef __cplusplus
}
#endif
#endif
