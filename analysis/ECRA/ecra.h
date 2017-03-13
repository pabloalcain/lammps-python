#ifndef ECRA_H
#define ECRA_H

#include "math.h"
#include "stdio.h"
#include "stdbool.h"
#include "stdlib.h"
#include "omp.h"

#define MIN(a, b) \
  ({__typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#ifdef __cplusplus
extern "C" {
#endif
  double enpart(double *x, double *v, int *type, int natoms,
                double size, double expansion, int *index);
#ifdef __cplusplus
  double binary_fusion(double *x, double *v, int *type, int natoms,
                       double size, double expansion, int *index);
}
#endif
double potential(double r, int t1, int t2);
#endif
