#ifndef CLUSTER_H
#define CLUSTER_H

#include "math.h"
#include "stdio.h"
#include "stdbool.h"
#include "stdlib.h"

#define MIN(a, b) \
  ({__typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#ifdef __cplusplus
extern "C" {
#endif
  void cluster(double *x, double *v, int *type, int natoms,
               double size, bool energy, int model, int *index);
  int connections(int *index, double *x, double *v, int *type, int natoms,
                  double size, double expansion, bool energy, int model,
                  int *connect);
  double enclus(double *x, double *v, int *type, int natoms,
                double size, double expansion, int model);

#ifdef __cplusplus
}
#endif
double potential(double r, int model);
double potentialij(double r, int t1, int t2, int model);
double distance(double *d, double size, int *idx);
static inline bool link(double rsq, double *dv, int idx, int t1, int t2, int model, double expansion);
#endif
