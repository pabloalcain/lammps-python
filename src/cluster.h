#ifndef CLUSTER_H
#define CLUSTER_H

#include "math.h"
#include "stdio.h"
#include "stdbool.h"

#define MIN(a, b) \
  ({__typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

#ifdef __cplusplus
extern "C" {
#endif
  void cluster(double *x, double *v, int *type, int natoms,
               double size, bool energy, bool pbc, int *index);
  int connections(int *index, double *x, double *v, int *type, int natoms,
		  double size, double expansion, bool energy, int *connect);

#ifdef __cplusplus
}
#endif
double potential(double r);
#endif
