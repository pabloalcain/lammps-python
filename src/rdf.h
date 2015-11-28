#ifndef RDF_H
#define RDF_H

#include "math.h"
#include "stdio.h"
#ifdef __cplusplus
extern "C" {
#endif
  void rdf(double *x, int *type, int natoms, int nbins, double size, double *gr);
#ifdef __cplusplus
}
#endif
double f1(double l);
double f2(double l);
double volume(double l);
#endif
