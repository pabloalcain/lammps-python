#ifndef MINKOWSKI_C_H
#define MINKOWSKI_C_H

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
  void minkowski(double *x, int natoms, double size, double rad, double rcell, double *out);
#ifdef __cplusplus
}
#endif
void add_pixel(bool *lattice, int n, int i, int j, int k, double *out);
void set_pixel(bool *lattice, int n, int i, int j, int k, bool value);
bool get_pixel(bool *lattice, int n, int i, int j, int k);


const int v_body =  1; // view M-DR algorithm for explanation
const int s_body = -6;
const int s_face =  2;
const int c_body =  3;
const int c_face = -2;
const int c_edge =  1;
const int e_body = -1;
const int e_face =  1;
const int e_edge = -1;
const int e_vert =  1;

#endif
