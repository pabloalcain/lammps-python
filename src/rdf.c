#include "rdf.h"

void rdf(double *x, int *type, int natoms, int nbins, double size, double *gr){
  /* Calculates rdf with PBC in 3D for a cubic box and returns on
     gr. x is an array that has the info in XYZ XYZ XYZ...fashion. The
     rdf is calculated for the whole box up to sqrt(3)/2.

     So far it assumes we have 2 different type of particles, and
     labels:

     k = 1 => all v all
     k = 2 =>  1  v  1
     k = 3 =>  1  v  2
     k = 4 =>  2  v  2

     on gr, returns a 2D array with this format:

     column 0: Distance
     column k: g(r) in the pair labeled by k
*/
  double rmax;
  double dr;
  double x1[3];
  double dx, r;
  int idx;

  int npairs = 4;
  int ncols = npairs + 1;
  int hist[ncols];
  int labels[3][3];
  labels[1][1] = 2;
  labels[1][2] = 3;
  labels[2][1] = 3;
  labels[2][2] = 4;

  for (int k = 0; k < ncols * nbins; k++) gr[k] = 0.0;
  for (int k = 0; k < ncols; k++) hist[k] = 0;

  rmax = size * sqrt(3.0)/2;
  dr = rmax/nbins;
  for (int i = 0; i < natoms; i++) {
    int itype = type[i];
    for (int k = 0; k < 3; k++) x1[k] = x[3*i + k];
    for (int j = i+1; j < natoms; j++) {
      int jtype = type[j];
      int idxpair = labels[itype][jtype];
      r = 0.0;
      for (int k = 0; k < 3; k++) {
        dx = x[3*j + k] - x1[k];
        while (dx > size/2)  dx -= size;
        while (dx < -size/2) dx += size;
        r += dx * dx;
      }
      r = sqrt(r);
      if (r > rmax) continue;
      idx = (int)(r*nbins/rmax);
      gr[ncols * idx + 1] += 1.0;
      hist[1] += 1;
      gr[ncols * idx + idxpair] += 1.0;
      hist[idxpair] += 1;
    }
  }
  for (int k = 0; k < nbins; k++) {
    double p;
    r = (k+0.5)*rmax/nbins;
    gr[ncols * k] = r;
    p = volume((r+0.5*dr)/size) - volume((r-0.5*dr)/size);
    for (int i = 1; i < ncols; i++) {
      gr[ncols * k + i] /= p * hist[i];
    }
  }
  return;
}

double f1(double l) {
  return atan(sqrt(4 * l * l - 2));
}

double f2(double l) {
  double ang = 2 * l * (4 * l * l - 3) / (sqrt(4 * l * l - 2) * (4 * l * l + 1));
  return 8 * l * atan(ang);
}

double volume(double l) {
  /* Given the radius gives the volume that is inside the cube of
     length 1 AND the sphere of radius l. See for instance
     https://www.cmu.edu/biolphys/deserno/pdf/gr_periodic.pdf
  */

  double v;
  if (l < 0.5d) {
    v = 4.0/3.0 * M_PI * l * l * l;
  } else if (l < 1.0/sqrt(2.0)) {
    v = - M_PI/12 * (3 - 36 * l * l +32 * l * l * l);
  }
  else if (l < 0.5*sqrt(3.0d)) {
    v = - M_PI/4 + 3* M_PI * l * l;
    v += sqrt(4* l * l - 2);
    v += (1 - 12 * l * l) * f1(l);
    v += 2.0/3.0 * l * l * f2(l);
  }
  else {
    fprintf(stderr, "Warning! l = %g; max value is %g\n", l, 0.5*sqrt(3.0d));
    v = 1;
  }
  return v;
}
