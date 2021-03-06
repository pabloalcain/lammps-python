#include "rdf.h"

void rdf(double *x, int *type, int natoms, int nbins, double size, int *pairs,
         int npairs, char pbc, double *gr) {
  /* Calculates rdf with PBC in 3D for a cubic box and returns on
     gr. x is an array that has the info in XYZ XYZ XYZ...fashion. The
     rdf is calculated for the whole box up to sqrt(3)/2.

     It gets the pairs at the array pairs, and puts the rdf of
     pair[2*k] and pair[2*k+1] in column k+1 of gr
  */
  double rmax;
  double dr;
  double x1[3];
  double dx, r;
  int idx;

  int ncols = npairs + 1;
  int hist[ncols];

  // Warning!! maximum two types of particles
  int labels[3][3];
  for (int i = 0; i < 3; i++) 
    for (int j = 0; j < 3; j++) 
      labels[i][j] = 0;

  for (int k = 0; k < npairs; k++) {
    int i = pairs[2*k];
    int j = pairs[2*k + 1];
    labels[i][j] = k + 1;
  }
  for (int k = 0; k < ncols * nbins; k++) gr[k] = 0.0;
  for (int k = 0; k < ncols; k++) hist[k] = 0;

  if (pbc) rmax = size * sqrt(3.0)/2;
  else rmax = size * sqrt(3.0);

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
        if (pbc) {
          while (dx > size/2)  dx -= size;
          while (dx < -size/2) dx += size;
        }
        r += dx * dx;
      }
      r = sqrt(r);
      idx = (int)(r*nbins/rmax);
      int all = labels[0][0];
      if (all != 0) {
        hist[all] += 1;
        gr[ncols * idx + all] += 1.0;
      }
      hist[idxpair] += 1;
      gr[ncols * idx + idxpair] += 1.0;
    }
  }
  for (int k = 0; k < nbins; k++) {
    double p;
    double r2, r1;
    r = (k+0.5)*rmax/nbins;
    r2 = (r+0.5*dr)/size;
    r1 = (r-0.5*dr)/size;
    gr[ncols * k] = r;
    if (pbc) 
      p = volume(r2) - volume(r1);
    else
      p = 4.0/3.0 * M_PI * (r2 * r2 * r2 - r1 * r1 * r1);
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
