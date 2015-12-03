#include "ssf.h"

void ssf(double *x, int *type, int natoms, double size, int npoints, double *k, double *sk){
  /* Calculates ssf with in 3D for a cubic box and returns on ssf. x is
     an array that has the info in XYZ XYZ XYZ...fashion. The rdf is
     calculated for the whole box up to sqrt(3)/2.

     So far it assumes we have 2 different type of particles, and
     labels:

     k = 1 => all v all
     k = 2 =>  1  v  1
     k = 3 =>  1  v  2
     k = 4 =>  2  v  2

     on sk, returns a 2D array with this format:

     column 0: k
     column k: s(k) in the pair labeled by k
*/
  double rmax, rmax2;
  double x1[3];
  double dx, r;

  int npairs = 4;
  int ncols = npairs + 1;
  int hist[ncols];
  int labels[3][3];
  labels[1][1] = 2;
  labels[1][2] = 3;
  labels[2][1] = 3;
  labels[2][2] = 4;
  for (int i = 0; i < ncols; i++)
    hist[i] = 0;
  for (int i = 0; i < ncols * npoints; i++)
    sk[i] = 0;
  FILE *fp = fopen("types.dat", "w");
  rmax = size/2;// * sqrt(3.0)/2;
  rmax2 = rmax * rmax;
  for (int ii = 0; ii < npoints; ii++) {
    double ki = k[ii];
    sk[ncols * ii] = ki;
    printf("%d/%d\n", ii, npoints);
    for (int i = 0; i < natoms; i++) {
      int itype = type[i];
      for (int l = 0; l < 3; l++) x1[l] = x[3*i + l];
      for (int j = i+1; j < natoms; j++) {
        int jtype = type[j];
        int idxpair = labels[itype][jtype];
        r = 0.0;
        for (int l = 0; l < 3; l++) {
          dx = x[3*j + l] - x1[l];
          while (dx > size/2)  dx -= size;
          while (dx < -size/2) dx += size;
          r += dx * dx;
        }
        //if (r > rmax2) continue;
        r = sqrt(r);
        sk[ncols * ii + 1] += sin(ki*r)/(ki*r); 
        sk[ncols * ii + idxpair] += sin(ki*r)/(ki*r);
        if (ii == 0) {
          hist[1] += 1;
          hist[idxpair] += 1;
        }
      }
    }
  }
  for (int i = 0; i < ncols; i++) {
    fprintf(fp, "%d, ", hist[i]);
  }
  fclose(fp);
  for (int i = 0; i < npoints; i++) {
    for (int j = 1; j < ncols; j++) {
      sk[ncols * i + j] /= natoms;
      sk[ncols * i + j] *= sqrt((double)hist[j]/hist[1]);
      sk[ncols * i + j] += 1;
    }
  }
  
  return;
}
