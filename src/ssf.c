#include "ssf.h"

void ssf(double *x, int *type, int natoms, double size, int npoints, int nrep, double *k, double *sk){
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
  
  /* 2*nrep is the number of actual replicas we calculate, 2*(2*nrep)
     + 1 is the possible distances among all the replicas */
  int ndist = 4*nrep+1;
  int *deg = (int *) malloc((ndist)*sizeof(int));
  int *dist = (int *) malloc((ndist)*sizeof(int));
  for (int ii = 0; ii < ndist; ii++) {
    dist[ii] = ii - 2*nrep;
    deg[ii] = 2*nrep + 1 - abs(dist[ii]);
  }
  
  #ifdef _OPENMP 
  #pragma omp parallel
  #endif
  {
    int tid = 0;
    int nthreads = 1;
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    nthreads = omp_get_num_threads();
    #endif
    for (int it = 0; it < npoints; it += nthreads) {
      int ii = it + tid;
      if (tid == 0) {
        printf(" %d%%...\r", (ii * 100)/npoints);
        fflush(stdout);
      }
      if (ii >= npoints) break;
      double ki = k[ii];
      sk[ncols * ii] = ki;
      for (int i = 0; i < natoms; i++) {
        double x1[3];
        double dx, r;
        
        int itype = type[i];
        for (int l = 0; l < 3; l++) x1[l] = x[3*i + l];
        for (int j = 0*(i+1); j < natoms; j++) {
          int jtype = type[j];
          int idxpair = labels[itype][jtype];
          for (int ix = 0; ix < ndist; ix ++) {
            for (int iy = 0; iy < ndist; iy ++) {
              for (int iz = 0; iz < ndist; iz ++) {
                int deg_tot = deg[ix] * deg[iy] * deg[iz];
                int this_dist[3] = {dist[ix], dist[iy], dist[iz]};
                double cont;
                r = 0.0;
                for (int l = 0; l < 3; l++) {
                  dx = x[3*j + l] - x1[l] + size * this_dist[l];
                  r += dx * dx;
                }
                r = sqrt(r);
                if (ki*r == 0.0d) cont = 1.0d;
                else cont = sin(ki*r)/(ki*r);
                cont = cont * deg_tot;
                sk[ncols * ii + 1] += cont;
                sk[ncols * ii + idxpair] += cont;
                if (ii == 0) {
                  hist[1] += deg_tot;
                  hist[idxpair] += deg_tot;
                }
              }
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < npoints; i++) {
    for (int j = 1; j < ncols; j++) {
      sk[ncols * i + j] /= natoms * ( 2 * nrep + 1) * ( 2 * nrep + 1) * ( 2 * nrep + 1);
      sk[ncols * i + j] *= sqrt((double)hist[j]/hist[1]);
    }
  }
  printf("Done!      \n");

  free(deg);
  free(dist);
  return;
}
