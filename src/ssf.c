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
  int nr = 2;
  int ndist = 2*nr+1;
  int *deg = (int *)malloc((ndist)*sizeof(int));
  int *dist = (int *)malloc((ndist)*sizeof(int));
  for (int ii = 0; ii < ndist; ii++)
    dist[ii] = ii - 2;
  deg[0] = 1;
  deg[1] = 2;
  deg[2] = 3;
  deg[3] = 2;
  deg[4] = 1;

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
      if (ii >= npoints) break;
      double ki = k[ii];
      sk[ncols * ii] = ki;
      printf("%d/%d\n", ii, npoints);
      for (int i = 0; i < natoms; i++) {
	int itype = type[i];
	for (int l = 0; l < 3; l++) x1[l] = x[3*i + l];
	for (int j = i+1; j < natoms; j++) {
	  int jtype = type[j];
	  int idxpair = labels[itype][jtype];
	  for (int ix = 0; ix < ndist; ix ++) {
	    for (int iy = 0; iy < ndist; iy ++) {
	      for (int iz = 0; iz < ndist; iz ++) {
		double deg_tot = deg[ix] * deg[iy] * deg[iz];
		double cont;
		r = 0.0;
		for (int l = 0; l < 3; l++) {
		  dx = x[3*j + l] - x1[l] + size * dist[l];
		  r += dx * dx;
		}
		r = sqrt(r);
		if (ki == 0.0d) cont = 1.0d;
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
      sk[ncols * i + j] /= natoms;
      sk[ncols * i + j] *= sqrt((double)hist[j]/hist[1]);
      sk[ncols * i + j] += 1;
    }
  }
  
  return;
}
