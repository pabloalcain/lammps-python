#include "ssf_angle.h"

void ssf_angle(double *x, int *type, int natoms, double size, int npoints,
               int naver, int nrep, int *pairs, int npairs, double *k, double *sk){
  /* Calculates ssf with in 3D for a cubic box and returns on ssf. x is
     an array that has the info in XYZ XYZ XYZ...fashion. The rdf is
     calculated for the whole box up to sqrt(3)/2.


     on sk, returns a 2D array with this format:

     column 0: k
     column k: s(k) in the pair labeled by k
  */
  int ncols = npairs + 1;
  for (int i = 0; i < ncols * npoints; i++)
    sk[i] = 0;

  // Warning!! maximum two types of particles
  int labels[3][3];
  for (int k = 0; k < npairs; k++) {
    int i = pairs[2*k];
    int j = pairs[2*k + 1];
    labels[i][j] = k + 1;
  }

  double *qx = (double *) malloc(naver*sizeof(double));
  double *qy = (double *) malloc(naver*sizeof(double));
  double *qz = (double *) malloc(naver*sizeof(double));
  double *qw = (double *) malloc(naver*sizeof(double));
  ld_by_order(naver, qx, qy, qz, qw);
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
        printf("\r%d%%...", (ii * 100)/npoints);
        fflush(stdout);
      }
      if (ii >= npoints) break;
      double s_cont = 0.0;
      double ki = k[ii];
      sk[ncols * ii] = ki;
      for (int j = 0; j < naver; j++) {
        double q1[3];
        double c_real, c_imag;
        q1[0] = ki * qx[j];
        q1[1] = ki * qy[j];
        q1[2] = ki * qz[j];
        c_real = 0.0;
        c_imag = 0.0;
        for (int i = 0; i < natoms; i++) {
          double cont;
          cont = 0.0;
          for (int l = 0; l < 3; l++) {
            cont += q1[l] * x[3*i+l];
          }
          c_real += cos(cont);
          c_imag += sin(cont);
        }
	
        double pbc_real = 0.0;
        double pbc_imag = 0.0;
        for (int i = 0; i < nrep; i++) {
          for (int l = 0; l < nrep; l++) {
            for (int m = 0; m < nrep; m++) {
              double cont;
              cont = q1[0] * size * i;
              cont += q1[1] * size * l;
              cont += q1[2] * size * m;
              pbc_real += cos(cont);
              pbc_imag += sin(cont);
            }
          }
        }
        s_cont += (c_real * c_real + c_imag * c_imag)*(pbc_real * pbc_real + pbc_imag * pbc_imag)*qw[j];
      }
      sk[ncols * ii + 1] = s_cont;
    }
  }
  for (int i = 0; i < npoints; i++) {
    for (int j = 1; j < ncols; j++) {
      sk[ncols * i + j] /= natoms*nrep*nrep*nrep;
    }
  }
  printf("\rDone!      \n");
  free(qx);
  free(qy);
  free(qz);
  free(qw);
  return;
}
