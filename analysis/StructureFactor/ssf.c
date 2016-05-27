#include "ssf.h"

void ssf(double *x, int *type, int natoms, double size, int npoints,
         int naver, int nrep, int *pairs, int npairs, double *k, double *sk){
  /* Calculates ssf with in 3D for a cubic box and returns on ssf. x is
     an array that has the info in XYZ XYZ XYZ...fashion. The rdf is
     calculated for the whole box up to sqrt(3)/2.


     on sk, returns a 2D array with this format:

     column 0: k
     column k: s(k) in the pair labeled by k
  */
  int ncols = npairs + 1;
  int hist[ncols];
  
  for (int i = 0; i < ncols * npoints; i++)
    sk[i] = 0;

  for (int i = 0; i < ncols; i++)
    hist[i] = 0;

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
      if (ii >= npoints) break;
      double s_cont[ncols];
      for (int i = 0; i < ncols; i++) s_cont[i] = 0;
      double ki = k[ii];
      sk[ncols * ii] = ki;
      /* average in the sphere */
      for (int j = 0; j < naver; j++) {
        double q1[3];
        q1[0] = ki * qx[j];
        q1[1] = ki * qy[j];
        q1[2] = ki * qz[j];

        /* sum over all atom positions */
        double cell_real[ncols], cell_imag[ncols];
        for (int i =0; i < ncols; i++) {
          cell_real[i] = 0.0;
          cell_imag[i] = 0.0;
        }

        for (int i = 0; i < natoms; i++) {
          int itype = type[i];
          int idxpair = labels[itype][itype];
          if (ii == 0 && j == 0) {
            int all = labels[0][0];
            if (all != 0) hist[all] += nrep * nrep * nrep;
            hist[idxpair] += nrep * nrep * nrep;
          }
          double accum;
          accum = q1[0] * x[3*i + 0];
          accum += q1[1] * x[3*i + 1];
          accum += q1[2] * x[3*i + 2];
          double this_real, this_imag;
          this_real = cos(accum);
          this_imag = sin(accum);
          int all = labels[0][0];
          if (all != 0) {
            cell_real[all] += this_real;
            cell_imag[all] += this_imag;
          }
          cell_real[idxpair] += this_real;
          cell_imag[idxpair] += this_imag;
        }

        /* sum over all cell images */
        double pbc_real, pbc_imag;
        pbc_real = 0.0;
        pbc_imag = 0.0;
        for (int i = 0; i < nrep; i++) {
          for (int l = 0; l < nrep; l++) {
            for (int m = 0; m < nrep; m++) {
              double accum;
              accum = q1[0] * size * i;
              accum += q1[1] * size * l;
              accum += q1[2] * size * m;
              pbc_real += cos(accum);
              pbc_imag += sin(accum);
            }
          }
        }
        for (int i = 1; i < ncols; i++) {
          double contr = cell_real[i] * cell_real[i] + cell_imag[i] * cell_imag[i];
          contr *= pbc_real * pbc_real + pbc_imag * pbc_imag;
          s_cont[i] += contr * qw[j];
        }
      }
      for (int i = 1; i < ncols; i++) {
        sk[ncols * ii + i] = s_cont[i];
      }
    }
  }
  for (int i = 0; i < npoints; i++) {
    for (int j = 1; j < ncols; j++) {
      sk[ncols * i + j] /= hist[j];
    }
  }
  free(qx);
  free(qy);
  free(qz);
  free(qw);
  return;
}
