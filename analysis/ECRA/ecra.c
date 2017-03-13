#include "ecra.h"

double enpart(double *x, double *v, int *type, int natoms,
              double size, double expansion, int *index){
  /* Calculates the energy of a cluster partition given by index*/
  double cutsq; // insert it after testing
  double mass;
  int *uniq = (int *) malloc(natoms * sizeof(int));
  int nuniq = 0;
  cutsq = 5.4 * 5.4;
  mass = 938.0;

  for (int i = 0; i < natoms; i++) {
    bool new = true;
    for (int j = 0; j < nuniq; j++) {
      if (index[i] == uniq[j]) {
        new = false;
        break;
      }
    }
    if (new) {
      uniq[nuniq++] = index[i];
    }
  }
  double e = 0.0;
  for (int icl = 0; icl < nuniq; icl++) {
    double pot = 0.0;
    double kin = 0.0;
    int nclus = 0;
    for (int i = 0; i < natoms; i++) {
      if (index[i] != uniq[icl]) continue;
      nclus += 1;
      for (int j = i + 1; j < natoms; j++) {
        if (index[j] != uniq[icl]) continue;
        double dv = 0;
        double dx = 0;
        for (int k = 0; k < 3; k++) {
          double tdx = 0;
          double tdv = 0;
          tdx = x[3*i + k] - x[3*j + k];
          tdv = v[3*i + k] - v[3*j + k];
          if (tdx > size/2) {
            dx += (tdx - size) * (tdx - size);
            dv += (tdv + 2 * expansion * size) * (tdv + 2 * expansion * size);
          } else if (tdx <  - size/2) {
            dx += (tdx + size) * (tdx + size);
            dv += (tdv - 2 * expansion * size) * (tdv - 2 * expansion * size);
          } else {
            dx += tdx * tdx;
            dv += tdv * tdv;
          }
        }
        if (dx < cutsq) {
          dx = sqrt(dx);
          pot += potential(dx, type[i], type[j]);
        }
        kin += dv * mass/2;
      }
    }
    e += pot + kin/nclus;
  }
  free(uniq);
  return e;
}

double binary_fusion(double *x, double *v, int *type, int natoms,
              double size, double expansion, int *index){
  /* Calculates and produces the most favorable binary fusion*/
  /* Makes a very large malloc! */
  double cutsq; // insert it after testing
  double mass;
  cutsq = 5.4 * 5.4;
  mass = 938.0;
  int nthreads;
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  //double *minde_tid = (double *) malloc(nthreads * sizeof(double));
  int *i_merge_tid = (int *) malloc(nthreads * sizeof(int));
  int *j_merge_tid = (int *) malloc(nthreads * sizeof(int));
  double *minde_tid = (double *) malloc(nthreads * sizeof(double));
  for (int i = 0; i < nthreads; i++) {
    minde_tid[i] = 0;
    i_merge_tid[i] = 0;
    j_merge_tid[i] = 0;
  }
  int *uniq = (int *) malloc(natoms * sizeof(int));
  int nuniq = 0;

  for (int i = 0; i < natoms; i++) {
    bool new = true;
    for (int j = 0; j < nuniq; j++) {
      if (index[i] == uniq[j]) {
        new = false;
        break;
      }
    }
    if (new) {
      uniq[nuniq++] = index[i];
    }
  }

  int **clus = (int **) malloc(natoms * sizeof(int*));
  int *size_clus = (int *) malloc(natoms * sizeof(int));

  for (int i = 0; i < natoms; i++) {
    clus[i] = (int *) malloc(natoms * sizeof(int));
    size_clus[i] = 0;
  }

  for (int i = 0; i < natoms; i++) {
    int c = index[i];
    clus[c][size_clus[c]] = i;
    size_clus[c] +=1;
  }
  double minde = 0;
  int i_merge = 0;
  int j_merge = 0;
#pragma omp parallel
  {
    for (int iicl = 0; iicl < natoms; iicl+=nthreads){
      int tid = omp_get_thread_num();
      int icl = iicl + tid;
      if (size_clus[icl] == 0) continue;
      if (icl >= natoms) break;
      for (int jcl = icl + 1; jcl < natoms; jcl++) {
        if (size_clus[jcl] == 0) continue;
        double dpot = 0.0;
        double dvi[3];
        for (int k = 0; k < 3; k++) dvi[k] = 0;

        for (int ii = 0; ii < size_clus[icl]; ii++) {
          int i = clus[icl][ii];
          double dvj[3];
          for (int k = 0; k < 3; k++) {
            dvi[k] += v[3*i + k];
          }
          for (int k = 0; k < 3; k++) dvj[k] = 0;
          for (int jj = 0; jj < size_clus[jcl]; jj++) {
            int j = clus[jcl][jj];
            double dx = 0;
            if (i == 0) {
              for (int k = 0; k < 3; k++) {
                dvj[k] += v[3*j + k];
              }
            }
            for (int k = 0; k < 3; k++) {
              double tdx = x[3*i + k] - x[3*j + k];
              if (tdx > size/2) {
                dx += (tdx - size) * (tdx - size);
              } else if (tdx <  - size/2) {
                dx += (tdx + size) * (tdx + size);
              } else {
                dx += tdx * tdx;
              }
            }
            if (dx < cutsq) {
              dx = sqrt(dx);
              dpot += potential(dx, type[i], type[j]);
            }
          }
          double dv = 0.0;
          for (int k = 0; k < 3; k++) {
            dv += (dvi[k]/size_clus[icl] - dvj[k]/size_clus[jcl]) * (dvi[k]/size_clus[icl] - dvj[k]/size_clus[jcl]);
          }
          double de = dpot + size_clus[icl] * size_clus[jcl] / (size_clus[icl] + size_clus[jcl]) * mass / 2 * dv;
          if (de < minde_tid[tid]) {
            minde_tid[tid] = de;
            i_merge_tid[tid] = icl;
            j_merge_tid[tid] = jcl;
          }
        }
      }
    }
  }

  for (int i = 0; i < nthreads; i++) {
    if (minde_tid[i] < minde) {
      minde = minde_tid[i];
      i_merge = i_merge_tid[i];
      j_merge = j_merge_tid[i];
    }
  }
  if (minde < 0) {
    for (int i = 0; i < natoms; i++){
      if (index[i] == i_merge) index[i] = j_merge;
    }
  }
  free(j_merge_tid);
  for (int i = 0; i < natoms; i++) {
    free(clus[i]);
  }
  free(uniq);
  free(i_merge_tid);
  free(minde_tid);
  free(clus);
  free(size_clus);
  return minde;
}

double potential(double r, int t1, int t2) {
  /*  double Vr = 3088.118;
  double Va = 2666.647;
  double V0 = 373.118;
  double Vc = 1.44;

  double ur = 1.7468;
  double ua = 1.6;
  double u0 = 1.5;
  double uc = 0.05;*/


  double Vr = 3097.0;
  double Va = 2696.0;
  double V0 = 379.5;
  double Vc = 1.44;

  double ur = 1.648;
  double ua = 1.528;
  double u0 = 1.628;
  double uc = 0.05;
  if (t1 == t2 && t1 == 1) {
    return V0 * exp(-u0 * r) / r;
  } else if (t1 == t2 && t1 == 2) {
    return V0 * exp(-u0 * r) / r + Vc * exp(-uc * r) / r;
  } else {
    return Vr * exp(-ur * r) / r - Va * exp(-ua * r) / r;
  }
}
