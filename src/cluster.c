#include "cluster.h"

void cluster(double *x, double *v, int *type, int natoms,
             double size, bool energy, bool pbc, int *index){
  /* Calculates clusters.

     energy: if true, calculate MSTE
     pbc: if true, use PBC

     On index, returns a 1D array with the cluster index*/
  double cutsq;
  double mu;
  cutsq = 5.4 * 5.4;
  mu = 938.0/2;
  for (int i = 0; i < natoms; i++) index[i] = i;
  
  while (true) {
    bool done;
    done = true;
    for (int i = 0; i < natoms; i++) {
      double x1[3];
      double v1[3];
      int type1;
      for (int k = 0; k < 3; k++) x1[k] = x[3*i + k];
      for (int k = 0; k < 3; k++) v1[k] = v[3*i + k];
      type1 = type[i];
      for (int j = 0; j < natoms; j++) {
        double rsq;
        if (j == i) continue;
        if (index[i] == index[j]) continue;
        rsq = 0;
        for (int k = 0; k < 3; k++) {
          double d = x[3*j + k] - x1[k];
          if (pbc) {
            if (d > size/2) d -= size;
            else if (d < -size/2) d += size;
          }
          rsq += d * d;
        }
        if (rsq < cutsq) {
          bool clust;
          clust = true;
          if (energy) {
            if (type1 == type[j]) {
              clust = false;
            } else {
              double vsq, eng;
              vsq = 0.0;
              for (int k = 0; k < 3; k++) {
                double dv = v[3*j + k] - v1[k];
                vsq += dv * dv;
              }
              eng = potential(sqrt(rsq)) + 1.0/2.0 * mu * vsq;
              clust = (eng < 0.0);
            }
          }
          if (clust) {
            index[i] = index[j] = MIN(index[i],index[j]);
            done = false;
          }
        }
      }
    }
    if (done) break;
  }
  return;
}

double potential(double r) {
  double Vr = 3088.118;
  double Va = 2666.647;
  double ur = 1.7468;
  double ua = 1.6;
  return Vr * exp(-ur * r) / r - Va * exp(-ua * r) / r;
}

  

