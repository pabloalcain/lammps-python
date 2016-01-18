#include "cluster.h"

void cluster(double *x, double *v, int *type, int natoms, double size, int *index){
  /* Calculates clusters. 

     On index, returns a 1D array with the cluster index*/
  double cutsq;
  cutsq = 5.4;
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
          if (d > size/2) d -= size;
          else if (d < -size/2) d += size;
          rsq += d * d;
        }
        if (rsq < cutsq) {
          index[i] = index[j] = MIN(index[i],index[j]);
          done = false;
        }
      }
    }
    if (done) break;
  }
  return;
}
