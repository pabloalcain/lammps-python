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

int connections(int *index, double *x, double *v, int *type, int natoms,
                double size, double expansion, bool energy, int *connect) {
  /* from clusters labeled by index, return connections in connect,
     structured as [i1, j1, wall1, i2, j2, wall2, ...]
     
     for walls, these are 1/-1 for connections in first dimension,
     2/-2 for connections in second dimension and 4/-4 for connections
     in third dimension.

     returns the number of connections found.
  */
  int count = 0;
  double cutsq = 5.4 * 5.4;
  double mu = 938.0/2;
  static int direct[] = {1, 2, 4};
  for (int i = 0; i < natoms; i++) {
    double x1[3];
    double v1[3];
    int type1;
    for (int k = 0; k < 3; k++) x1[k] = x[3*i + k];
    for (int k = 0; k < 3; k++) v1[k] = v[3*i + k];
    type1 = type[i];
    for (int j = i + 1; j < natoms; j++) {
      int idx = 0;
      int c = 0;
      double rsq;
      rsq = 0;
      double min = size;
      for (int k = 0; k < 3; k++) {
        double d = x[3*j + k] - x1[k];
        if (d > size/2) {
          d -= size; /* d is < 0! */
          if (-d < min) {
            min = -d;
            idx = -direct[k];
            c++;
          }
        }
        else if (d < -size/2) {
          d += size;
          if (d < min) {
            min = d;
            idx = direct[k];
            c++;
          }
        }
        rsq += d * d;
      }
      /* if we connect through only one boundary! */
      bool clust = false;
      if (rsq < cutsq && c == 1) {
        clust = true;
        if (energy) {
          if (type[j] == type1) {
            clust = false;
          }
          else {
            double vsq = 0.0;
            for (int k = 0; k < 3; k++) {
              double dv = v[3*j + k] - v1[k];
	      if (idx == direct[k]) {
		dv += expansion;
	      }
	      else if (idx == -direct[k]){
		dv -= expansion;
	      }
              vsq += dv * dv;
            }
            double eng = potential(sqrt(rsq)) + 1.0/2.0 * mu * vsq;
            clust = (eng < 0.0);
          }
        }
      }
      if (clust) {
        bool add = true;
        for (int l = 0; l < count; l++) {
          if (connect[3*l] == index[i] &&
              connect[3*l + 1] == index[j] &&
              connect[3*l + 2] == idx) {
            add = false;
            break;
          }
        }
        if (add) {
          connect[3*count] = index[i];
          connect[3*count + 1] = index[j];
          connect[3*count + 2] = idx;
          connect[3*count + 3] = index[j];
          connect[3*count + 4] = index[i];
          connect[3*count + 5] = - idx;
          count += 2;
        }
      }
    }
  }
  return count;
}
  

double potential(double r) {
  double Vr = 3088.118;
  double Va = 2666.647;
  double ur = 1.7468;
  double ua = 1.6;
  return Vr * exp(-ur * r) / r - Va * exp(-ua * r) / r;
}
