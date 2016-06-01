#include "cluster.h"

static int direct[] = {1, 2, 4, 8, 16, 32};
static int wallidx[] = {1, 2, 4, -1, -2, -4};

void cluster(double *x, double *v, int *type, int natoms,
             double size, bool energy, int *index){
  /* Calculates clusters.

     energy: if true, calculate MSTE

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
  for (int i = 0; i < natoms; i++) {
    double x1[3];
    double v1[3];
    for (int k = 0; k < 3; k++) x1[k] = x[3*i + k];
    for (int k = 0; k < 3; k++) v1[k] = v[3*i + k];
    for (int j = i + 1; j < natoms; j++) {
      int idx;
      double rsq, d[3], dv[3];
      bool clust;
      for (int k = 0; k < 3; k++) {
        d[k] = x[3*j + k] - x1[k];
        dv[k] = v[3*j + k] - v1[k];
      }
      idx = 0;
      rsq = distance(d, size, &idx);
      /* check if the connection fits energy criterion! */
      if (energy) {
        clust = link(rsq, dv, idx, type[i], type[j], expansion);
      } else {
        clust = (idx != 0);
      }
      if (!clust) continue;
      /* to avoid duplicate entries */
      for (int k = 0; k < 6; k++) {
        bool conn;
        conn = (idx >> k) & 1;
        if (conn) {
          bool add = true;
          for (int l = 0; l < count; l++) {
            int wall;
            wall = wallidx[k];
            if (connect[3*l] == index[i] &&
                connect[3*l + 1] == index[j] &&
                connect[3*l + 2] == wall) {
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
  }
  return count;
}
  

double distance(double *d, double size, int *idx) {
  /* If across boundaries the distance is less than the cut distance,
     return the boundaries through which it connects */
  double min, rsq;
  rsq = 0.0;
  min = 5.4;
  for (int k = 0; k < 3; k++) {
    if (d[k] > size/2) {
      d[k] -= size; /* d is < 0! */
      if (-d[k] < min) {
        *idx += direct[k+3];
      }
    }
    else if (d[k] < -size/2) {
      d[k] += size;
      if (d[k] < min) {
        *idx += direct[k];
      }
    }
    rsq += d[k] * d[k];
  }
  if (rsq > min * min) {
    *idx = 0;
  }
  return rsq;
}

double potential(double r) {
  double Vr = 3088.118;
  double Va = 2666.647;
  double ur = 1.7468;
  double ua = 1.6;
  return Vr * exp(-ur * r) / r - Va * exp(-ua * r) / r;
}

bool link(double rsq, double *dv, int idx, int t1, int t2, double expansion){
  double mu = 938.0/2;
  bool clust = false;
  if (idx != 0) {
    clust = true;
    if (t1 == t2) {
      clust = false;
    }
    else {
      double vsq = 0.0;
      for (int k = 0; k < 3; k++) {
        bool conn;
        conn = (idx >> k) & 1;
        if (conn) {
          dv[k] += expansion;
        }
        conn = (idx >> (k + 3)) & 1;
        if (conn) {
          dv[k] -= expansion;
        }
        vsq += dv[k] * dv[k];
      }
      double eng = potential(sqrt(rsq)) + 1.0/2.0 * mu * vsq;
      clust = (eng < 0.0);
    }
  }
  return clust;
}
