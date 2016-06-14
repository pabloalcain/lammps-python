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

double potential(double r, int t1, int t2) {
  double Vr = 3088.118;
  double Va = 2666.647;
  double V0 = 373.118;
  double Vc = 1.44;

  double ur = 1.7468;
  double ua = 1.6;
  double u0 = 1.5;
  double uc = 0.05;

  if (t1 == t2 && t1 == 1) {
    return V0 * exp(-u0 * r) / r;
  } else if (t1 == t2 && t1 == 2) {
    return V0 * exp(-u0 * r) / r + Vc * exp(-uc * r) / r;
  } else {
    return Vr * exp(-ur * r) / r - Va * exp(-ua * r) / r;
  }
}
