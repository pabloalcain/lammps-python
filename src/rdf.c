#include "rdf.h"

void rdf(double *x, int natoms, int nbins, double size, double *gr){
  /* Calculates rdf with PBC in 3D for a cubic box and returns on
     gr. x is an array that has the info in XYZ XYZ XYZ...fashion. The
     rdf is calculated for the whole box up to sqrt(3)/2.*/

  double rmax;
  double p;
  double dr;

  double x1[3];
  double dx, r;
  int idx;
  
  for (int k = 0; k < 2 * nbins; k++) gr[k] = 0.0;
  rmax = size * sqrt(3.0)/2;
  dr = rmax/nbins;
  for (int i = 0; i < natoms; i++) {
    for (int k = 0; k < 3; k++) x1[k] = x[3*i + k];
    for (int j = i+1; j < natoms; j++) {
      r = 0.0;
      for (int k = 0; k < 3; k++) {
        dx = x[3*j + k] - x1[k];
        while (dx > size/2)  dx -= size; 
        while (dx < -size/2) dx += size; 
        r += dx * dx;
      }
      r = sqrt(r);
      idx = (int)(r*nbins/rmax);
      gr[2 * idx + 1] +=1.0;
    }
  }
  for (int k = 0; k < nbins; k++) {
    r = (k+0.5)*rmax/nbins;
    gr[2 * k] = r;
    p = volume((r+dr)/size) - volume(r/size);
    gr[2 * k + 1] /= p * (natoms * natoms) / 2;
  }
  return;
}

double f1(double l) {
  return atan(sqrt(4 * l * l - 2));
}

double f2(double l) {
  double ang = 2 * l * (4 * l * l -3) / (sqrt(4 * l * l - 2) * (4 * l * l + 1));
  return 8 * l * atan(ang);
}
double volume(double l) {
  /* Given the radius gives the volume that is inside the cube of
     length 1 AND the sphere of radius l. See for instance
     https://www.cmu.edu/biolphys/deserno/pdf/gr_periodic.pdf
  */
  
  double v;
  if (l < 0.5) {
    v = 4.0/3.0 * M_PI * l * l * l;
  } else if (l < 1.0/sqrt(2.0)) {
    v = - M_PI/12 * (3 - 36 * l * l +32 * l * l * l);
  }
  else if (l < sqrt(3.0)/2.0) {
    v = - M_PI/4 + 3* M_PI * l * l;
    v += sqrt(4* l * l - 2);
    v += (1 - 12 * l * l) * f1(l);
    v += 2.0/3.0 * l * l * f2(l);
  }
  else v = 1;
  return v;
}
                        
