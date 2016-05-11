#include "minkowski.h"

void minkowski(double *x, int natoms, double size, double rad, double rcell, double *out){
  /* We will assume that the box is always cubic and centered, hence
     we only have to communicate one value for the whole box.
  */

  int nc, discr, d2;
  bool *lattice;
  double invrcell;
  nc = ceil(size/rcell);
  lattice = (bool *)malloc(nc * nc * nc * sizeof(bool));
  for (int i = 0; i < nc * nc * nc; i++) lattice[i] = false;
  for (int i = 0; i < 4; i++) out[i] = 0.0;
  rcell = size/nc; /* update so cell fits in box */
  invrcell = 1.0/rcell;

  discr = rad / rcell;
  d2 = discr * discr;
  for (int ii = 0; ii < natoms; ii++) {
    int dx, dy, dz;
    dx = (x[ii*3+0] + size/2) * invrcell;
    dy = (x[ii*3+1] + size/2) * invrcell;
    dz = (x[ii*3+2] + size/2) * invrcell;
    for (int i = -discr; i < discr + 1; i++) {
      for (int j = -discr; j < discr + 1; j++) {
        for (int k = -discr; k < discr + 1; k++) {
          if ((i * i + j * j + k * k) > d2) continue;
          if (get_pixel(lattice, nc, dx + i, dy + j, dz + k)) continue;
          set_pixel(lattice, nc, dx + i, dy + j, dz + k, true);
          add_pixel(lattice, nc, dx + i, dy + j, dz + k, out);
        }
      }
    }
  }
  free(lattice);
  /* scaling of magnitudes */
  out[0] *= rcell * rcell * rcell;
  out[1] *= rcell * rcell;
  out[2] *= rcell;
  return;
}

void add_pixel(bool *lattice, int nc, int i, int j, int k, double *out) {
  int nfaces = 0;
  int nedges = 0;
  int nvert  = 0;

  for (int l = -1; l < 3; l+=2) {
    int li, lj, lk;
    int fxl, fyl, fzl;
    li = i + l;
    lj = j + l;
    lk = k + l;

    fxl = 1 - get_pixel(lattice, nc, li, j, k); // 6 faces adjacent: x, y and z
    fyl = 1 - get_pixel(lattice, nc, i, lj, k);
    fzl = 1 - get_pixel(lattice, nc, i, j, lk);

    nfaces+= fxl + fyl + fzl;

    for (int m = -1; m < 3; m+=2) {
      int mj, mk;
      int fym, fzm;
      int exylm, eyzlm, exzlm;
      mj = j + m;
      mk = k + m;

      fym = 1 - get_pixel(lattice, nc, i, mj, k);
      fzm = 1 - get_pixel(lattice, nc, i, j, mk);

      exylm = 1 - get_pixel(lattice, nc, li, mj, k); // 12 edges adjacent: xy, yz, xz
      eyzlm = 1 - get_pixel(lattice, nc, i, lj, mk);
      exzlm = 1 - get_pixel(lattice, nc, li, j, mk);

      nedges += fxl * fym * exylm;
      nedges += fyl * fzm * eyzlm;
      nedges += fxl * fzm * exzlm;

      for (int n = -1; n < 3; n+=2) {
        int nk;
        int fzn;
        int eyzmn, exzln;
        int vxyzlmn;
        nk = k + n;

        fzn = 1 - get_pixel(lattice, nc, i, j, nk);

        eyzmn = 1 - get_pixel(lattice, nc, i, mj, nk);
        exzln = 1 - get_pixel(lattice, nc, li, j, nk);

        vxyzlmn = 1 - get_pixel(lattice, nc, li, mj, nk); // 8 vertices adjacent: xyz

        nvert += fxl * fym * fzn * exylm * eyzmn * exzln * vxyzlmn;
      }
    }
  }

  out[0]+= v_body;
  out[1]+= s_body + s_face * nfaces;
  out[2]+= c_body + c_face * nfaces + c_edge * nedges;
  out[3]+= e_body + e_face * nfaces + e_edge * nedges + e_vert * nvert;
  return;
}

bool get_pixel(bool *lattice, int nc, int i, int j, int k){
  if (i < 0) i+=nc;
  if (j < 0) j+=nc;
  if (k < 0) k+=nc;

  if (i >= nc) i-=nc;
  if (j >= nc) j-=nc;
  if (k >= nc) k-=nc;

  return lattice[i * nc * nc + j * nc + k];
}

void set_pixel(bool *lattice, int nc, int i, int j, int k, bool value){
  if (i < 0) i+=nc;
  if (j < 0) j+=nc;
  if (k < 0) k+=nc;

  if (i >= nc) i-=nc;
  if (j >= nc) j-=nc;
  if (k >= nc) k-=nc;

  lattice[i * nc * nc + j * nc + k] = value;
  return;
}
