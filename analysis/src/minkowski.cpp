/* -------------------------------------------------------------------
   Esta rutina calcula las funcionales de minkowski, partiendo de la
   salida de LAMMPS para x, indistintamente de los tipos (se podría
   eventualmente hacer algo similar para calcular los minkowski de los
   neutrones y los protones por separado).

   Por ahora requiere sí o sí una caja cúbica y de dimensión 3.  Hay
   una pequeña modificación respecto de la interpretación que puede
   llevar a pequeñas incongruencias con las aplicaciones anteriores
   del algoritmo: Ahora pasamos como información el tamaño de la celda
   que vamos a tomar, pero después lo reajustamos ligeramente para
   asegurarnos que entre un número entero de esas celdas en la caja.

   No se me ocurre cómo paralelizarlo, pero igual creo que esta
   implementación es bastante más rápida que la que veníamos
   usando. Básicamente, a medida que digitalizamos, vamos calculando
   las funcionales todas de una.

------------------------------------------------------------------- */

#include <iostream>
#include <math.h>
#include "mpi.h"

#ifndef DIM
#define DIM 3
#endif

#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
#include "minkowski.h"
#include "domain.h"
#include "memory.h"

using namespace LAMMPS_NS;

double *array;
double ***lattice;
double *scratch;
int nx, ny, nz;


void minkowski(LAMMPS *lmp, double rad, double rcell, double *array)
{  
  int natoms = static_cast<int> (lmp->atom->natoms);
  double *x = new double[DIM*natoms];
  int discr;
  
  // This works for non-cubic boxes, but it introduces an error
  // of order cell/boxsz
  nx = ( lmp->domain->subhi[0] - lmp->domain->sublo[0] ) / rcell + 1;
  ny = ( lmp->domain->subhi[1] - lmp->domain->sublo[1] ) / rcell + 1;
  nz = ( lmp->domain->subhi[2] - lmp->domain->sublo[2] ) / rcell + 1;
  
  // Update value of cell size. See previous comment
  rcell = ( lmp->domain->boxhi[0] - lmp->domain->boxlo[0] ) / nx;


  // Get discrete size of particles
  
  discr = rad / rcell;
  
  lmp->memory->create(lattice,nx,ny,nz,"minkowski:lattice");
  
  lammps_gather_atoms(lmp,"x",1,3,x);
  
  lmp->memory->create(scratch,4,"minkowski:scratch");
  
  for (int i = 0; i < 4; i++)
    scratch[i] = 0;
  
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
	lattice[i][j][k] = false;
  


   int discr2;
   int dx, dy, dz; // discrete coordinates
  
  discr2 = discr * discr;
  
  /* Each proc takes care of his own subspace, eventually adding up
     the lattice points next to his. This way, communication is
     reduced to a minimum. */
  
  for (int ii = 0; ii < natoms; ii++) {
    dx = ( x[ii*3+0] - lmp->domain->sublo[0] ) / rcell;
    dy = ( x[ii*3+1] - lmp->domain->sublo[1] ) / rcell;
    dz = ( x[ii*3+2] - lmp->domain->sublo[2] ) / rcell;
    for (int i = -discr; i < discr + 1; i++) {
      for (int j = -discr; j < discr + 1; j++) {
	for (int k = -discr; k < discr + 1; k++) {
	  if ((i * i + j *j +  k *k ) > discr2) continue;
	  if (get_pixel(dx + i, dy + j,dz + k)) continue;
	  set_pixel(dx + i, dy + j,dz + k, true);
	  add_pixel(dx + i, dy + j,dz + k);
	}
      }
    }
  }

  /* Rescale magnitudes */
  
  scratch[0] *= rcell * rcell * rcell;
  scratch[1] *= rcell * rcell;
  scratch[2] *= rcell;
  

  /* Sum values across procs */
  MPI_Allreduce(scratch,array,4,MPI_DOUBLE,MPI_SUM,lmp->world);

  lmp->memory->destroy(scratch);
  return;
}

void add_pixel(int i, int j, int k)
{
  int li, lj, lk, mj, mk, nk;
  
  int fxl, fyl, fzl, fym, fzm, fzn;

  int exylm, eyzlm, exzlm, eyzmn, exzln;
  
  int vxyzlmn;


  const int v_body =  1; // view M-DR algorithm for explanation
  const int s_body = -6;
  const int s_face =  2;
  const int c_body =  3;
  const int c_face = -2;
  const int c_edge =  1;
  const int e_body = -1;
  const int e_face =  1;
  const int e_edge = -1;
  const int e_vert =  1;

 
  int nfaces = 0;
  int nedges = 0;
  int nvert  = 0;

  for (int l = -1; l < 3; l+=2) {
    li = i + l;
    lj = j + l;
    lk = k + l;

    fxl = 1 - get_pixel(li, j, k); // 6 faces adjacent: x, y and z
    fyl = 1 - get_pixel(i, lj, k); 
    fzl = 1 - get_pixel(i, j, lk);
    
    nfaces+= fxl + fyl + fzl;

    for (int m = -1; m < 3; m+=2) {
      mj = j + m;
      mk = k + m;

      fym = 1 - get_pixel(i, mj, k); 
      fzm = 1 - get_pixel(i, j, mk);

      exylm = 1 - get_pixel(li, mj, k); // 12 edges adjacent: xy, yz, xz
      eyzlm = 1 - get_pixel(i, lj, mk);
      exzlm = 1 - get_pixel(li, j, mk); 
      
      nedges += fxl * fym * exylm;
      nedges += fyl * fzm * eyzlm;
      nedges += fxl * fzm * exzlm;
    
      for (int n = -1; n < 3; n+=2) {
	nk = k + n;
	
	fzn = 1 - get_pixel(i, j, nk);
	
	eyzmn = 1 - get_pixel(i, mj, nk);
	exzln = 1 - get_pixel(li, j, nk);
	
	vxyzlmn = 1 - get_pixel(li, mj, nk); // 8 vertices adjacent: xyz
	
	nvert += fxl * fym * fzn * exylm * eyzmn * exzln * vxyzlmn;
      }
    }
  }
  scratch[0]+= v_body;
  scratch[1]+= s_body + s_face * nfaces;
  scratch[2]+= c_body + c_face * nfaces + c_edge * nedges;
  scratch[3]+= e_body + e_face * nfaces + e_edge * nedges + e_vert * nvert;
  return;
}

/* ---------------------------------------------------------------------- */

bool get_pixel(int i, int j, int k){
  if (i < 0) i+=nx;
  if (j < 0) j+=ny;
  if (k < 0) k+=nz;

  if (i >= nx) i-=nx;
  if (j >= ny) j-=ny;
  if (k >= nz) k-=nz;

  return lattice[i][j][k];
}

/* ---------------------------------------------------------------------- */

void set_pixel(int i, int j, int k, bool val){
  if (i < 0) i+=nx;
  if (j < 0) j+=ny;
  if (k < 0) k+=nz;

  if (i >= nx) i-=nx;
  if (j >= ny) j-=ny;
  if (k >= nz) k-=nz;

  lattice[i][j][k] = val;
}
