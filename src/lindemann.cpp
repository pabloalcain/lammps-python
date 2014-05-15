/* -------------------------------------------------------------------
   Esta rutina calcula g(r) para TODOS los pares, partiendo de la
   salida de LAMMPS para x y el type que puede ser 1 o 2. Está
   paralelizado, pero por ahora muy burdamente. Falta encontrar una
   forma de mapear un par único (i,j) a un índice k.

   Por ahora requiere sí o sí una caja cúbica y de dimensión 3 (aunque
   creo que eso se puede cambiar trivialmente cambiando DIM)
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

using namespace LAMMPS_NS;

double lindemann(LAMMPS *lmp)
{  
  
  int mystart, myend; // MPI arguments for memory distribution
  int me,nprocs;

  
  double boxlo, boxhi, boxsz;
  int npart = static_cast<int> (lmp->atom->natoms);
  int nprot=0;
  int nprotall=0;
  int nneutall=0;

  double r;
  double r2;
  
  double *x = new double[DIM*npart];
  int *type = new int[npart];
  int *id = new int[natoms];
  
  boxlo = *((double *)lammps_extract_global(lmp,"boxxlo"));
  boxhi = *((double *)lammps_extract_global(lmp,"boxxhi"));
  boxsz=boxhi-boxlo;
  
  lammps_gather_atoms(lmp,"x",1,3,x);
  lammps_gather_atoms(lmp,"type",0,1,type);
  lammps_gather_atoms(lmp,"id",0,1,id);
  
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  
  mystart = (npart / nprocs) * me;
  if (npart % nprocs > me){
    mystart += me;
    myend = mystart + (npart / nprocs) + 1;
  }else{
    mystart += npart % nprocs;
    myend = mystart + (npart / nprocs);
  }

  int i;
  for (i=mystart; i<myend; i++){
    if (type[i] == 1)
      nprot++;
    
    for (int k= 0; k<DIM; k++){
      double r=x[i*DIM+k];
      double r2+=r*r;
    }
    dist=sqrt(dist);
    if ( dist > rmax ) continue;
    // choose the pair type:
    // all             -> 0
    // proton-proton   -> 1
    // neutron-neutron -> 2
    // proton-neutron  -> 3
    if (type[i] == 1 && type[j] == 1)
      type_pair=1;
    else if (type[i] == 2 && type[j] == 2)
	type_pair=2;
      else
	type_pair=3;    
      bin=(int)((dist*nbins)/rmax);
      hist[bin][0]++;
      hist[bin][type_pair]++;
    }
  }
  
  MPI_Allreduce(hist[0], histall[0], npair*nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nprot, &nprotall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  nneutall=npart-nprotall;
  
  if (me==0){
    for (int i = 0; i < nbins; i++){
      array[i][0]=(i+0.5)*delr;
      for (int j = 0; j < npair; j++)
	array[i][j+1] = (boxsz * boxsz * boxsz)*histall[i][j]/(array[i][0]*array[i][0])/(4*3.1416*delr);
    }
    for (int i = 0; i < nbins; i++){
      array[i][1]=2 * array[i][1] / ( npart * npart );
      array[i][2]=2 * array[i][2] / ( nprotall * nprotall );
      array[i][3]=2 * array[i][3] / ( nneutall * nneutall );
      array[i][4]=1 * array[i][4] / ( nprotall * nneutall );
    }
    
    // for (int i = 0; i < nbins; i++) {
    //   for (int j = 0; j < npair+1; j++) 
    // 	std::cout<<array[i][j]<<" ";
    //   std::cout<<std::endl;
    // }
  }

  return array;
}
