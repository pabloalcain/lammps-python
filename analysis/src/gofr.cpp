/* -------------------------------------------------------------------
   Esta rutina calcula g(r) para TODOS los pares, partiendo de la
   salida de LAMMPS para x y el type que puede ser 1 o 2. Está
   paralelizado, pero por ahora muy burdamente. Falta encontrar una
   forma de mapear un par único (i,j) a un índice k.

   Por ahora requiere sí o sí una caja cúbica y de dimensión 3 (aunque
   creo que eso se puede cambiar trivialmente cambiando DIM)

   Otro "inconveniente" a tener en cuenta es la no-portabilidad: asume
   que hay dos tipos de partículas y las llama protones y neutrones.


   TODO: 
   
   [1 ]: Buscar alguna forma (parecida a la de LAMMPS?) de keep tabs
   on the número de pares que se consideran, y automatizarlo un poco.
   
   [5 ]: Permitir cajas de otros tamaños/otras formas?

------------------------------------------------------------------- */

#include "gofr.h"

using namespace LAMMPS_NS;
using namespace LAMMPS_NS::MathConst;

void rdf(LAMMPS *lmp, int nbins, double rmax, double *array)
{  
  int npairs;            // # of rdf pairs
  double delr,delrinv;   // bin width and its inverse
  int ***rdfpair;        // map 2 type pair to rdf pair for each histo
  int **nrdfpair;        // # of histograms for each type pair
  int *ilo,*ihi,*jlo,*jhi;
  double **hist;         // histogram bins
  double **histall;      // summed histogram bins across all procs
  
  int *typecount;
  int *icount,*jcount;



  int natoms = static_cast<int> (lmp->atom->natoms);
  int ntypes = static_cast<int> (lmp->atom->ntypes);

  double *x = new double[3*natoms];
  lammps_gather_atoms(lmp,"x",1,3,x);
  int *type = new int[natoms];
  lammps_gather_atoms(lmp,"type",0,1,type);
  
  if (ntypes == 1) npairs = 1;
  else npairs = 1 + (ntypes * (ntypes + 1) ) / 2;
  
  lmp->memory->create(rdfpair,npairs,ntypes+1,ntypes+1,"rdfall:rdfpair");
  lmp->memory->create(nrdfpair,ntypes+1,ntypes+1,"rdfall:nrdfpair");

  ilo = new int[npairs];
  ihi = new int[npairs];
  jlo = new int[npairs];
  jhi = new int[npairs];

  ilo[0] = 1; ihi[0] = ntypes;
  jlo[0] = 1; jhi[0] = ntypes;
  
  if (npairs > 1) {
    int idx = 1;
    for (int i = 1; i <= ntypes; i++) {
      for (int j = i; j <= ntypes; j++) {
	ilo[idx] = i; ihi[idx] = i;
	jlo[idx] = j; jhi[idx] = j;
	idx++;
      }
    }
  }

  int i,j;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      nrdfpair[i][j] = 0;

  for (int m = 0; m < npairs; m++)
    for (i = ilo[m]; i <= ihi[m]; i++)
      for (j = jlo[m]; j <= jhi[m]; j++)
        rdfpair[nrdfpair[i][j]++][i][j] = m;

  lmp->memory->create(hist,npairs,nbins,"rdfall:hist");
  lmp->memory->create(histall,npairs,nbins,"rdfall:histall");
  typecount = new int[ntypes+1];
  icount = new int[npairs];
  jcount = new int[npairs];

  delr = rmax / nbins;
  
  delrinv = 1.0/delr;

  // set 1st column of output array to bin coords

  for (i = 0; i < nbins; i++)
    array[i*(2*npairs+1)+0] = (i+0.5) * delr;

  // count atoms of each type that are also in group

  
  for (i = 1; i <= ntypes; i++) typecount[i] = 0;
  for (i = 0; i < natoms; i++)
    typecount[type[i]]++;

  // icount = # of I atoms participating in I,J pairs for each histogram
  // jcount = # of J atoms participating in I,J pairs for each histogram

  for (int m = 0; m < npairs; m++) {
    icount[m] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++) icount[m] += typecount[i];
    jcount[m] = 0;
    for (i = jlo[m]; i <= jhi[m]; i++) jcount[m] += typecount[i];
  }

  int *scratch = new int[npairs];
  MPI_Allreduce(icount,scratch,npairs,MPI_INT,MPI_SUM,lmp->world);
  for (i = 0; i < npairs; i++) icount[i] = scratch[i];
  MPI_Allreduce(jcount,scratch,npairs,MPI_INT,MPI_SUM,lmp->world);
  for (i = 0; i < npairs; i++) jcount[i] = scratch[i];
  delete [] scratch;

  
  int itype,jtype,ipair,jpair,ibin,ihisto;
  double xtmp,ytmp,ztmp,delx,dely,delz,r;
  
  double rmax2 = rmax * rmax;
  
  for (i = 0; i < npairs; i++)
    for (j = 0; j < nbins; j++)
      hist[i][j] = 0.0;
  
  int mystart, myend; // MPI arguments for memory distribution
  int me,nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  
  mystart = (natoms / nprocs) * me;
  if (natoms % nprocs > me){
    mystart += me;
    myend = mystart + (natoms / nprocs) + 1;
  }else{
    mystart += natoms % nprocs;
    myend = mystart + (natoms / nprocs);
  }

  for (i = mystart; i < myend; i++) {
    xtmp = x[i*3+0];
    ytmp = x[i*3+1];
    ztmp = x[i*3+2];
    itype = type[i];
    for (j = i+1; j < natoms; j++) {
      jtype = type[j];
      ipair = nrdfpair[itype][jtype];
      jpair = nrdfpair[jtype][itype];
      if (!ipair && !jpair) continue;

      delx = xtmp - x[j*3+0];
      dely = ytmp - x[j*3+1];
      delz = ztmp - x[j*3+2];
      lmp->domain->minimum_image(delx,dely,delz);
      r = delx*delx;
      if (r > rmax2) continue;
      r+= dely*dely;
      if (r > rmax2) continue;
      r+= delz*delz;
      if (r > rmax2) continue;
      r = sqrt(r);
      ibin = static_cast<int> (r*delrinv);

      if (ipair)
        for (ihisto = 0; ihisto < ipair; ihisto++)
          hist[rdfpair[ihisto][itype][jtype]][ibin] += 1.0;
        if (jpair)
          for (ihisto = 0; ihisto < jpair; ihisto++)
            hist[rdfpair[ihisto][jtype][itype]][ibin] += 1.0;
	
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],npairs*nbins,MPI_DOUBLE,MPI_SUM,lmp->world);

  // convert counts to g(r) and coord(r) and copy into output array
  // nideal = # of J atoms surrounding single I atom in a single bin
  //   assuming J atoms are at uniform density

  double constant,nideal,gr,ncoord,rlower,rupper;

  if (lmp->domain->dimension == 3) {
    constant = 4.0*MY_PI / (3.0*lmp->domain->xprd*lmp->domain->yprd*lmp->domain->zprd);

    for (int m = 0; m < npairs; m++) {
      ncoord = 0.0;
      for (ibin = 0; ibin < nbins; ibin++) {
        rlower = ibin*delr;
        rupper = (ibin+1)*delr;
        nideal = constant *
          (rupper*rupper*rupper - rlower*rlower*rlower) * jcount[m];
        if (icount[m]*nideal != 0.0)
          gr = histall[m][ibin] / (icount[m]*nideal);
        else gr = 0.0;
        ncoord += gr*nideal;
        array[ibin*(2*npairs+1)+1+2*m] = gr;
        array[ibin*(2*npairs+1)+2+2*m] = ncoord;
      }
    }

  } else {
    constant = MY_PI / (lmp->domain->xprd*lmp->domain->yprd);

    for (int m = 0; m < npairs; m++) {
      ncoord = 0.0;
      for (ibin = 0; ibin < nbins; ibin++) {
        rlower = ibin*delr;
        rupper = (ibin+1)*delr;
        nideal = constant * (rupper*rupper - rlower*rlower) * jcount[m];
        if (icount[m]*nideal != 0.0)
          gr = histall[m][ibin] / (icount[m]*nideal);
        else gr = 0.0;
        ncoord += gr*nideal;
        array[ibin*(2*npairs+1)+1+2*m] = gr;
        array[ibin*(2*npairs+1)+2+2*m] = ncoord;
      }
    }
  }
  return;
}
