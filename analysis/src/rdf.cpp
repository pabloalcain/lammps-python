/* -------------------------------------------------------------------
   This routine takes a LAMMPS object and calculates the radial
   distribution function, with an O(N^2) algorithm in order to be able
   to get as long r as boxsize/2, with all type pairs.


   TODO: 

   [ ] Add some wrappers, so gofr can be called without a pointer to
   the whole LAMMPS object, but also with the actual information: x,
   type, boxsize, etc.

   [ ] Consider adding multi-replica information.

   [ ] Parallelization is too naive, and can be changed if we find
   an intelligent way to map the pair (i, j) to a unique id k.

   
   ------------------------------------------------------------------- */


#ifndef DIM
#define DIM 3
#endif

using namespace LAMMPS_NS;

void rdf(LAMMPS *lmp, int nbins, double rmax, double **rdf)
{  
  
  int mystart, myend; // MPI arguments for memory distribution
  int me,nprocs;


  double boxlo, boxhi, boxsz;
  int npart = static_cast<int> (lmp->atom->natoms);
  int nprot=0;
  int nneut=0;
  int npair=4; // different type of pairs we can make: if there are T
	       // type of particles, npair = 1 + T*(T+1)/2 
               // [1 is for "all v all"]
  double *x = new double[DIM*npart];
  int *type = new int[npart];

  double rmax;
  double delr;
  double hist[nbins][npair];
  double histall[nbins][npair];
  double array[nbins][npair+1];

  boxlo = *((double *)lmp->lammps_extract_global(lmp,"boxxlo"));
  boxhi = *((double *)lmp->lammps_extract_global(lmp,"boxxhi"));
  boxsz=boxhi-boxlo;
  rmax=boxsz/2;
  
  delr=rmax/nbins;
  lmp->lammps_gather_atoms(lmp,"x",1,3,x);
  lmp->lammps_gather_atoms(lmp,"type",0,1,type);
  
  for (int i = 0; i < nbins; i++)
    for (int j = 0; j < npair; j++)
      hist[i][j] = 0.0;
  
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

  int i, j;
  for (i=mystart; i<myend; i++){
    if (type[i] == 1)
      nprot++;
    
    for (j=i+1; j<npart; j++){
      double dist=0.0;
      int bin, type_pair;
      for (int k= 0; k<DIM; k++){
	double dr=x[i*DIM+k]-x[j*DIM+k];
	if (dr > boxsz/2)  dr-=boxsz;
	if (dr < -boxsz/2) dr+=boxsz;
	dist+=dr*dr;
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
  
  int scratch;
  MPI_Allreduce(&nprot, &scratch, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  nprot=scratch;
  nneut=npart-nprot;

  
  if (me==0){
    for (int i = 0; i < nbins; i++){
      rdf[i][0]=(i+0.5)*delr;
      for (int j = 0; j < npair; j++)
	rdf[i][j+1] = (boxsz * boxsz * boxsz)*histall[i][j]/(array[i][0]*array[i][0])/(4*3.1416*delr);
    }
    for (int i = 0; i < nbins; i++){
      rdf[i][1]=2 * rdf[i][1] / ( npart * npart );
      rdf[i][2]=2 * rdf[i][2] / ( nprot * nprot );
      rdf[i][3]=2 * rdf[i][3] / ( nneut * nneut );
      rdf[i][4]=1 * rdf[i][4] / ( nprot * nneut );
    }
    
  }
  return
}
