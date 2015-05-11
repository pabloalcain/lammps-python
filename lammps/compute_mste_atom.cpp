/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_mste_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSTEAtom::ComputeMSTEAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute cluster/atom command");

  double cutoff = force->numeric(FLERR,arg[3]);
  cutsq = cutoff*cutoff;

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_forward = 1;

  nmax = 0;
  clusterID = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMSTEAtom::~ComputeMSTEAtom()
{
  memory->destroy(clusterID);
}

/* ---------------------------------------------------------------------- */

void ComputeMSTEAtom::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use compute cluster/atom unless atoms have IDs");
  if (force->pair == NULL)
    error->all(FLERR,"Compute cluster/atom requires a pair style be defined");
  if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,
               "Compute cluster/atom cutoff is longer than pairwise cutoff");

  // need an occasional full neighbor list
  // full required so that pair of atoms on 2 procs both set their clusterID

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"cluster/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute cluster/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMSTEAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeMSTEAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double delpx,delpy,delpz,psq;
  double eng,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow clusterID array if necessary

  if (atom->nlocal+atom->nghost > nmax) {
    memory->destroy(clusterID);
    nmax = atom->nmax;
    memory->create(clusterID,nmax,"cluster/mste:clusterID");
    vector_atom = clusterID;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // if group is dynamic, insure ghost atom masks are current

  if (group->dynamic[igroup]) {
    commflag = 0;
    comm->forward_comm_compute(this);
  }

  // every atom starts in its own cluster, with clusterID = atomID

  tagint *tag = atom->tag;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) clusterID[i] = tag[i];
    else clusterID[i] = 0;
  }

  // loop until no more changes on any proc:
  // acquire clusterIDs of ghost atoms
  // loop over my atoms, checking distance to neighbors
  // if both atoms are in cluster, assign lowest clusterID to both
  // iterate until no changes in my atoms
  // then check if any proc made changes

  commflag = 1;
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  int change,done,anychange;

  while (1) {
    comm->forward_comm_compute(this);

    change = 0;
    while (1) {
      done = 1;
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (!(mask[i] & groupbit)) continue;

        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
	itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
	  factor_lj = special_lj[sbmask(j)];
	  factor_coul = special_coul[sbmask(j)];
          j &= NEIGHMASK;
          if (!(mask[j] & groupbit)) continue;

          if (clusterID[i] == clusterID[j]) continue;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          jtype = type[j];
          rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < cutsq) {
            delpx = mass[itype] * v[tag[i]][0] - mass[jtype] * v[tag[j]][0];
            delpy = mass[itype] * v[tag[i]][1] - mass[jtype] * v[tag[j]][1];
            delpz = mass[itype] * v[tag[i]][2] - mass[jtype] * v[tag[j]][2];
            psq = delpx*delpx + delpy*delpy + delpz*delpz;
            eng = force->pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
            eng += psq / 2.0 * (mass[itype] + mass[jtype])/(mass[itype] * mass[jtype]);
	    if (eng < 0.0) {
              clusterID[i] = clusterID[j] = MIN(clusterID[i],clusterID[j]);
              done = 0;
	    }
          }
        }
      }
      if (!done) change = 1;
      if (done) break;
    }

    // stop if all procs are done

    MPI_Allreduce(&change,&anychange,1,MPI_INT,MPI_MAX,world);
    if (!anychange) break;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeMSTEAtom::pack_forward_comm(int n, int *list, double *buf,
				       int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  if (commflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = clusterID[j];
    }
  } else {
    int *mask = atom->mask;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(mask[j]).d;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeMSTEAtom::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  if (commflag)
    for (i = first; i < last; i++) clusterID[i] = buf[m++];
  else {
    int *mask = atom->mask;
    for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMSTEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
