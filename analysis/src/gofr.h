#ifndef GOFR_H
#define GOFR_H

#include <iostream>
#include <math.h>
#include "mpi.h"
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
using namespace LAMMPS_NS;
extern "C" {
  void rdf(LAMMPS *, int, double, double*);
}
#endif
