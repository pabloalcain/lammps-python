#ifndef RDF_H
#define RDF_H

#include <iostream>
#include <math.h>
#include "mpi.h"

#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

using namespace LAMMPS_NS;
void rdf(LAMMPS, int, double, double);
#endif
