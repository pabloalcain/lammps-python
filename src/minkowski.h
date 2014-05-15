#ifndef MINKOWSKI_H
#define MINKOWSKI_H

#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

using namespace LAMMPS_NS;

extern "C" {
  void minkowski(LAMMPS *lmp, double rmin, double rcell, double *array);
}
void add_pixel(int, int, int);
bool get_pixel(int, int, int);
void set_pixel(int, int, int, bool);
#endif
