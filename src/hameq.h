#ifndef HAMEQ_H_INCLUDED
#define HAMEQ_H_INCLUDED
#include "ion.h"

void run_hameq(ion &ion, lattice &lattice, double dt, bool at_predicted, double*xp, double*yp, double*zp);

void apply_hameq(ion &ion, lattice &lattice, double dt);

#endif // HAMEQ_H_INCLUDED
