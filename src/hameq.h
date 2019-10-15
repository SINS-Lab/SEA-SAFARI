#ifndef HAMEQ_H_INCLUDED
#define HAMEQ_H_INCLUDED
#include "ion.h"

void run_hameq(ion &ion, lattice &lattice, double dt, double *force);

void apply_hameq(ion &ion, lattice &lattice, double dt);

#endif // HAMEQ_H_INCLUDED
