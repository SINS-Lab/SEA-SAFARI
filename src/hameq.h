#ifndef HAMEQ_H_INCLUDED
#define HAMEQ_H_INCLUDED
#include "ion.h"

/**
 * This computes the forces to apply on the lattice
 *
 *
 * @param ion - the ion to interact with the lattice
 * @param lattice - the lattice to scatter off
 * @param dt - timestep for interaction
 */
void run_hameq(ion &ion, lattice &lattice, double dt);

/**
 * @param ion - the ion to interact with the lattice
 * @param lattice - the lattice to scatter off
 * @param dt - timestep for interaction
 */
void apply_hameq(ion &ion, lattice &lattice, double dt);

#endif // HAMEQ_H_INCLUDED
