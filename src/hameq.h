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
 * @param predicted - if true, this does the calculations for
 *                    after time dt has passed, otherwise is
 *                    done for the current location
 */
void run_hameq(Ion &ion, Lattice &lattice, double dt, bool predicted);

/**
 * @param ion - the ion to interact with the lattice
 * @param lattice - the lattice to scatter off
 * @param dt - timestep for interaction
 */
void apply_hameq(Ion &ion, Lattice &lattice, double dt);

#endif // HAMEQ_H_INCLUDED
