#ifndef TRAJ_H_INCLUDED
#define TRAJ_H_INCLUDED
#include "lattice.h"
#include "ion.h"

/**
 * Sets the flags indicated by the various pointers, by checking
 * conditions such as z-position, energy, and x/y.
 * 
 * @param ion - the ion to validate.
 * @param E - kinetic energy of the ion.
 * 
 */ 
bool validate(Ion &ion, bool *buried, bool *off_edge, bool *stuck,
                        bool *froze, bool *left, double E);

/**
 * Computes the trajectory of the ion, if log, it prints extra
 * information about the trajectory and status if the ion.
 * 
 * if xyz, it outputs the coordinates of the ion and lattice for
 * each fram of computation
 * 
 * @param ion - the ion to scatter off the lattice
 * @param lattice - the lattice to scatter off.
 */ 
void traj(Ion &ion, Lattice &lattice, bool log, bool xyz);

#endif // TRAJ_H_INCLUDED