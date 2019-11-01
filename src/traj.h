#pragma once
#include "lattice.h"
#include <iostream>
#include "ion.h"

/**
 * Sets the flags indicated by the various pointers, by checking
 * conditions such as z-position, energy, and x/y.
 * 
 * @param ion - the ion to validate
 * @param E - total energy of the ion
 * 
 * @param buried - deeper than BDIST into the surface
 * @param off_edge - ran off the edge of the crystal
 * @param stuck - has less energy than SENRGY
 * @param froze - took too many steps to compute
 * @param left - left the surface, ie got above Z0
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
 * @param log - log the trajectory to .traj file
 * @param xyz - log the xyz trajectory
 */ 
void traj(Ion &ion, Lattice &lattice, bool log, bool xyz);
