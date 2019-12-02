#pragma once
#include "lattice.h"
#include <iostream>
#include "ion.h"



/**
 * This updates the value of near for this ion. it will search for
 * atoms within radius of the ion, it will only find up to the 
 * first target_number of ions found.
 * 
 * @param ion_ptr - pointer to the involved ion, this is null for inter-lattice checks
 * @param site - site to look for near to.
 * @param lattice - the lattice containing atoms to look for
 * @param radius - max search distance for atoms
 * @param target_number - max value of near to allow
 * @param max_rr - maximum distance squared to include (checked after radius)
 * @param re_sort - whether the list needs to be re-sorted by nearest
 */ 
int fill_nearest(Ion* ion_ptr, Site &site, Lattice &lattice, int radius, int target_number, double max_rr, bool re_sort);

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
                        bool *froze, bool *left, double& E);

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
void traj(Ion &ion, Lattice &lattice, bool& log, bool& xyz);
