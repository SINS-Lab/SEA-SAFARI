#pragma once
#include "ion.h"     // for Ion
#include "lattice.h" // for Lattice

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
void run_hameq(Ion &ion, Lattice *lattice, double dt, bool predicted, double *dr_max);
void run_hameq(std::vector<Ion *> &ions, Lattice *lattice, double dt, bool predicted, double *dr_max);

/**
 * @param ion - the ion to interact with the lattice
 * @param lattice - the lattice to scatter off
 * @param dt - timestep for interaction
 */
void apply_hameq(Ion &ion, Lattice *lattice, double dt);
void apply_hameq(std::vector<Ion*> &ions, Lattice *lattice, double dt);

/**
 * Checks if this site sputtered for this ion, if so, 
 * it will add it to the list of sputter for the ion.
 */ 
bool check_sputter(Ion &ion, Site *s);

/**
 * This updates the current location/momentum of
 * the particle, to the corrected values based on the
 * forces at the current and predicted locations.
 * 
 * Positions are corrected by the error in the
 * values from the two forces, but momenta are
 * updated with just the average of the two forces.
 * 
 * @param s - the particle to update
 * @param dt - time step for this update.
 */
void update_site(Site &s, double dt);

/**
 * This predicts the location of the site, assuming the
 * current value of the momentum and force.
 * 
 * @param s - the particle to predict.
 * @param dt - the time step for the prediction
 */
void predict_site_location(Site &s, double dt);

/**
 * This updates the current location/momentum of
 * the particle, to the corrected values based on the
 * forces at the current and predicted locations.
 * 
 * Positions are corrected by the error in the
 * values from the two forces, but momenta are
 * updated with just the average of the two forces.
 * 
 * This also then updates all of the subsites for this site,
 * note that this function does not get optimized as well as
 * the one that only does a single site.
 * 
 * @param s - the particle to update
 * @param 
 * @param dt - time step for this update.
 */
void update_sites(Site &s, int last_update, double dt);


double compute_error(Site &site, double dt);

void apply_ion_lattice(Ion &ion, Site *s, double *F_at, double *r_i, 
                       double ax, double ay, double az, 
                       double dt, bool predicted, double *F_ion);

void apply_lattice_lattice(Site *s, Site *s2, Ion &ion, double *F_at, 
                           double atomk, double dt, double ax, double ay, double az,
                           bool predicted, bool recoil, bool useLennardJones, bool doubleCount);

void apply_dynamic_lattice(Site *s, Lattice *lattice, int start, Ion &ion, 
                           double *F_at, double atomk, double dt, 
                           double ax, double ay, double az,
                           bool predicted, bool recoil, bool useLennardJones);