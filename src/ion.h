#ifndef ION_H_INCLUDED
#define ION_H_INCLUDED
#include "lattice.h"

/**
 * This is a particle which can fly by the lattice.
 * It inherits the site class, which gives it the same
 * position/momentum variables. it then also track
 * additional information, such as charge and nearby things
 */
class Ion : public Site
{
public:
    //Charge of the ion
    int q;

    int index = 0;

    int steps = 0;
    double time = 0;

    int max_n = 0;

    int near = 0;
    //Nearest distance
    double r_min = 1e3;
    //This is the last place the near things were updated.
    int last_index = -1;
    //Last distance checked for nearby atoms.
    int last_radius = 1;

    Site near_sites[100];

    //potential energy the particle is in.
    double V = 0;
    //potential energy the particle is in.
    double V_t = 0;

    double max_site_displacement = 0;
    double max_site_momentum = 0;
    /**
     * Sets the initial conditions for the ion.
     *
     *
     * @param eV - energy in eV
     * @param theta0 - Theta angle for incoming beam
     * @param phi0 - Phi angle for the incoming beam
     * @param x - x impact parameter
     * @param y - y impact parameter
     */
    void set_KE(double eV, double theta0, double phi0,double x, double y);

    //Returns the number of nearby atom found in the given lattice.
    //This also updates the near_sites and near_atoms arrays.
    //Call check_distances() before using near_dists.
    int fill_nearest(Lattice &lattice, int radius, int target_number);
};

/**
 * This calculates the trajectory for the ion, 
 * scattering off the given lattice
 *
 * @param ion - the ion to interact with the lattice
 * @param lattice - the lattice to scatter off
 * @param log - if true, will log the individual trajectory
 */
void traj(Ion &ion, Lattice &lattice, bool log);

#endif // ION_H_INCLUDED
