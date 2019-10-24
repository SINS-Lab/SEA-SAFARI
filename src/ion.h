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
    int q = 1;
    //Index of the ion, each one fired should be different.
    int index = 0;
    //The initial energy of this ion
    double E0 = 0;

    //Number of traj steps done, note that this is not
    //the same as the number of time steps, as this is
    //incremented even if it reduces timestep and tries over.
    int steps = 0;
    //Total flight time of the ion
    double time = 0;

    //Maximum number of nearby atoms for this ion
    int max_n = 0;
    //Current number of nearby atoms.
    int near = 0;
    //Nearest distance ever
    double r_min = 1e3;
    //Nearest distance when finding nearby, squared
    double rr_min_find = 1e3;
    //This is the last place the near things were updated.
    //Setting this back to -1 will force a re-check of near.
    int last_index = -1;
    //Last distance checked for nearby atoms.
    int last_radius = 1;

    //All of the sites nearby, only guarenteed to be filled
    //with site up to near
    Site *near_sites[100];

    //potential energy the particle is in.
    double V = 0;
    //potential energy the particle is in.
    double V_t = 0;

    /**
     * Sets the initial conditions for the ion.
     *
     *
     * @param theta0 - Theta angle for incoming beam
     * @param phi0 - Phi angle for the incoming beam
     * @param x - x impact parameter
     * @param y - y impact parameter
     */
    void set_KE(double theta0, double phi0,double x, double y);

    /**
     * This updates the value of near for this ion. it will search for
     * atoms within radius of the ion, it will only find up to the 
     * first target_number of ions found.
     * 
     * @param lattice - the lattice containing atoms to look for
     * @param radius - max search distance for atoms
     * @param target_number - max value of near to allow
     */ 
    int fill_nearest(Lattice &lattice, int radius, int target_number);
};

#endif // ION_H_INCLUDED
