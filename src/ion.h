#pragma once
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
    //The initial energy of this ion
    double E0 = 0;

    //Number of traj steps done, note that this is not
    //the same as the number of time steps, as this is
    //incremented even if it reduces timestep and tries over.
    int steps = 0;
    //Total flight time of the ion
    double time = 0;

    //This is a weighting factor for the ion's detectability.
    double weight = 1;

    //potential energy the particle is in.
    double V = 0;

    Ion()
    {
        // In the parent class, this is left null, unless
        // lattice springs are used.
        // So for Ion, we initialize it at this size.
        near_sites = new Site*[256];
    }

    ~Ion()
    {
        delete near_sites;
    }

    /**
     * Overrides the reset for Site, this is more specific
     * for the ion
     */ 
    void reset();

    /**
     * Sets the initial conditions for the ion.
     *
     * @param E0 - The initial energy of the ion.
     * @param theta0 - Theta angle for incoming beam
     * @param phi0 - Phi angle for the incoming beam
     * @param x - x impact parameter
     * @param y - y impact parameter
     */
    void set_KE(double E0, double theta0, double phi0,double x, double y);
};