#pragma once
#include "particles.h" // for Site

/**
 * This is a particle which can fly by the lattice.
 * It inherits the site class, which gives it the same
 * position/momentum variables. it then also track
 * additional information, such as charge and nearby things
 */
class Ion : public Site
{
public:
    // The initial energy of this ion
    double E0 = 0;
    // potential energy the particle is in.
    double V = 0;
    // Used to track maximum number of active sites when in dynamic
    // neighbour tracking mode.
    int max_active = 0;
    // Things with might sputtur off the surface
    Site** sputter = NULL;
    int sputtered = 0;
    
    Ion()
    {
        // In the parent class, this is left null, unless
        // lattice springs are used.
        // So for Ion, we initialize it at this size.
        near_sites = new Site *[MAX_NEAR];
        q = 1;
    }

    ~Ion()
    {
        if (near_sites != NULL)
            delete near_sites;
        near_sites = NULL;
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
    void set_KE(double E0, double theta0, double phi0, double x, double y);
};