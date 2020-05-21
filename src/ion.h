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

    bool done = false;

    // exit conditions

    // Below -BDIST
    bool buried = false;
    // Took longer than MAX_STEPS
    bool froze = false;
    // Energy lower than SENRGY
    bool stuck = false;
    // Left RAX/RAY
    bool off_edge = false;
    // Too much change in energy over last 2 ticks
    bool discont = false;
    // Above Z1
    bool left = false;
    
    bool resort = true;
    bool reindex = true;
    
    Ion()
    {
        // In the parent class, this is left null, unless
        // lattice springs are used.
        // So for Ion, we initialize it at this size.
        near_sites = new Site *[MAX_NEAR];
        // near_dists = new double *[MAX_NEAR * 6];
        // near_forces = new double *[MAX_NEAR * 6];
        q = 1;
    }

    ~Ion()
    {
        if(sputter!=NULL)
            delete sputter;
        sputter = NULL;
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