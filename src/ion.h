#ifndef ION_H_INCLUDED
#define ION_H_INCLUDED
#include "lattice.h"

class ion
{
public:
    //Initial "position", this is the impact parameter
    double r_0[3];
    //Position of the ion
    double r[3];
    //Momentum of the ion
    double p[3];

    //Atom for this ion
    atom atom;
    //Charge of the ion
    int q;

    //Forces on the ion here
    double dp_dt[3];
    //Forces on the ion next step
    double dp_dt_t[3];

    int index = 0;

    int steps = 0;
    double time = 0;

    int max_n = 0;

    int near = 0;
    //Nearest distance
    double r_min = 1e3;
    int last_index = -1;

    site near_sites[100];
    double near_dists[100];
    int near_atoms[100];

    double v_total = 0;

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

    double &operator[](int index)
    {
        return r[index];
    }

    ion();  // Constructor
    ~ion(); // Destructor

    //Returns the number of nearby atom found in the given lattice.
    //This also updates the near_sites and near_atoms arrays.
    //Call check_distances() before using near_dists.
    int fill_nearest(lattice &lattice);
    //Updates the values in near_dists for near sites.
    //The coordinates of the ion are the given arguments,
    //This is for running it at possible predicted locations instead.
    void check_distances(double x, double y, double z);
};

/**
 * This calculates the trajectory for the ion, 
 * scattering off the given lattice
 *
 * @param ion - the ion to interact with the lattice
 * @param lattice - the lattice to scatter off
 * @param log - if true, will log the individual trajectory
 */
void traj(ion &ion, lattice &lattice, bool log);

#endif // ION_H_INCLUDED
