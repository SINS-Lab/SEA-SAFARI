#pragma once
#include "vec_math.h" // for diff_sqr
#include <string>     // for the atom symbol

#define MAX_NEAR 2048

class Atom
{
public:
    double mass = 1;
    double mass_inv = 1;
    double two_mass = 2;
    double mass_inv_2 = 0.5;
    double charge = 0;
    int index = -1;
    std::string symbol = "n";
    // Spring constants for x, y, z directions
    double spring[3];
    // Std deviations in rx,ry,rz for the current temperature
    double dev_r[3];
    // Std deviations in px,py,pz for the current temperature
    double dev_p[3];

    void init(double mass, double charge, std::string symbol)
    {
        this->mass = mass;
        this->mass_inv = 1 / mass;
        this->two_mass = mass * 2;
        this->mass_inv_2 = 1 / two_mass;
        this->charge = charge;
        this->symbol = symbol;
    }

    void init_pots(std::string &filename);
};

class Site
{
public:
    // Original location
    double *r_0 = NULL;
    // Original Momentum
    double *p_0 = NULL;
    // The atom here
    Atom *atom = NULL;
    // Charge for image considerations
    int q = 0;

    // position
    double r[3];
    // momentum
    double p[3];
    // forces
    double dp_dt[3];

    // position after time dt
    double r_t[3];
    // momentum after dt
    double p_t[3];
    // forces after dt
    double dp_dt_t[3];

    // Position during last neighbour check
    double r_u[3];

    // potential the particle is in.
    double V = 0;
    // Kinetic energy of the particle
    double T = 0;

    // potential the particle is in, at next time step
    double V_t = 0;
    // Kinetic energy of the particle, at next time step
    double T_t = 0;

    // Things here for saving to data files

    // Number of traj steps done, note that this is not
    // the same as the number of time steps, as this is
    // incremented even if it reduces timestep and tries over.
    int steps = 0;
    // Total flight time of the ion
    double time = 0;
    // This is a weighting factor for the ion's detectability.
    double weight = 1;
    // Maximum integration error for this ion.
    double Eerr_max = 0;
    // Number of times a site-site intersection occurs.
    int site_site_intersects = 0;

    // Used to track the last ion which has interacted with us.
    int last_ion = -1;
    // This is used for randomly position the site, normally
    // this is the same value as last_ion, but adaptive grid differs
    int thermal_seed = -1;
    // This is a unique identifier for this particle,
    // The lattice and the ions use different counters.
    int index = -1;
    // this is used for xyz output of nearest.
    bool near_check = 0;
    // This is last integration step performed, 
    // if this is the same for a lattice site and an ion, 
    // that ion will not be considered for force calculations
    int last_step = -1;
    // This is the last tick that the particle positions were updated.
    // This is used only for the more complex updating when considering
    // lattice-lattice correlations.
    int last_update = -1;

    // Start of block of values related to nearby sites

    // All of the sites nearby, only guarenteed to be filled
    // with site up to near
    Site **near_sites = NULL;
    // This is the nearby sites when at rest, this is
    // set during the initial setting of the springs.
    Site **rest_near_sites = NULL;
    // Number nearby at rest
    int rest_near_count = 0;
    // Maximum number of nearby atoms for this ion
    int max_n = 0;
    // Nearest distance ever
    double r_min = 1e3;
    // Nearest distance when finding nearby, squared
    double rr_min_find = 1e3;
    // This is the total number of Site*s stored in near_sites
    int total_near = 0;
    // Current number of nearby atoms.
    int near = 0;
    // This is the last place the near things were updated.
    // Setting this back to -1 will force a re-check of near.
    int last_index = -1;

    // These are for what cell we are in, they are not used for ions
    int cell_number = -1;
    // This is where in the cell's array we are.
    int cell_index = -1;

    // Whether we have left our original location
    bool left_origin = false;
    bool unbound = false;
    int sputter_tick = -1;


    // This is invalidated if we are turned into a projectile via
    // a cascade event, it is set to the ion index that invalidated us
    int valid = -1;

    // For the ion, this is min depth z, for the sites, it is rest position
    double log_z = 1e3;

    // This is used to flagging whether the lattice site should check nearby forces
    // For runs with multiple ions, each ion should have the same value for this.
    int hameq_tick = -1;
    // This is used simlarly to hameq_tick, however is only incremented before the
    // first call to hameq in the sequence
    int force_reset_tick = -1;

    // Some arrays for stuffing positions and momenta in for easier processing

    // This is stuffed with the position pointers for
    // the nearby sites's current and future positions
    // double **near_dists = NULL;
    // double **near_forces = NULL;

    Site()
    {
        r_0 = new double[6];
        p_0 = r_0 + 3;
    }

    ~Site()
    {
        if (near_sites != NULL)
            delete near_sites;
        if (rest_near_sites != NULL)
            delete rest_near_sites;
        // if (near_dists != NULL)
        //     delete near_dists;
        // if (near_forces != NULL)
        //     delete near_forces;
        near_sites = NULL;
        rest_near_sites = NULL;
        // near_dists = NULL;
        // near_forces = NULL;
    }

    Site(const Site &other)
    {
        r_0 = other.r_0;
        p_0 = other.p_0;
        atom = other.atom;
        index = other.index;
        cell_index = other.cell_index;
        cell_number = other.cell_number;
    }

    bool compare_sites(const Site *a, const Site *b)
    {
        double dist_a = diff_sqr(a->r, this->r);
        double dist_b = diff_sqr(b->r, this->r);
        return dist_a < dist_b;
    }

    void invalidate(int index)
    {
        valid = index;
    }

    void validate()
    {
        valid = -1;
    }

    /**
     * Resets the site to r_0 and p_0.
     * This also then calls thermalize, so
     * if an ion index is assigned previously,
     * that willl be used for seeding the RNG
     */
    void reset();

    void reset_forces()
    {
        dp_dt[0] = 0;
        dp_dt[1] = 0;
        dp_dt[2] = 0;
        V = 0;

        dp_dt_t[0] = 0;
        dp_dt_t[1] = 0;
        dp_dt_t[2] = 0;
        V_t = 0;
    }

    /**
     * Prints some info for this site to debug_file
     */
    void write_info();
};