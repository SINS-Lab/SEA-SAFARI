#pragma once
#include "vec_math.h"
#include <string>

struct Atom
{
    double mass = 1;
    double charge = 0;
    int index = -1;
    std::string symbol = "n";
    // Spring constants for x, y, z directions
    double spring[3];
    // Std deviations in rx,ry,rz for the current temperature
    double dev_r[3];
    // Std deviations in px,py,pz for the current temperature
    double dev_p[3];

    //Lennard Jones potential related
    /**
     * These are pairs of parameters, 
     * for each sub-array,
     * index 0 is ε, index 1 is σ
     * 
     * The index in the main array is the atom index,
     * ie L_J_params[atom.index] is for a like pair.
     * 
     * TODO actually implement it this way later.
     */
    //double** L_J_params;
};

class Site
{
public:
    //Original location
    double *r_0 = NULL;
    //Original Momentum
    double *p_0 = NULL;
    //The atom here
    Atom *atom = NULL;

    //position
    double r[3];
    //momentum
    double p[3];
    //forces
    double dp_dt[3];

    //position after time dt
    double r_t[3];
    //forces after dt
    double dp_dt_t[3];

    // Position during last neighbour
    double r_u[3];

    //Used to track the last ion which has interacted with us.
    int last_ion = -1;
    //This is used for randomly position the site, normally
    //this is the same value as last_ion, but adaptive grid differs
    int thermal_seed = -1;
    //This is a unique identifier for this particle,
    //The lattice and the ions use different counters.
    int index = -1;
    //this is used for xyz output of nearest.
    bool near_check = 0;
    //This is last integration step performed, this
    //allows for faster checks of whether force needs
    //to be reset.
    int last_step = -1;
    //This is the last tick that the particle positions were updated.
    //This is used only for the more complex updating when considering
    //lattice-lattice correlations.
    int last_update = -1;

    //Start of block of values related to nearby sites

    //All of the sites nearby, only guarenteed to be filled
    //with site up to near
    Site **near_sites = NULL;
    // This is the nearby sites when at rest, this is
    // set during the initial setting of the springs.
    Site **rest_near_sites = NULL;
    int rest_near_count = 0;
    //Maximum number of nearby atoms for this ion
    int max_n = 0;
    //Nearest distance ever
    double r_min = 1e3;
    //Nearest distance when finding nearby, squared
    double rr_min_find = 1e3;
    //This is the total number of Site*s stored in near_sites
    int total_near = 0;
    //Current number of nearby atoms.
    int near = 0;
    //This is the last place the near things were updated.
    //Setting this back to -1 will force a re-check of near.
    int last_index = -1;

    //These are for what cell we are in, they are not used for ions
    int cell_number = -1;
    //This is where in the cell's array we are.
    int cell_index = -1;

    uint64_t hameq_tick = 0;
    uint64_t update_tick = 0;

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
        near_sites = NULL;
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

    /**
     * Resets the site to r_0 and p_0.
     * This also then calls thermalize, so
     * if an ion index is assigned previously,
     * that willl be used for seeding the RNG
     */
    void reset();

    /**
     * Prints some info for this site to debug_file
     */
    void write_info();
};