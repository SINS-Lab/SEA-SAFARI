#ifndef SAFIO_H_INCLUDED
#define SAFIO_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>


struct Atom
{
    double mass;
    double charge;
    int index;
    std::string symbol;
    double spring[3];
};

class Site
{
public:
    //Original location
    double r_0[3];
    //Original Momentum
    double p_0[3];
    //The atom here
    Atom atom;

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

    //Used to track the last ion which has interacted with us.
    int last_ion = -1;
    //This is a unique identifier for this particle
    int index;

    Site()
    {
        reset();
    }

    double &operator[](int index)
    {
        return r[index];
    }

    void reset();

    void write_info();
};

struct Safio
{
    //Parameters for the beam
    double E0;
    double THETA0;
    double PHI0;
    double MASS;
    char* SYMION;

    //Parameters for the detector
    double EMIN;
    double EMAX;
    double ESIZE;
    double ASIZE;

    int NDTECT;
    double* DTECTPAR;

    //Integration parameters
    double DELLOW;
    double DELT0;
    double DEMAX;
    double DEMIN;
    double ABSERR;
    int NPART;
    bool RECOIL;
    double Z1;
    int MAX_STEPS;

    double R_MAX; //Maximum R value for interactions
    double rr_max;// ^ squared.
    double DR_MIN_TAB; //Step for cached values

    double ZMIN; //Start point for values caches
    double ZSTEP; //Step for cached values

    int MAXDIV; //Unused
    int MINDIV; //Unused

    int NWRITX;
    int NWRITY;

    //This is how many AX and AY
    //to build the lattice for,
    //The "radius" of the slab
    double RAY;
    double RAX;

    //Site Potential Values
    int NPAR;
    int IPOT;
    double* POTPAR;

    //Image Potential Values
    int NIMPAR;
    int IIMPOT;
    double* PIMPAR;

    double TEMP;
    double SEED;

    int NITER;//Unused

    //Whether to use image potentials
    bool IMAGE;

    //Energy to consider "stuck"
    double SENRGY;
    //Distance to consider "buried"
    double BDIST;

    //Parameters for the target
    //Lattice Constants
    double AX;
    double AY;
    double AZ;

    //Basis Coordinates
    int NBASIS;
    std::vector<Site> BASIS;

    //Surface face
    double* face;

    //Basis Atoms
    int NTYPES;
    std::vector<Atom> ATOMS;
    //Whether the lattice sites are springy
    bool CORR;
    double ATOMK;
    double RNEIGH;

    //Number of trajectories to run.
    //If this is 1, it will run a single gridscat,
    //located at XSTART, YSTART
    int NUMCHA;
    //X-range for impact parameters
    double XSTART;
    double XSTEP;
    double XSTOP;
    //Y-range for impact parameters
    double YSTART;
    double YSTEP;
    double YSTOP;
    
    //These were ZBL-potential related in the old 
    //version of safari, not currently used.
    int NBZ;
    double TOL;
    int NBG;
    double GTOL;
    double* GMAX;
    double* NG;
    double* NZ;
    double* ZMAX;

    void load(std::string safio_file);
};

//Global settings variable, we only need this loaded once anyway.
extern Safio settings;
//This is the .data file to store the outputs.
extern std::ofstream out_file;
//This is the .traj file to ion trajectory to.
extern std::ofstream traj_file;
//This is the .traj file to ion xyz to.
extern std::ofstream xyz_file;
//This is the .dbug file to store extra info.
extern std::ofstream debug_file;
//This is the .dbug file to lattice info in
extern std::ofstream crystal_file;

#endif // SAFIO_H_INCLUDED
