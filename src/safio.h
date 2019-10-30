#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "particles.h"

struct Safio
{
    //Parameters for the beam
    double E0;         //Beam Energy
    double THETA0;     //Beam Theta
    double PHI0;       //Crystal Phi
    double MASS;       //Ion Mass
    char* SYMION;      //Ion Atomic Symbol

    Atom ion;          //The Atom object for the ions

    //Parameters for the detector
    double EMIN;       //Min Energy
    double EMAX;       //Max Energy
    double ESIZE;      //Energy Resolution
    double ASIZE;      //Angular Resolution

    int NDTECT;        //Detector index
    double* DTECTPAR;  //Parameters for detector

    //Integration parameters
    double DELLOW;     //Low time step
    double DELT0;      //High time step
    double DEMAX;      //Exponent on error check
    double DEMIN;      //unused at present
    double ABSERR;     //Reference value for error check

    int NPART;         //Target interaction number
    bool RECOIL;       //Whether ions recoil on impact
    double Z1;         //Initial height of ions
    int MAX_STEPS;     //Maximum integration steps before failing

    double R_MAX;      //Maximum R value for interactions
    double rr_max;     // ^ squared.
    double DR_MIN_TAB; //Step for cached values

    double ZMIN;       //Start point for values caches
    double ZSTEP;      //Step for cached values

    int DIST_SEARCH;   //Unused
    int FAILED_DE;     //Unused

    int SCAT_FLAG;     //If this is 666, will run scats, otherwise tests
    int SCAT_TYPE;     // 666 = montecarlo, 777 = gridscat 888 = chainscat

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

    //initial index of ion, setting this allows same-trajectory thermal runs.
    int ion_index;

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
    //If true, will load the crystal from appropriate file name.
    bool load_crystal;
    //The surface face of the loaded crystal
    double* loaded_face;

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

    //Parameters for electronic frictional forces.
    double F_a = 0;
    double F_b = 0;
    
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
