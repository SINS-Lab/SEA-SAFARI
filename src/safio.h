#ifndef SAFIO_H_INCLUDED
#define SAFIO_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "lattice.h"

struct Safio
{
    //Parameters for the beam
    double E0;
    double THETA0;
    double PHI0;
    double MASS;
    std::string SYMION;

    //Parameters for the detector
    double EMIN;
    double EMAX;
    double ESIZE;
    double ASIZE;

    int NDTECT;
    std::vector<double>DTECTPAR;

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
    double RRMIN;
    double RRSTEP;
    double ZMIN;
    double ZSTEP;
    int MAXDIV;
    int MINDIV;

    int NWRITX;
    int NWRITY;
    double FAX;
    double FAY;
    int NPAR;
    int IPOT;
    std::vector<double>POTPAR;
    int NIMPAR;
    int IIMPOT;
    std::vector<double>PIMPAR;
    double TEMP;
    double SEED;
    int NITER;
    bool IMAGE;
    double SENRGY;
    double BDIST;

    //Parameters for the target
    //Lattice Constants
    double AX;
    double AY;
    double AZ;

    //Basis Coordinates
    int NBASIS;
    std::vector<Site> BASIS;
    int NTYPES;

    double face[3];

    //Basis Atoms
    std::vector<Atom> ATOMS;
    bool CORR;
    double ATOMK;
    double RNEIGH;

    //No idea on good defaults for these;
    int NUMCHA;
    double XSTART;
    double XSTEP;
    double XSTOP;
    double YSTART;
    double YSTEP;
    double YSTOP;
    int NBZ;
    double TOL;
    int NBG;
    double GTOL;
    std::vector<double>GMAX;
    std::vector<double>NG;
    std::vector<double>NZ;
    std::vector<double>ZMAX;

    void load();
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
