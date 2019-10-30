#pragma once
#include <string>

struct Atom
{
    double mass = 1;
    double charge = 0;
    int index = -1;
    std::string symbol = "n";
    //Spring constants for x, y, z directions 
    double spring[3];
    //Std deviations in rx,ry,rz for the current temperature
    double dev_r[3];
    //Std deviations in px,py,pz for the current temperature
    double dev_p[3];
};

class Site
{
public:
    //Original location
    double *r_0 = NULL;
    //Original Momentum
    double *p_0 = NULL;
    //The atom here
    Atom* atom;

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
    //This is a unique identifier for this particle, the ion would be 0.
    int index = -1;
    //this is used for xyz output of nearest.
    bool near_check = 0;

    Site()
    {
        r_0 = new double[6];
        p_0 = r_0 + 3;
    }

    Site(const Site& other)
    {
        r_0 = other.r_0;
        p_0 = other.p_0;
        atom = other.atom;
        index = other.index;
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