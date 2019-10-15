#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED
#include <iostream>
#include "vec_math.h"
#include <vector>
#include <unordered_map>

struct site
{
    //site position
    double r[3];
    //Original site location
    double r_0[3];
    //site momentum
    double p[3];    
    //Forces on the ion here
    double dp_dt[3];
    //Forces on the ion next step
    double dp_dt_t[3];

    site()
    {
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
    }

    double &operator[](int index)
    {
        return r[index];
    }
    int index;
};

class cell
{
public:
    //Sites in this current cell
    site sites[100];
    //Number of sites in this cell
    int num;
    //Index for this in hash map
    int pos_hash;
};

struct atom
{
    double mass;
    double charge;
    int index;
    std::string symbol;
    double spring[3];
};

struct lattice
{
    //All sites in the lattice
    std::vector<site> sites;
    //Atoms for the sites,
    // site index points to here.
    std::vector<atom> atoms;
    //Basis used for this lattice,
    //Site indices for this are the originals
    std::vector<site> basis;

    //All cells in the lattice
    std::unordered_map<int,cell> cell_map;

    //Rotation matrix for lattice
    mat3d R;
    //Inverse of the rotation matrix
    mat3d R_inv;

    //Basis vectors for lattice.
    vec3d ex;
    vec3d ey;
    vec3d ez;

    void build_lattice();

    void reset();

};

#endif // LATTICE_H_INCLUDED
