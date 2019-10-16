#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED
#include <iostream>
#include "vec_math.h"
#include <vector>
#include <unordered_map>

struct atom
{
    double mass;
    double charge;
    int index;
    std::string symbol;
    double spring[3];
};

class site
{
public:
    //Original location
    double r_0[3];
    //The atom here
    atom atom;

    //position
    double r[3];
    //momentum
    double p[3];    
    //forces
    double dp_dt[3];

    //position after time dt
    double r_t[3];
    //momentum after dt
    double p_t[3];    
    //forces after dt
    double dp_dt_t[3];

    site()
    {
        reset();
    }

    double &operator[](int index)
    {
        return r[index];
    }

    double distance(site &other, bool predicted);

    void reset();

    void write_info();

    int index;
};

class cell
{
public:
    //Sites in this current cell
    site *sites;
    //Number of sites in this cell
    int num;
    //Index for this in hash map
    int pos_hash;

    cell()
    {
        sites = new site[100];
    }

    ~cell()
    {
        delete []sites;
    }
};

struct lattice
{
    //All sites in the lattice
    std::vector<site> sites;
    //Basis used for this lattice,
    //Site indices for this are the originals
    std::vector<site> basis;

    //All cells in the lattice
    std::unordered_map<int,cell*> cell_map;

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
