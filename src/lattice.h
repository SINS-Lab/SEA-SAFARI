#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED
#include <iostream>
#include "vec_math.h"
#include <vector>
#include <unordered_map>

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
    //momentum after dt
    double p_t[3];    
    //forces after dt
    double dp_dt_t[3];

    int last_ion = -1;

    Site()
    {
        reset();
    }

    double &operator[](int index)
    {
        return r[index];
    }

    double distance(Site &other, bool predicted);

    void reset();

    void write_info();

    int index;
};

class Cell
{
public:
    //Sites in this current cell
    Site *sites;
    //Number of sites in this cell
    int num;
    //Index for this in hash map
    int pos_hash;
    //Used to check if this cell has been checked
    int check_stamp = -1;

    Cell()
    {
        sites = new Site[100];
    }

    ~Cell()
    {
        delete []sites;
    }
};

struct Lattice
{
    //All sites in the lattice
    std::vector<Site> sites;
    //Basis used for this lattice,
    //Site indices for this are the originals
    std::vector<Site> basis;

    //All cells in the lattice
    std::unordered_map<int,Cell*> cell_map;

    //Rotation matrix for lattice
    Mat3d R;
    //Inverse of the rotation matrix
    Mat3d R_inv;

    //Basis vectors for lattice.
    Vec3d ex;
    Vec3d ey;
    Vec3d ez;

    void build_lattice();

    void reset();

};

#endif // LATTICE_H_INCLUDED
