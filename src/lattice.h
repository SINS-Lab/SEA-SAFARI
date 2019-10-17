#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED
#include <iostream>
#include "vec_math.h"
#include "safio.h"
#include <vector>
#include <unordered_map>

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
    std::vector<Site*> sites;
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

};

#endif // LATTICE_H_INCLUDED
