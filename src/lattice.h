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
    //Used to check if this cell has been checked
    int ion_stamp = -1;

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
    void add_site(Atom a, double px, double py, double pz);
    Cell* get_cell(double x, double y, double z);
    Cell* make_cell(double x, double y, double z);

};

#endif // LATTICE_H_INCLUDED
