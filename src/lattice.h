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
    int num = 0;
    //Index for this in hash map
    int pos_hash = -1;
    //Used to check if this cell has been checked
    int check_stamp = -1;
    //Used to check if this cell has been checked
    int ion_stamp = -1;

    Cell()
    {
        sites = new Site[100];
    }

    Cell(const Cell& other);

    ~Cell()
    {
        delete []sites;
    }
};

struct Lattice
{
    //All sites in the lattice, 
    //This is not actually used much.
    std::vector<Site*> sites;

    //All cells in the lattice, This is
    //what is used for any lookups
    std::unordered_map<int,Cell> cell_map;

    //Rotation matrix for lattice
    Mat3d R;
    //Inverse of the rotation matrix
    Mat3d R_inv;

    //Basis vectors for lattice.
    Vec3d ex;
    Vec3d ey;
    Vec3d ez;

    //Default constructor
    Lattice(){}
    //Copy constructor
    Lattice(const Lattice& other);
    //Destructor
    ~Lattice();
    //Constructs the lattice based on the settings
    void build_lattice();
    //Adds an atom of type a, at location x, y, z;
    void add_site(Atom& a, double x, double y, double z);
    //Retrieves the cell for the given coordinates,
    //NULL if no cell is found
    Cell* get_cell(double x, double y, double z);
    //Version that uses the hash directly
    Cell* get_cell(int pos_hash);
    //Makes the cell for the given coordinate, if the
    //cell already exists, it retrieves old one instead.
    Cell* make_cell(double x, double y, double z);

};

#endif // LATTICE_H_INCLUDED
