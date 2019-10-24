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
    Site *sites = NULL;
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

    int id = 1;

    //Default constructor
    Lattice(){}
    //Copy constructor
    Lattice(const Lattice& other);
    //Destructor
    ~Lattice(){}
    //Constructs the lattice based on the settings
    void build_lattice();
    //Loads a lattice from the given input stream
    void load_lattice(std::ifstream& input);
    //Adds an atom of type a, at location x, y, z;
    void add_site(Atom* a, double x, double y, double z);
    //Adds the given site to the lattice
    void add_site(Site& site);
    //Retrieves the cell for the given coordinates,
    //NULL if no cell is found
    Cell* get_cell(double x, double y, double z);
    //Version that uses the hash directly
    Cell* get_cell(int pos_hash);
    //Makes the cell for the given coordinate, if the
    //cell already exists, it retrieves old one instead.
    Cell* make_cell(double x, double y, double z);

    /**
     * This function populates the contents of the given vector of sites, with
     * the rotated equivalent of the sites given.
     * 
     * @param dir - Reference direction, should generally be 0,0,1
     * @param face - the face of the crystal that should be pointing in dir
     * @param ex_basis - direction of first lattice coordinate reference
     * @param ey_basis - direction of second lattice coordinate reference
     * @param ez_basis - direction of third lattice coordinate reference
     * 
     * @param ex - new direction of first lattice reference
     * @param ey - new direction of second lattice reference
     * @param ez - new direction of third lattice reference
     * 
     * @param scale_basis - If true, the sites_in will be scaled by lattice parameters
     * 
     * @param sites_in - the vector to populate with rotated sites
     * @param sites_out - the vector to populate with rotated sites
     * @param maxZI - the index of the dir-most point in the new basis
     * 
     */ 
    void rotate_sites(Vec3d& dir, Vec3d& face, Vec3d& ex_basis, Vec3d& ey_basis, Vec3d& ez_basis,
                    Vec3d* ex, Vec3d* ey, Vec3d* ez, bool scale_basis,
                    std::vector<Site>* sites_out, std::vector<Site>& sites_in, int* maxZI);

};

#endif // LATTICE_H_INCLUDED
