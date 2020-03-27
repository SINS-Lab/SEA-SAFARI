#pragma once
#include <iostream>
#include "space_math.h" // vec_math.h and mask stuff
#include "particles.h"
#include <vector>
#include <unordered_map>

class Cell
{
    int lastSize = -1;

public:
    //Sites in this current cell
    Site **sites = NULL;
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
    }

    Cell(const Cell &other);

    ~Cell()
    {
        delete[] sites;
    }

    void addSite(Site *site);

    void removeSite(Site *site);
};

struct Lattice
{
    //All sites in the lattice,
    //This is not actually used much.
    std::vector<Site *> sites;

    //All cells in the lattice, This is
    //what is used for any lookups
    std::unordered_map<int, Cell *> cell_map;

    int id = 1;

    //If relevant, this is used to determine
    //the range of locations used for the
    //xyz trajectories, the 6 entries are:
    //x1,y1,z1,x2,y2,z2
    //This is relevant if settings.SCAT_TYPE
    double xyz_bounds[6];

    Mask mask;

    //Here we have some values which are recored as the run occurs.
    //These are recorded to the debug file after the run is completed.
    int undetectable_num = 0; //Code -5
    int trapped_num = 0;      //Code -10
    int stuck_num = 0;        //Code -100
    int buried_num = 0;       //Code -200
    int froze_num = 0;        //Code -300
    int left_num = 0;         //Code -400
    int err_num = 0;          //Code -500
    int out_of_mask = 0;      //No code for this.

    //Default constructor
    Lattice() {}
    //Copy constructor
    Lattice(const Lattice &other);
    //Destructor
    ~Lattice() {}
    //Constructs the lattice based on the settings
    void build_lattice();
    //Loads a lattice from the given input stream
    void load_lattice(std::ifstream &input);
    //Adds an atom of type a, at location x, y, z;
    void add_site(Atom *a, double x, double y, double z);
    //Adds the given site to the lattice
    void add_site(Site *site);
    //Retrieves the cell for the given coordinates,
    //NULL if no cell is found
    Cell *get_cell(double x, double y, double z);
    //Version that uses the hash directly
    Cell *get_cell(int pos_hash);
    //Makes the cell for the given coordinate, if the
    //cell already exists, it retrieves old one instead.
    Cell *make_cell(double x, double y, double z);
    //Version that uses the hash directly
    Cell *make_cell(int pos_hash);

    void clear_stats()
    {
        undetectable_num = 0; //Code -5
        trapped_num = 0;      //Code -10
        stuck_num = 0;        //Code -100
        buried_num = 0;       //Code -200
        froze_num = 0;        //Code -300
        left_num = 0;         //Code -400
        err_num = 0;          //Code -500
    }

    void add_stats(Lattice *other)
    {
        undetectable_num += other->undetectable_num;
        trapped_num += other->trapped_num;
        stuck_num += other->stuck_num;
        buried_num += other->buried_num;
        froze_num += other->froze_num;
        left_num += other->left_num;
        err_num += other->err_num;
        out_of_mask += other->out_of_mask;
    }

    /**
     * This initializes the nearby neighbours for the sites.
     * It should call traj.h's fill_nearest for the sites, and
     * then do some cleanup.
     * 
     * @param nearest - Number of nearest neighbours to find.
     */
    void init_springs(int nearest);

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
    void rotate_sites(Vec3d &dir, Vec3d &face, Vec3d &ex_basis, Vec3d &ey_basis, Vec3d &ez_basis,
                      Vec3d *ex, Vec3d *ey, Vec3d *ez, bool scale_basis,
                      std::vector<Site> *sites_out, std::vector<Site> &sites_in, int *maxZI);
};

void moveSite(Site *site, Cell *from, Cell *to);