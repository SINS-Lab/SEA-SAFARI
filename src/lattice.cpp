#include "lattice.h"
#include <cmath>
#include <cstdio>
#include <algorithm>    //std::sort
#include <vector>
#include "safio.h"
#include "space_math.h"
#include "string_utils.h"

void Lattice::rotate_sites(Vec3d& dir, Vec3d& face, Vec3d& ex_basis, Vec3d& ey_basis, Vec3d& ez_basis,
                Vec3d* ex, Vec3d* ey, Vec3d* ez, bool scale_basis,
                std::vector<Site>* sites_out, std::vector<Site>& sites_in, int* maxZI)
{
    //Rotation matrix for lattice
    Mat3d R;
    //Inverse of the rotation matrix
    Mat3d R_inv;

    //Rotation matrix and inverse.
    R = make_rot_matrix(dir, face);
    R_inv = R.invert();

    //Basis vectors in the rotated coordinate system
    *ex = R_inv * ex_basis * settings.AX;
    *ey = R_inv * ey_basis * settings.AY;
    *ez = R_inv * ez_basis * settings.AZ;

    double maxZ = -1e20;
    int num = sites_in.size();
    sites_out->resize(num);

    for(int i = 0; i<num; i++)
    {
        Site s;
        Site &old = sites_in[i];
        Vec3d tmp;
        //Make basis of correct size.
        tmp.set(old.r_0);
        if(scale_basis)
        {
            tmp[0] *= settings.AX;
            tmp[1] *= settings.AY;
            tmp[2] *= settings.AZ;
        }
        Vec3d v = R * tmp;
        s.r_0[0] = v[0];
        s.r_0[1] = v[1];
        s.r_0[2] = v[2];
        s.atom = old.atom;
        s.index = sites_in[i].index;
        (*sites_out)[i] = s;
        double dot = v*dir;
        if(dot > maxZ)
        {
            maxZ = dot;
            *maxZI = i;
        }
    }
}

void Lattice::build_lattice()
{
    //Use maximum, with some extra room, for the radius.
    int n = std::max(settings.RAX * 3, settings.RAY * 3);
    // This adds some leeway to account for
    // floating point error in the matrix multiplications
    double zTop = settings.AZ * 0.1;
    //Set the bottom of the slab to twice the buried distance
    double zBottom = -settings.BDIST * 2;
    //Basis for the lattice.
    std::vector<Site> basis;

    Vec3d dir;
    //We this is 001 if not specified
    dir.set(settings.loaded_face);
    Vec3d face;
    face.set(settings.face);

    //Basis vectors in the basis coordinates 
    //TODO configurable this
    Vec3d ex_basis;
    ex_basis.set(1,0,0);
    Vec3d ey_basis;
    ey_basis.set(0,1,0);
    Vec3d ez_basis;
    ez_basis.set(0,0,1);

    //Basis vectors for lattice.
    Vec3d ex;
    Vec3d ey;
    Vec3d ez;

    //Index of topmost site in the basis
    int maxZI = 0;

    //Populate the basis vectors
    rotate_sites(dir, face, ex_basis, ey_basis, ez_basis,
                &ex, &ey, &ez, true,
                &basis, settings.BASIS, &maxZI);

    Vec3d cell_pos;

    int ns = -n;
    int ne = n;
    double px, py, pz;

    double x_max = settings.AX * settings.RAX;
    double y_max = settings.AY * settings.RAY;

    for(int x = ns; x<= ne; x++)
    {
        for(int y = ns; y<= ne; y++)
        {
            for(int z = ns; z<= ne; z++)
            {
                cell_pos[2] = ez[0] * x + ez[1] * y + ez[2] * z;
                //Check if entire basis cell will fit
                if (cell_pos[2] + basis[maxZI].r_0[2] > zTop)
                    continue;

                cell_pos[0] = ex[0] * x + ex[1] * y + ex[2] * z;
                cell_pos[1] = ey[0] * x + ey[1] * y + ey[2] * z;

                for(int i = 0; i<settings.NBASIS; i++)
                {
                    Site old = basis[i];
                    pz = cell_pos[2] + old.r_0[2];

                    //Cut off bottom of the crystal at some point.
                    if(pz < zBottom)
                        continue;

                    px = cell_pos[0] + old.r_0[0];

                    //Out of bounds in x
                    if(px > x_max || px < -x_max)
                        continue;

                    py = cell_pos[1] + old.r_0[1];

                    //Out of bounds in y
                    if(py > y_max || py < -y_max)
                        continue;
                    
                    //This is the atom for this site.
                    Atom* a = &settings.ATOMS[old.index-1];
                    add_site(a, px, py, pz);
                }
            }
        }
    }
    std::cout << "built lattice" << std::endl;
    debug_file << "built lattice" << std::endl;
}


void Lattice::add_site(Site& s)
{
    //This makes the cell if it doesn't exist, otherwise gets old one.
    Cell *cell = make_cell(s.r_0[0], s.r_0[1], s.r_0[2]);

    //Get the values from the cell.
    int *cel_num = &(cell->num);
    int num = *cel_num;
    Site *cel_sites = (cell->sites);
    
    //Sites are indexed to size, so that they
    //can be looked up to find their atom later.
    s.index = sites.size();
    sites.push_back(&s);
    cel_sites[num] = s;
    *cel_num = num + 1;
    //Initializes the site
    s.reset();
}

Site* make_site(Atom* a, double px, double py, double pz)
{
    //Make the given site
    Site *s = new Site();
    //Initialize rest location
    s->r_0[0] = px;
    s->r_0[1] = py;
    s->r_0[2] = pz;

    //TODO instead do some thermal distribution.
    s->p_0[0] = 0;
    s->p_0[1] = 0;
    s->p_0[2] = 0;

    //Assign the atom for it
    s->atom = a;
    return s;
}

void Lattice::add_site(Atom* a, double px, double py, double pz)
{
    Site* s = make_site(a, px, py, pz);
    //Add the site to the lattice
    add_site(*s);
}

void Lattice::load_lattice(std::ifstream& input)
{
    std::string line;
    double* values;
    std::vector<Site> loaded_sites;
    std::cout << "Loading lattice" << std::endl;
    while(getline(input, line))
    {
        //Split line, and convert to array of doubles.
        //format is: x y z charge mass
        values = to_double_array(split(line), 0, 4);
        //Third column in the crys file is mass
        int charge = (int) values[3];
        Atom *atom;
        bool found = false;
        //Lookup the atom
        for(Atom &a: settings.ATOMS)
        {
            //Assume it is this one, TODO account for isotopes
            double diff = a.charge - charge;
            //Use this as the floating point errors might make
            //an == check fail.
            if(fabs(diff) < 1)
            {
                atom = &a;
                found = true;
                break;
            }
        }
        //If somehow charge didn't match
        if(!found)
        {
            //Assume it is default atom, and log an error.
            debug_file<<"No Atom found for charge "<<charge;
            //TODO maybe make a new atom for it instead?
            atom = &settings.ATOMS[0];
            debug_file << " setting to: "<<atom->symbol << ", " << atom->charge<<std::endl;
        }
        //Make a new site
        Site* s = make_site(atom, values[0], values[1], values[2]);
        //Add it to our list
        loaded_sites.push_back(*s);
        //Cleanup the values array.
        delete []values;
    }

    Vec3d dir;
    //We this is 001 if not specified
    dir.set(settings.loaded_face);
    Vec3d face;
    face.set(settings.face);

    //Basis vectors in the basis coordinates 
    //TODO configurable this
    Vec3d ex_basis;
    ex_basis.set(1,0,0);
    Vec3d ey_basis;
    ey_basis.set(0,1,0);
    Vec3d ez_basis;
    ez_basis.set(0,0,1);

    //Basis vectors for lattice.
    Vec3d ex;
    Vec3d ey;
    Vec3d ez;

    //Index of topmost site in the basis
    int maxZI = 0;

    std::vector<Site> processed_sites;

    if(dir * face.normalize() != 1)
    {
        std::cout << "Rotating Lattice" << std::endl;
        //Populate the rotated vectors
        rotate_sites(dir, face, ex_basis, ey_basis, ez_basis,
                    &ex, &ey, &ez, false,
                    &processed_sites, loaded_sites, &maxZI);
    }
    else
    {
        std::cout << "No need to Rotate Lattice" << std::endl;
        std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
        std::cout << face[0] << " " << face[1] << " " << face[2] << std::endl;
        processed_sites = loaded_sites;
    }

    int num = processed_sites.size();
    for(int i = 0; i<num; i++)
    {
        //We want a copy to add, rather than the original
        Site site = processed_sites[i];
        //Add the site
        add_site(&settings.ATOMS[site.atom->index-1], 
                 site.r_0[0], site.r_0[1], site.r_0[2]);
    }
    debug_file<< "Loaded " << sites.size() << " sites from file" << std::endl;
}

Cell* Lattice::get_cell(double x, double y, double z)
{
    int pos_hash = to_hash(x,y,z);
    return get_cell(pos_hash);
}

Cell* Lattice::get_cell(int pos_hash)
{
    if(cell_map.find(pos_hash) == cell_map.end()) return NULL;
    return &cell_map[pos_hash];
}

Cell* Lattice::make_cell(double x, double y, double z)
{
    Cell* cell = get_cell(x,y,z);
    //Return the old cell we had.
    if(cell != NULL) return cell;
    //Otherwise, make a new one.
    int pos_hash = to_hash(x,y,z);
    cell = new Cell();
    cell_map[pos_hash] = *cell;
    cell->num = 0;
    cell->pos_hash = pos_hash;
    cell->sites = new Site[100];
    return cell;
}

Lattice::Lattice(const Lattice& other)
{
    //This is not a deep copy, we don't care though.
    //sites = other.sites; so we don't even bother copying it.

    //This should be a deep copy, so long as cells copy correctly.
    cell_map = other.cell_map;
}

Cell::Cell(const Cell& other)
{
    num = other.num;
    //The copy knows how many it needs to store!
    sites = new Site[num];
    for (int i = 0; i < num; i++)
    {
        Site &original = other.sites[i];
        //Copy it over
        Site site = original;
        site.atom = original.atom;
        //Deallocate these two
        delete []site.r_0;
        delete []site.p_0;
        site.r_0 = original.r_0;
        site.p_0 = original.p_0;
        site.index = original.index;
        sites[i] = site;
    }
}