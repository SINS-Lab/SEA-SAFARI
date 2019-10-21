#include "lattice.h"
#include <cmath>
#include <cstdio>
#include <algorithm>    //std::sort
#include <vector>
#include "space_math.h"

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
    //We like "Up"
    dir.set(0,0,1);
    Vec3d axis;
    axis.set(settings.face);

    //Basis vectors in the basis coordinates
    Vec3d ex_basis;
    ex_basis.set(1,0,0);
    Vec3d ey_basis;
    ey_basis.set(0,1,0);
    Vec3d ez_basis;
    ez_basis.set(0,0,1);

    //Rotation matrix and inverse.
    R = make_rot_matrix(dir, axis);
    R_inv = R.invert();

    //Basis vectors in the rotated coordinate system
    ex = R_inv * ex_basis * settings.AX;
    ey = R_inv * ey_basis * settings.AY;
    ez = R_inv * ez_basis * settings.AZ;

    double maxZ = -1e20;
    int maxZI = 0;

    for(int i = 0; i<settings.NBASIS; i++)
    {
        Site s;
        Vec3d tmp;

        //Make basis of correct size.
        tmp.set(settings.BASIS[i].r);
        tmp[0] *= settings.AX;
        tmp[1] *= settings.AY;
        tmp[2] *= settings.AZ;
        Vec3d v = R * tmp;
        s[0] = v[0];
        s[1] = v[1];
        s[2] = v[2];
        s.index = settings.BASIS[i].index;
        basis.push_back(s);

        if(s[2] > maxZ)
        {
            maxZ = s[2];
            maxZI = i;
        }
    }

    Vec3d cell_pos;
    cell_pos.set(0,0,0);

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
                if (cell_pos[2] + basis[maxZI][2] > zTop)
                    continue;

                cell_pos[0] = ex[0] * x + ex[1] * y + ex[2] * z;
                cell_pos[1] = ey[0] * x + ey[1] * y + ey[2] * z;

                for(int i = 0; i<settings.NBASIS; i++)
                {
                    Site old = basis[i];
                    pz = cell_pos[2] + old[2];

                    //Cut off bottom of the crystal at some point.
                    if(pz < zBottom)
                        continue;

                    px = cell_pos[0] + old[0];

                    //Out of bounds in x
                    if(px > x_max || px < -x_max)
                        continue;

                    py = cell_pos[1] + old[1];

                    //Out of bounds in y
                    if(py > y_max || py < -y_max)
                        continue;
                    
                    //This is the atom for this site.
                    Atom &a = settings.ATOMS[old.index-1];
                    add_site(a, px, py, pz);
                }
            }
        }
    }
    std::cout << "built lattice" << std::endl;
    debug_file << "built lattice" << std::endl;
}

void Lattice::add_site(Atom& a, double px, double py, double pz)
{
    //This makes the cell if it doesn't exist, otherwise gets old one.
    Cell *cell = make_cell(px, py, pz);

    //Get the values from the cell.
    int *cel_num = &(cell->num);
    int num = *cel_num;
    Site *cel_sites = (cell->sites);
    
    Site *s = new Site();
    s->r_0[0] = px;
    s->r_0[1] = py;
    s->r_0[2] = pz;

    //TODO instead do some thermal distribution.
    s->p_0[0] = 0;
    s->p_0[1] = 0;
    s->p_0[2] = 0;

    //Initializes the site
    s->reset();

    //Sites are indexed to size, so that they
    //can be looked up to find their atom later.
    s->index = sites.size();
    s->atom = &a;
    sites.push_back(s);
    cel_sites[num] = *s;
    *cel_num = num + 1;
}

Cell* Lattice::get_cell(double x, double y, double z)
{
    int pos_hash = to_hash(x,y,z);
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

Lattice::~Lattice()
{
    for (auto x : cell_map)
    {
        delete &x.second;
    }
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
        site.r_0 = original.r_0;
        site.p_0 = original.p_0;
        sites[i] = site;
    }
}