#include "lattice.h"
#include "safio.h"
#include <cmath>
#include <cstdio>
#include <algorithm>    //std::sort
#include <vector>
#include "vec_math.h"

void lattice::build_lattice()
{
    //TODO add some way to define this in safio.
    int n = 20;
    // This adds some leeway to account for
    // floating point error in the matrix multiplications
    double zTop = settings.AZ * 0.1;
    double zBottom = -settings.AZ * 4;

    vec3d dir;
    //We like "Up"
    dir.set(0,0,1);
    vec3d axis;
    axis.set(settings.face);

    //Basis vectors in the basis coordinates
    vec3d ex_basis;
    ex_basis.set(1,0,0);
    vec3d ey_basis;
    ey_basis.set(0,1,0);
    vec3d ez_basis;
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
        site s;
        vec3d tmp;

        //Make basis of correct size.
        tmp.set(settings.BASIS[i].r);
        tmp[0] *= settings.AX;
        tmp[1] *= settings.AY;
        tmp[2] *= settings.AZ;
        vec3d v = R * tmp;
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

    vec3d cell_pos;
    cell_pos.set(0,0,0);

    int max_in_cell = 0;

    int ns = -n;
    int ne = n;
    double px, py, pz;

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
                    site old = basis[i];
                    pz = cell_pos[2] + old[2];

                    //Cut off bottom of the crystal at some point.
                    if(pz < zBottom)
                        continue;

                    px = cell_pos[0] + old[0];
                    py = cell_pos[1] + old[1];

                    //This hashes to the current cube
                    int pos_hash = to_hash(px, py, pz);

                    int *cel_num;
                    int num = 0;
                    site *cel_sites;

                    if(cell_map.find(pos_hash) == cell_map.end())
                    {
                        cell *cel = new cell();
                        cell_map[pos_hash] = cel;
                        cel->num = num;
                        cel_num = &(cel->num);
                        cel->pos_hash = pos_hash;
                        cel_sites = (cel->sites);
                    }
                    else
                    {
                        cel_num = &(cell_map[pos_hash]->num);
                        num = *cel_num;
                        cel_sites = (cell_map[pos_hash]->sites);
                    }

                    site s;
                    atom a;
                    s.r_0[0] = px;
                    s.r_0[1] = py;
                    s.r_0[2] = pz;
                    s.reset();
                    a = settings.ATOMS[old.index-1];
                    //Sites are indexed to size, so that they
                    //can be looked up to find their atom later.
                    s.index = sites.size();
                    s.atom = a;
                    sites.push_back(s);
                    cel_sites[num] = s;
                    num++;
                    *cel_num = num;
                    if(num > max_in_cell) max_in_cell = num;
                }
            }
        }
    }
    std::cout << "built lattice" << std::endl;
    debug_file << "built lattice" << std::endl;
}

void site::reset()
{
    double zeros[3] = {0,0,0};
    std::copy(std::begin(r_0), std::end(r_0), r);
    std::copy(std::begin(r_0), std::end(r_0), r_t);

    std::copy(std::begin(zeros), std::end(zeros), dp_dt);
    std::copy(std::begin(zeros), std::end(zeros), dp_dt_t);

    //TODO set the momentum for the site
    //base on thermal stuff
    std::copy(std::begin(zeros), std::end(zeros), p);
    std::copy(std::begin(zeros), std::end(zeros), p_t);
}

double site::distance(site &other, bool predicted)
{
    if(predicted)
    {
        return sqrt((r_t[0]-other.r_t[0])*(r_t[0]-other.r_t[0])+
                    (r_t[1]-other.r_t[1])*(r_t[1]-other.r_t[1])+
                    (r_t[2]-other.r_t[2])*(r_t[2]-other.r_t[2]));
    }
    else
    {
        return sqrt((r[0]-other.r[0])*(r[0]-other.r[0])+
                    (r[1]-other.r[1])*(r[1]-other.r[1])+
                    (r[2]-other.r[2])*(r[2]-other.r[2]));
    }
}

void site::write_info()
{
    debug_file <<"Atom: "<< atom.symbol << std::endl;
    debug_file <<"r  : "<< r[0] <<" "<< r[1] <<" "<< r[2] << std::endl;
    debug_file <<"p  : "<< p[0] <<" "<< p[1] <<" "<< p[2] << std::endl;
    debug_file <<"r_t: "<< r_t[0] <<" "<< r_t[1] <<" "<< r_t[2] << std::endl;
    debug_file <<"p_t: "<< p_t[0] <<" "<< p_t[1] <<" "<< p_t[2] << std::endl;
    debug_file <<"F: "<< dp_dt[0] <<" "<< dp_dt[1] <<" "<< dp_dt[2] << std::endl;
    debug_file <<"F_t: "<< dp_dt_t[0] <<" "<< dp_dt_t[1] <<" "<< dp_dt_t[2] << std::endl;
}

void lattice::reset()
{

}
