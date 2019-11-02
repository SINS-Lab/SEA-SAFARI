#include "ion.h"
#include "potentials.h"
#include "safio.h"
#include "hameq.h"
#include "lattice.h"
#include "temps.h"
#include "space_math.h"

#include <functional> // std::minus 
#include <algorithm> // std::transform 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cstdio>
#include <cmath>

int Ion::fill_nearest(Lattice &lattice, int radius, int target_num, bool re_sort)
{
    Ion &ion = *this;

    //Initial locations are where ion is.
    double cell_x = ion.r[0];
    double cell_y = ion.r[1];
    double cell_z = ion.r[2];

    int pos_hash = to_hash(cell_x, cell_y, cell_z);
    //New site, so we need to re-calculate things
    if(pos_hash != last_index)
    {
        //We need to resort this if we find new index.
        re_sort = true;

        //Update last index checked.
        last_index = pos_hash;

        //Number in current cell being checked.
        int num;

        //Reset total nearby, as we are now looking for them.
        total_near = 0;

        //Initialize this big.
        ion.rr_min_find = 1e6;

        //Location of mask
        Vec3d loc;

        //Only particles closer than this are considered.
        double near_distance_sq = settings.rr_max;
        //Initialize this to max allowed distance.
        double rr_min = near_distance_sq;

        //volume of the mask to check.
        int nmax = pow(2*radius+1, 3);

        int last_check = -1;

        //Loop over the mask, this is a radial loop.
        for(int n = 0; n < nmax; n++)
        {
            //Gets x,y,z for the centered cube.
            index_to_loc(n, loc);
            //Translates to ion coordinate
            double x = loc[0] * CELL_SIZE + cell_x;
            double y = loc[1] * CELL_SIZE + cell_y;
            double z = loc[2] * CELL_SIZE + cell_z;

            pos_hash = to_hash(x,y,z);
            //Due to mask and cell sizes differing, this
            //can happen easily, so lets skip if it did.
            if(last_check == pos_hash) continue;
            last_check = pos_hash;

            Cell* cell = lattice.get_cell(x,y,z);
            //No cell here? skip.
            if(cell == NULL) continue;
            //Reset check stamp if new ion.
            if(cell->ion_stamp != ion.index) cell->check_stamp = -1;
            //Already seen this cell this tick? skip.
            if(cell->check_stamp == ion.steps) continue;
            //Stamp cell so it gets skipped if seen again.
            cell->check_stamp = ion.steps;
            cell->ion_stamp = ion.index;
            //Get number of sites from the cell
            num = cell->num;
            //Check each site in this cell.
            for(int i = 0; i<num; i++)
            {
                Site *s = &cell->sites[i];
                //If some other ion has seen the site, reset it here.
                //This reset puts it back to where it should be.
                if(s->last_ion!=ion.index)
                {
                    s->last_ion = ion.index;
                    s->reset();
                }
                //Check if site is close enough
                double rr = diff_sqr(ion.r, s->r);
                if(rr > near_distance_sq) continue;
                rr_min = std::min(rr_min, rr);
                ion.rr_min_find = std::min(rr, ion.rr_min_find);
                //Add the site to our tracked sites.
                near_sites[total_near] = s;
                total_near++;
                //If we have enough, goto end.
                if(total_near > 255) goto end;
            }
        }
    }
end:
    if(re_sort)
    {
        //Reset number nearby.
        near = 0;
        //Sort the sites, and update near to min of near or target_num
        std::sort(ion.near_sites, ion.near_sites + total_near, [this](Site* a, Site* b){return compare_sites(a, b);});
        //Update nearby number, taking min of these two.
        near = std::min(total_near, target_num);
    }
    //Sets this to 0, so that the max check later is fine.
    if(ion.rr_min_find == 1e6) ion.rr_min_find = 0;
    return near;
}

double zeros[3] = { 0,0,0 };

void Ion::set_KE(double E0,  double theta0, double phi0, double x, double y)
{
    //Start by resetting the ion to its initial state.
    reset();
    //Set the energy of the ion
    this->E0 = E0;
    //Attempt to thermalize the ion
    if(settings.ESIZE > 0)
    {
        thermaize_ion(*this);
    }

    atom = &settings.ion;
    //TODO lookup table for atomic symbols...
    //TODO charge config somewhere.
    q = 1;

    //Convert incoming angle to radians.
    theta0 = theta0 * M_PI / 180;
    phi0 = phi0 * M_PI / 180;

    double p0 = sqrt(2 * atom->mass * E0);
    double p_trans = p0 * sin(theta0);

    //This is the initial momentum, before surface effects.
    double p_x0 = p_trans * cos(phi0);
    double p_y0 = p_trans * sin(phi0);
    double p_z0 = -p0 * cos(theta0);

    //Impact parameters offsets, this aims at the impact point.
    r[0] = r_t[0] = -settings.Z1 * tan(theta0) * cos(phi0) + x;
    r[1] = r_t[1] = -settings.Z1 * tan(theta0) * sin(phi0) + y;
    r[2] = r_t[2] = settings.Z1;
    
    //Reset the forces, this cleans up some of the debug output.
    std::copy(zeros, zeros + 3, dp_dt);
    std::copy(zeros, zeros + 3, dp_dt_t);

    //Set the "initial" location to the targetted impact point.
    r_0[0] = x;
    r_0[1] = y;
    r_0[2] = settings.Z1;

    //If we have image effect, account for that here.
    if(settings.use_image)
    {
        p_z0 = -sqrt((p_z0 * p_z0) - (2 * atom->mass * Vi_z(settings.Z1, q)));
    }

    //Set the initial momentum of the ion
    p[0] = p_x0;
    p[1] = p_y0;
    p[2] = p_z0;
}

void Ion::reset()
{
    E0 = 0;
    steps = 0;
    time = 0;
    max_n = 0;
    near = 0;
    r_min = 1e3;
    rr_min_find = 1e3;
    last_index = -1;
    total_near = 0;
    V = 0;
}

void Site::reset()
{
    //Reset positions and momenta
    std::copy(r_0, r_0 + 3, r);
    std::copy(p_0, p_0 + 3, p);
    
    //Thermalize the site
    thermaize(*this);
}

void Site::write_info()
{
    debug_file << "Atom: " << atom->symbol << std::endl;
    debug_file << "r  : " << r[0] << " " << r[1] << " " << r[2] << std::endl;
    debug_file << "p  : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    debug_file << "r_t: " << r_t[0] << " " << r_t[1] << " " << r_t[2] << std::endl;
    debug_file << "F  : " << dp_dt[0] << " " << dp_dt[1] << " " << dp_dt[2] << std::endl;
    debug_file << "F_t: " << dp_dt_t[0] << " " << dp_dt_t[1] << " " << dp_dt_t[2] << std::endl;
}