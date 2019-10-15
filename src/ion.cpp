#include "ion.h"
#include <math.h>
#include "potentials.h"
#include "safio.h"
#include "hameq.h"
#include "lattice.h"
#include "vec_math.h"

#include <functional> // std::minus 
#include <algorithm> // std::transform 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cstdio>

#define M_PI           3.14159265358979323846  /* pi */

double sqr(double *V)
{
    return V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
}

ion::~ion()
{
}

ion::ion()
{
    near = 0;
}

void ion::check_distances(double x, double y, double z, bool predicted)
{
    //r_min = 1e3;
    double ax, ay, az;

    for(int i = 0; i<near; i++)
    {
        site s = near_sites[i];

        ax = predicted? s.r_t[0] : s.r[0];
        ay = predicted? s.r_t[1] : s.r[1];
        az = predicted? s.r_t[2] : s.r[2];

        double dx = x - ax;
        double dy = y - ay;
        double dz = z - az;

        double r = sqrt(dx*dx + dy*dy + dz*dz);
        near_dists[i] = r;
        r_min = std::min(r_min, r);
    }
}

int ion::fill_nearest(lattice &lattice, int radius, int target_num)
{
    ion ion = *this;

    int cell_x = ion[0];
    int cell_y = ion[1];
    int cell_z = ion[2];

    int pos_hash = to_hash(cell_x, cell_y, cell_z);
    if(pos_hash == last_index)
        return near;
    last_index = pos_hash;

    int r = radius;
    int num;
    int near = 0;

    //This ensures we at least track the surface atoms.
    int k_max = std::min(cell_z + radius, radius);
    int k_min = k_max - 2*radius;

    int i_min = cell_x - r;
    int i_max = cell_x + r;

    int j_min = cell_y - r;
    int j_max = cell_y + r;

    for(int i = i_min; i<=i_max; i++)
    {
        for(int j = j_min; i<=j_max; i++)
        {
            for(int k = k_min; k<=k_max; k++)
            {
                pos_hash = to_hash(i, j, k);
                cell cel;

                if(lattice.cell_map.find(pos_hash) == lattice.cell_map.end())
                {
                    continue;
                }
                else
                {
                    cel = lattice.cell_map[pos_hash];
                }

                num = cel.num;

                for(int i = 0; i<num; i++)
                {
                    site s = cel.sites[i];

                    //Reset the site to original location.
                    s[0] = s.r_0[0];
                    s[1] = s.r_0[1];
                    s[2] = s.r_0[2];

                    //TODO set the momentum for the site
                    s.p[0] = 0;
                    s.p[1] = 0;
                    s.p[2] = 0;

                    near_sites[near] = s;
                    near_atoms[near] = s.atom.index;
                    near++;
                }
            }
        }
    }

    return near;
}


void ion::set_KE(double eV, double theta0, double phi0,double x, double y)
{
    atom.mass = settings.MASS;
    atom.symbol = settings.SYMION;
    //TODO lookup table for atomic symbols...
    //TODO charge config somewhere.
    q = 1;

    double p0 = sqrt(2 * atom.mass * eV);
    double p_trans = p0 * sin(theta0);

    //This is the initial momentum, before surface effects.
    double p_x0 = p_trans * cos(phi0);
    double p_y0 = p_trans * sin(phi0);
    double p_z0 = -p0 * cos(theta0);

    //Impact parameters offsets.
    r[0] = -settings.Z1 * tan(theta0) * cos(phi0) + x;
    r[1] = -settings.Z1 * tan(theta0) * sin(phi0) + y;
    r_0[2] = r[2] = settings.Z1;

    r_0[0] = x;
    r_0[1] = y;

    //If we have image effect, account for that here.
    if(settings.IMAGE)
    {
        p_z0 = -sqrt((p_z0 * p_z0) - (2 * atom.mass * Vi_z(settings.Z1, q)));
    }

    p[0] = p_x0;
    p[1] = p_y0;
    p[2] = p_z0;
}

bool validate(ion &ion, bool *buried, bool *off_edge, bool *stuck, bool *froze, bool *left, double *dt, double E)
{
    //Verify time step is in range.
    *dt = std::min(std::max(*dt, settings.DELLOW), settings.DELT0);

    //left crystal
    if(ion[2] > settings.Z1)
    {
        *left = true;
        return false;
    }
    //Fell of the edge
    if(ion[0] > 100 || ion[0] < -100 || ion[1] > 100 || ion[1] < -100)
    {
        *off_edge = true;
        return false;
    }
    //Buried
    if(ion[2] < -settings.BDIST)
    {
        *buried = true;
        return false;
    }
    //Too many steps
    if(ion.steps > 4000)
    {
        *froze = true;
        return false;
    }

    if((E + ion.v_total) < settings.SENRGY)
    {
        *stuck = true;
        return false;
    }
    return true;
}

void traj(ion &ion, lattice &lattice, bool log)
{
    //Time step
    double dt = 0.01;
    //Kinetic Energy of the ion
    double E;
    //Magnitude squared of the ion momentum
    double psq;
    //Ion mass
    double mass = ion.atom.mass;

    double pp = 0, pzz = 0, theta, phi;

    //Control conditions
    bool buried = false;
    bool froze = false;
    bool stuck = false;
    bool off_edge = false;
    bool left = false;

    //Parameters for checking
    //how much things are changing by.
    double V, V_t, dV;
    double change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = sqr(ion.p);
    E = psq * 0.5 / mass;
    ion.steps = 0;

start:

    //Find nearby lattice atoms
    ion.fill_nearest(lattice, 2, settings.NPART);
    //Increment counter for how many steps we have taken
    ion.steps++;

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, &dt, E);
    if(froze || buried || stuck || left || off_edge)
    {
        if(log)
            std::cout << "Exited" << std::endl;
        goto end;
    }

    //Find the forces at the current location
    run_hameq(ion, lattice, 0);
    V = ion.v_total;
    //Find the forces at the next location,
    run_hameq(ion, lattice, dt);
    V_t = ion.v_total;

    //Now we do some checks to see if the timestep needs to be adjusted
    dV = V_t - V;

    //Compute amount for time to step by.
    if(dV!=0)
    {
        dV = fabs(dV);
        //If we have exceeded max change, drop timestep.
        if(dV > settings.DEMAX && dt > settings.DELLOW)
        {
            dt = settings.DELLOW;
            //Reset potential counter to previous loctation
            ion.v_total = V;
            //Reset these, as it means we are passing into the surface.
            ion.last_index = -1;
            ion.near = 0;
            goto start;
        }
        if(V == 0)
        {
            //Reset these, as it means we are passing into the surface.
            ion.last_index = -1;
            ion.near = 0;
            //Set to low time step for entering this region
            change = 0.1;
        }
        else
        {
            //get ratio of change in f to f.
            dV /= V;
            //Change based on relative value of that to abserr option
            change = settings.ABSERR/dV;
        }

        if(change >= 2)
            change = 2;
        if(change > 1 && change < 2)
            change = 1;
        
        //Large change in energy, try to reduce timestep.
        if(change < .2 && dt > settings.DELLOW)
        {
            dt *= change;
            if(dt < settings.DELLOW)
                dt = settings.DELLOW;
            //Reset potential counter to previous loctation
            ion.v_total = V;
            goto start;
        }
    }
    else
    {
        change = 2;
    }

    //Apply changes, this updates energy and lattice loctations.
    apply_hameq(ion, lattice, dt);

    //Update some parameters for saving.
    ion.r_0[2] = fmin(ion.r[2], ion.r_0[2]);
    ion.time += dt;

    //Update kinetic energy.
    psq = sqr(ion.p);
    E = psq * 0.5 / mass;

    //Log things if needed
    if(log)
    {
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\t\%f\n",
                ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],
                ion.time,ion.steps,E,ion.v_total,(E+ion.v_total), ion.r_min,dt,dV);
        debug_file << buffer;
    }

    //Apply differences to timestep;
    dt *= change;

    //Back up we go.
    goto start;

end:

    if(froze)
    {
        theta = 0;
        phi = 90;
        E = -300;
    }
    else if(stuck)
    {
        theta = 0;
        phi = 90;
        E = -100;
    }
    else if(buried)
    {
        theta = 0;
        phi = 90;
        E = -200;
    }
    else if(off_edge)
    {
        theta = 0;
        phi = 90;
        E = -400;
    }
    else
    {
        double px = ion.p[0];
        double py = ion.p[1];
        double pz = ion.p[2];
        psq = sqr(ion.p);
        //Find the momentum at infinity
        if(settings.IMAGE)
        {
            pp = (pz*pz) - (2*mass*Vi_z(settings.Z1, ion.q));
            pzz = pp < 0 ? -sqrt(-pp) : sqrt(pp);
        }
        else
        {
            pp = pz * pz;
            pzz = pz;
        }
        //Ion is not escaping.
        if(pzz < 0)
        {
            theta = 0;
            phi = 90;
            E = -10;
        }
        else
        {
            E = 0.5 * psq / mass;
            theta = acos(pzz/sqrt(psq)) * 180/M_PI;
            if(px == 0 && py==0)
            {
                phi = 90;
            }
            else
            {
                phi = atan2(py, px) * 180/M_PI;
            }
            // sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%f\n",
            //         ion.r_0[0],ion.r_0[1],ion.r_0[2],
            //         E,theta,phi,ion.time,ion.steps,ion.max_n,ion.r_min);
            // out_file << buffer;
        }
    }
   sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%f\n",
           ion.r_0[0],ion.r_0[1],ion.r_0[2],
           E,theta,phi,ion.time,ion.steps,ion.max_n,ion.r_min);
   out_file << buffer;
    return;
}
