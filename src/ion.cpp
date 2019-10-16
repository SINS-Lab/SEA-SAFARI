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

void ion::check_distances(double x, double y, double z, bool predicted)
{
    //TODO this should cull the sites if they are out of range?
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

    near = 0;

    double k_max;
    double k_min;

    double i_min;
    double i_max;

    double j_min;
    double j_max;

    i_min = cell_x - r;
    i_max = cell_x + r;
    j_min = cell_y - r;
    j_max = cell_y + r;

    //This ensures we at least track the surface atoms.
    k_min = std::min(cell_z - r, -r);
    k_max = cell_z + r;

    double near_distance_sq = 10*10;

    for(int i = i_min; i<=i_max; i++)
    {
        for(int j = j_min; j<=j_max; j++)
        {
            for(int k = k_min; k<=k_max; k++)
            {
                pos_hash = to_hash(i, j, k);
                
                site *cel_sites;
                num = 0;

                if(lattice.cell_map.find(pos_hash) == lattice.cell_map.end())
                {
                    continue;
                }
                else
                {
                    cell *cel = lattice.cell_map[pos_hash];
                    num = cel->num;
                    cel_sites = cel->sites;
                }

                for(int i = 0; i<num; i++)
                {
                    site s = cel_sites[i];
                    double rr = (s[0] - ion[0])*(s[0] - ion[0])
                               +(s[1] - ion[1])*(s[0] - ion[1])
                               +(s[2] - ion[2])*(s[0] - ion[2]);

                    if(rr > near_distance_sq) continue;

                    //Reset the site to original location.
                    s.reset();

                    near_sites[near] = s;
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
    r[2] = settings.Z1;

    r_0[0] = x;
    r_0[1] = y;
    r_0[2] = settings.Z1;

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

    if((E + ion.V) < settings.SENRGY)
    {
        *stuck = true;
        return false;
    }
    return true;
}

void traj(ion &ion, lattice &lattice, bool log)
{
    //Time step
    double dt = 0.1;
    //Energy of the ion
    double E, dE;
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
    //how much things are changing by.
    double T, T_t, dT;
    double change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = sqr(ion.p);
    E = psq * 0.5 / mass;
    ion.steps = 0;
    ion.V = 0;

    if(log)
    {
        debug_file << "\n\nStarting Ion Trajectory Output\n";
        ion.write_info();
        traj_file << "\n\nx\ty\tz\tpx\tpy\tpz\tt\tn\tT\tV\tE\tr_min\tdt\tdV\n";
    }

start:

    //Find nearby lattice atoms
    ion.fill_nearest(lattice, 1, settings.NPART);
    //Increment counter for how many steps we have taken
    ion.steps++;

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, &dt, E);

    if(froze || buried || stuck || left || off_edge)
    {
        if(log) debug_file << "Exited" << std::endl;
        goto end;
    }

    //Find the forces at the current location
    V = ion.V;
    run_hameq(ion, lattice, &T, 0);

    if(log)
    {
        debug_file << "\nBetween hameq" << dt << std::endl;
        ion.write_info();
    }

    //Find the forces at the next location,
    run_hameq(ion, lattice, &T_t, dt);
    V_t = ion.V;

    if(log)
    {
        debug_file << "\nAfter hameq " << dt << std::endl;
        ion.write_info();
    }

    //Now we do some checks to see if the timestep needs to be adjusted
    dV = V_t - V;
    dT = T_t - T;

    E = T + V;
    dE = std::max(fabs(dV), fabs(dT));

    //Compute amount for time to step by.
    if(dE != 0)
    {
        //check if energy has changed too much.
        change = settings.DEMAX/dE;
        
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
            ion.V = V;
            goto start;
        }
    }
    else
    {
        change = 2;
    }

    //Check if displacement has changed too much


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
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t\%f\n",
                ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],
                ion.time,ion.steps,E,ion.V,(E+ion.V), ion.near,dt,dV);
        traj_file << buffer;
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
