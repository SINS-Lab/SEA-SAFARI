#include "ion.h"
#include <math.h>
#include "potentials.h"
#include "safio.h"
#include "hameq.h"
#include "lattice.h"
#include "vec_math.h"
#include "space_math.h"

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

double diff_sqr(double *V, double *Y)
{
    double dx = V[0] - Y[0];
    double dy = V[1] - Y[1];
    double dz = V[2] - Y[2];
    return dx*dx + dy*dy + dz*dz;
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

    vec3d loc;

    double near_distance_sq = 5*5;
    double rr_min = near_distance_sq;

    int nmax = pow(2*r+1, 3);

    //Centre the cell on the surface if it is above it.
    cell_z = std::min(cell_z, 0);

    for(int n = 0; n < nmax; n++)
    {
        index_to_loc(n, loc);
        double x = loc[0] + cell_x;
        double y = loc[1] + cell_y;
        double z = loc[2] + cell_z;

        pos_hash = to_hash(x, y, z);
        
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
            double rr = diff_sqr(ion.r, s.r);
            if(rr > near_distance_sq) continue;
            rr_min = std::min(rr_min, rr);
            near_sites[near] = s;
            near++;
            if(near > 90) goto end;
        }
    }
end:
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

bool validate(ion &ion, bool *buried, bool *off_edge, bool *stuck, bool *froze, bool *left, double E)
{
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

    //Potential Energy
    double V, dV;
    //Kinetic Energy
    double T;
    
    //Parameters for checking if things need re-calculating
    vec3d r;
    vec3d dr;
    double diff;
    double dxp, dyp, dzp;
    double dpx, dpy, dpz;
    double dr_max;

    //Distance in angstroms to consider far enough moved.
    //After it moves this far, it will re-calculate nearest.
    double r_reset = 0;

    //Multiplier on timestep.
    double change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;
    ion.steps = 0;
    ion.V = 0;

    if(log)
    {
        debug_file << "\n\nStarting Ion Trajectory Output\n";
        ion.write_info();
        traj_file << "x\ty\tz\tpx\tpy\tpz\tt\tn\tT\tV\tE\tnear\tdt\tdr_max\tdV\n";
    }
    r.set(ion.r);

start:

    //Find nearby lattice atoms
    ion.fill_nearest(lattice, 4, settings.NPART);
    //Increment counter for how many steps we have taken
    ion.steps++;
    
    //Verify time step is in range.
    dt = std::min(std::max(dt, settings.DELLOW), settings.DELT0);

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, T);

    if(froze || buried || stuck || left || off_edge)
    {
        if(log) debug_file << "Exited" << std::endl;
        goto end;
    }

    ion.V = 0;
    ion.V_t = 0;

    //Find the forces at the current location
    run_hameq(ion, lattice, dt, false);

    dpx = ion.dp_dt[0];
    dpy = ion.dp_dt[1];
    dpz = ion.dp_dt[2];

    if(log)
    {
        debug_file << "\nBetween hameq, dt: " << dt << std::endl;
        ion.write_info();
    }

    //Find the forces at the next location,
    run_hameq(ion, lattice, dt, true);
    V = ion.V;

    dpx = dpx - ion.dp_dt_t[0];
    dpy = dpy - ion.dp_dt_t[1];
    dpz = dpz - ion.dp_dt_t[2];

    dxp = 0.25 * dt * dt * dpx / mass;
    dyp = 0.25 * dt * dt * dpy / mass;
    dzp = 0.25 * dt * dt * dpz / mass;

    dr_max = std::max(fabs(dxp), std::max(fabs(dyp), fabs(dzp)));

    if(log)
    {
        debug_file << "\nAfter hameq " << dt << std::endl;
        debug_file << "V: " << V << std::endl;
        debug_file << "T: " << T << std::endl;
        debug_file << "dR: " << dr_max << std::endl;
        ion.write_info();
    }

    //Now we do some checks to see if the timestep needs to be adjusted
    if(dr_max != 0)
    {
        //check if energy has changed too much.
        change = pow(settings.ABSERR/dr_max, settings.DEMAX);
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
            if(log)
            {
                debug_file << "change_down: " << change << std::endl;
            }
            ion.last_index = -1;
            goto start;
        }
    }
    else
    {
        change = 2;
    }

    //Apply changes, this updates energy and lattice loctations.
    apply_hameq(ion, lattice, dt);
    if(log)
    {
        debug_file << "Applied hameq, dt: " << dt << std::endl;
        debug_file << "change: " << change << std::endl;
        ion.write_info();
    }

    //check if we have gone too far, and need tp re-calculate nearby atoms
    dr = r - ion.r;
    diff = sqr(dr.v);
    if(diff > r_reset)
    {
        //Update ion location for last check
        r.set(ion.r);
        //Reset the index considered as "last found",
        //This forces it to re-calculate nearest atoms
        ion.last_index = -1;
    }

    //Update some parameters for saving.
    ion.r_0[2] = fmin(ion.r[2], ion.r_0[2]);
    ion.time += dt;

    //Update kinetic energy.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;
    V = ion.V;

    //Log things if needed
    if(log)
    {
        V = (ion.V + ion.V_t) / 2;
        dV = ion.V_t - ion.V;
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",
                ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],
                ion.time,ion.steps,T,V,(T+V),ion.near,dt,dr_max,dV);
        traj_file << buffer;
    }

    //Apply differences to timestep;
    dt *= change;

    //Back up we go.
    goto start;

end:
    double E = T;
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
