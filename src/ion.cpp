#include "ion.h"
#include <math.h>
#include "potentials.h"
#include "safio.h"
#include "hameq.h"
#include "lattice.h"
#include "vec_math.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <cstdio>

#define M_PI           3.14159265358979323846  /* pi */

ion::~ion()
{
}

ion::ion()
{
}

void ion::check_distances(double x, double y, double z)
{
    r_min = 1e3;
    for(int i = 0; i<near; i++)
    {
        site s = near_sites[i];

        double dx = x - s[0];
        double dy = y - s[1];
        double dz = z - s[2];

        double r = sqrt(dx*dx + dy*dy + dz*dz);
        near_dists[i] = r;
        r_min = std::min(r_min, r);
    }
}

int ion::fill_nearest(lattice &lattice)
{
    ion ion = *this;

    int cell_x = ion[0];
    int cell_y = ion[1];
    int cell_z = ion[2];

    int pos_hash = to_hash(cell_x, cell_y, cell_z);
    if(pos_hash == last_index)
        return near;
    last_index = pos_hash;

    int r = 3;
    int num;
    near = 0;
    for(int i = -r; i<=r; i++)
    {
        for(int j = -r; j<=r; j++)
        {
            for(int k = -r; k<=r; k++)
            {
                pos_hash = to_hash(cell_x + i, cell_y + j,cell_z + k);
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
                    near_atoms[near] = lattice.atoms[s.index].index;
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
    double dt = 0.001;
    //Kinetic Energy of the ion
    double E;
    //Magnitude squared of the ion momentum
    double psq;
    //Ion mass
    double mass = ion.atom.mass;

    double dx, dy, dz;

    //Momenta changes before time t
    double dpz_dt, dpy_dt, dpx_dt;
    //Momenta of particle
    double px, py, pz;

    //Predicted position of ion
    double xp, yp, zp;
    //Predicted Momenta changes
    double dpzp_dt, dpyp_dt, dpxp_dt;

    double pp = 0, pzz = 0, theta, phi;

    //Control conditions
    bool buried = false;
    bool froze = false;
    bool stuck = false;
    bool off_edge = false;
    bool left = false;

    //Parameters for checking
    //how much things are changing by.
    double err_max, ex, ey, ez, change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = ion.p[0]*ion.p[0] + ion.p[1]*ion.p[1] + ion.p[2]*ion.p[2];
    E = psq * 0.5 / mass;
    ion.steps = 0;

    bool dynamic_dt = false;

start:

    //Find nearby lattice atoms
    ion.fill_nearest(lattice);
    //Increment counter for how many steps we have taken
    ion.steps++;

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, &dt, E);
    if(froze || buried || stuck || left || off_edge)
    {
        if(log) std::cout << "Exited" << std::endl;
        goto end;
    }

    //reset some values.
    if(!dynamic_dt) dt = 0.001;

    //Initial values for these.
    px = ion.p[0];
    py = ion.p[1];
    pz = ion.p[2];

    //Find the time derivatives.
    run_hameq(ion, lattice, dt, false, &xp, &yp, &zp);

    //These are needed for error calculations.
    //after the second pass of run_hameq
    dpx_dt = ion.dp_dt[0];
    dpy_dt = ion.dp_dt[1];
    dpz_dt = ion.dp_dt[2];

    //Check time derivatives at predicted position,
    //This is to check if we need to reduce time steps.
    run_hameq(ion, lattice, dt, true, &xp, &yp, &zp);

    dpxp_dt = ion.dp_dt[0];
    dpyp_dt = ion.dp_dt[1];
    dpzp_dt = ion.dp_dt[2];

    //Compute differences between initial and predicted
    //momenta from the two runs of run_hameq
    //uses average forces
    px += 0.5 * dt * (dpxp_dt + dpx_dt);
    px += 0.5 * dt * (dpyp_dt + dpy_dt);
    px += 0.5 * dt * (dpzp_dt + dpz_dt);

    //Compute the displacements resulting from the forces.
    ex = 0.25 * dt * dt * (dpxp_dt - dpx_dt) / mass;
    ey = 0.25 * dt * dt * (dpyp_dt - dpy_dt) / mass;
    ez = 0.25 * dt * dt * (dpzp_dt - dpz_dt) / mass;

    //Compute maximum change.
    err_max = std::max(fabs(ex), std::max(fabs(ez), fabs(ey)));

    //TODO also include lattice recoil for err_max here.

    //Compute amount for time to step by.
    if(dynamic_dt && err_max!=0)
    {
        change = pow(settings.ABSERR/err_max, settings.DEMAX);
        if(change >= 2) change = 2;
        if(change > 1 && change < 2) change = 1;
        //Large change in energy, try to reduce timestep.
        if(change < .2 && dt > settings.DELLOW)
        {
            dt *= change;
            if(dt < settings.DELLOW) dt = settings.DELLOW;
            goto start;
        }
    }
    else
    {
        change = 2;
    }

    //New Ion Location
    xp += ex;
    yp += ey;
    zp += ez;

    //Change in position
    dx = xp - ion[0];
    dy = yp - ion[1];
    dz = zp - ion[2];

    ion[0] = xp;
    ion[1] = yp;
    ion[2] = zp;

    //Update ion momenta
    ion.p[0] = mass * dx / dt;
    ion.p[1] = mass * dy / dt;
    ion.p[2] = mass * dz / dt;

    //Apply changes, this updates energy and lattice loctations.
    apply_hameq(ion, lattice, dt);

    //Update some parameters for saving.
    ion.r_0[2] = fmin(ion.r[2], ion.r_0[2]);
    ion.time += dt;

    //Update kinetic energy.
    psq = px*px + py*py + pz*pz;
    E = psq * 0.5 / mass;

    //Log things if needed
    if(log)
    {
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n",
                ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],
                ion.time,ion.steps,E,ion.v_total,(E+ion.v_total), ion.r_min,dt);
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
        double  py = ion.p[1];
        double  pz = ion.p[2];
        psq = px*px + py*py + pz*pz;
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
//            sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",ion.r_0[0],ion.r_0[1],ion.r_0[2],E,theta,phi,ion.time,ion.steps,ion.max_n);
//            out_file << buffer;

        }
    }
    sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",
            ion.r_0[0],ion.r_0[1],ion.r_0[2],
            E,theta,phi,ion.time,ion.steps,ion.max_n);
    out_file << buffer;
    return;
}
