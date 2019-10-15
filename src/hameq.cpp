#include "safio.h"
#include "hameq.h"
#include "potentials.h"
#include <iostream>

void apply_hameq(ion &ion, lattice &lattice, double dt)
{
    int nearby = ion.near;
    if(settings.RECOIL)
        for(int i = 0; i<nearby; i++)
        {
            site s = ion.near_sites[i];
            atom a = lattice.atoms[ion.near_atoms[i]];
            s.p[0] += s.dp_dt[0] * dt;
            s.p[1] += s.dp_dt[1] * dt;
            s.p[2] += s.dp_dt[2] * dt;

            s[0] += s.p[0] * dt / a.mass;
            s[1] += s.p[1] * dt / a.mass;
            s[2] += s.p[2] * dt / a.mass;

            s.dp_dt[0] = s.dp_dt[1] = s.dp_dt[2] = 0;
        }

    //Update distances to the ion.
    ion.check_distances(ion[0], ion[1], ion[2]);
    //Add the ion-atom potential to the v_total.
    ion.v_total += Vr_r(ion.near_dists, ion.near_atoms, nearby);
}

void run_hameq(ion &ion, lattice &lattice, double dt)
{
    //Some useful variables.
    double dx, dy, dz, fx, fy, fz, ftx = 0, fty = 0, ftz = 0, rmin = 1e20;

    //Reset v_total
    ion.v_total = 0;

    //Force to populate
    double *force;

    force = ion.dp_dt;

    //Force on the target.
    double *F_at;

    //Ion coordinates
    double x, y, z;

    //Lattice atom coordintates.
    double ax, ay, az;

    //Initialize ion coordinates
    x = ion[0];
    y = ion[1];
    z = ion[2];

    //if not 0, we are computing for the predicted location.
    if(dt != 0)
    {
        force = ion.dp_dt_t;

        double mass = ion.atom.mass;
        //Ion velocities
        double dx_dt, dy_dt, dz_dt;
        //Ion Accelerations
        double dvx_dt, dvy_dt, dvz_dt;

        //v = p/m
        dx_dt = ion.p[0] / mass;
        dy_dt = ion.p[1] / mass;
        dz_dt = ion.p[2] / mass;
        //a = F/m, F = dp/dt
        dvx_dt = ion.dp_dt[0] / mass;
        dvy_dt = ion.dp_dt[1] / mass;
        dvz_dt = ion.dp_dt[2] / mass;

        x += dt * (dx_dt + 0.5*dvx_dt*dt);
        y += dt * (dy_dt + 0.5*dvy_dt*dt);
        z += dt * (dz_dt + 0.5*dvz_dt*dt);
    }

    //Initialize the force array.
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;
    //This was set earlier when looking for nearby sites.
    if(ion.near)
    {
        //Check this counter, used for statistics later.
        if(ion.near > ion.max_n)
            ion.max_n = ion.near;

        //Update distances to the ion.
        ion.check_distances(x,y,z);

        double *forces = dVr_dr(ion.near_dists, ion.near_atoms, ion.near);
        for(int i = 0; i<ion.near; i++)
        {
            //Atom-Ion distance.
            double r = ion.near_dists[i];

            //No force if ion is on an atom.
            if(r==0)
            {
                debug_file << "Ion intersected with atom?" <<std::endl;
                continue;
            }
            //Record this for tracking later.
            if(r < rmin)
                rmin = r;

            //Magnitude of force for this location.
            //Scaled by 1/r for use in vector later.
            double dV_dr = *(forces + i);

            dV_dr /= r;

            //Site near us.
            site s = ion.near_sites[i];
            //Select the force array to populate.
            F_at = s.dp_dt;

            //Site location
            ax = s[0];
            ay = s[1];
            az = s[2];

            if(dt!=0)
            {
                F_at = s.dp_dt_t;

                //TODO adjust ax,ay,az if recoil.
            }

            //Distances from site to atom
            dx = ax - x;
            dy = ay - y;
            dz = az - z;

            //Convert from magnitude to components of vector
            fx = dV_dr * dx;
            fy = dV_dr * dy;
            fz = dV_dr * dz;

            //Add force components to net force.
            ftx += fx;
            fty += fy;
            ftz += fz;

            //Set the initial momentum changes for the sites here.
            //Later other things might adjust them, but this is where
            //they are initialized.
            F_at[0] = fx;
            F_at[1] = fy;
            F_at[2] = fz;
        }

        //debug_file << ftx << " " << fty << " " << ftz << std::endl;

        //Superposition of all of the atom-ion interaction
        force[0] -= ftx;
        force[1] -= fty;
        force[2] -= ftz;
        //TODO add atom-atom interactions here.
    }
}
