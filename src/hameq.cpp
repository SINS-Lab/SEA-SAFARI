#include "safio.h"
#include "hameq.h"
#include "potentials.h"
#include <iostream>

void apply_hameq(ion &ion, lattice &lattice, double dt)
{

    //velocities
    double dr_dt[3];
    //Accelerations
    double d2r_dt2[3];

    double mass = ion.atom.mass;

    //set the force to average of the here and next
    ion.dp_dt[0] = (ion.dp_dt[0] + ion.dp_dt_t[0]) * 0.5;
    ion.dp_dt[1] = (ion.dp_dt[1] + ion.dp_dt_t[1]) * 0.5;
    ion.dp_dt[2] = (ion.dp_dt[2] + ion.dp_dt_t[2]) * 0.5;

    //v = p/m
    dr_dt[0] = ion.p[0] / mass;
    dr_dt[1] = ion.p[1] / mass;
    dr_dt[2] = ion.p[2] / mass;

    //a = F/m, F = dp/dt
    d2r_dt2[0] = ion.dp_dt[0] / mass;
    d2r_dt2[1] = ion.dp_dt[1] / mass;
    d2r_dt2[2] = ion.dp_dt[2] / mass;

    //New Ion Location
    ion[0] += dt * (dr_dt[0] + 0.5*d2r_dt2[0]*dt);
    ion[1] += dt * (dr_dt[1] + 0.5*d2r_dt2[1]*dt);
    ion[2] += dt * (dr_dt[2] + 0.5*d2r_dt2[2]*dt);

    //New Ion Momentum
    ion.p[0] += ion.dp_dt[0] * dt;
    ion.p[1] += ion.dp_dt[1] * dt;
    ion.p[2] += ion.dp_dt[2] * dt;

    int nearby = ion.near;
    if(settings.RECOIL)
        for(int i = 0; i<nearby; i++)
        {
            site s = ion.near_sites[i];
            atom a = s.atom;

            mass = a.mass;

            //set the force to average of the here and next
            s.dp_dt[0] = (s.dp_dt[0] + s.dp_dt_t[0]) * 0.5;
            s.dp_dt[1] = (s.dp_dt[1] + s.dp_dt_t[1]) * 0.5;
            s.dp_dt[2] = (s.dp_dt[2] + s.dp_dt_t[2]) * 0.5;

            //v = p/m
            dr_dt[0] = s.p[0] / mass;
            dr_dt[1] = s.p[1] / mass;
            dr_dt[2] = s.p[2] / mass;

            //a = F/m, F = dp/dt
            d2r_dt2[0] = s.dp_dt[0] / mass;
            d2r_dt2[1] = s.dp_dt[1] / mass;
            d2r_dt2[2] = s.dp_dt[2] / mass;

            //New Ion Location
            s[0] += dt * (dr_dt[0] + 0.5*d2r_dt2[0]*dt);
            s[1] += dt * (dr_dt[1] + 0.5*d2r_dt2[1]*dt);
            s[2] += dt * (dr_dt[2] + 0.5*d2r_dt2[2]*dt);

            //New Ion Momentum
            s.p[0] += s.dp_dt[0] * dt;
            s.p[1] += s.dp_dt[1] * dt;
            s.p[2] += s.dp_dt[2] * dt;
        }
}

void run_hameq(ion &ion, lattice &lattice, double dt)
{
    //Some useful variables.
    double dx, dy, dz, 
           fx, fy, fz, 
           ftx = 0, fty = 0, ftz = 0, 
           rmin = 1e20;

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

    ion.v_total = 0;

    //This was set earlier when looking for nearby sites.
    if(ion.near)
    {
        //Check this counter, used for statistics later.
        if(ion.near > ion.max_n)
            ion.max_n = ion.near;

        if(settings.RECOIL)
            for(int i = 0; i<ion.near; i++)
            {
                //Site near us.
                site s = ion.near_sites[i];
                atom atom = s.atom;
                double mass = atom.mass;
                //Ion velocities
                double dr_dt[3];
                //Ion Accelerations
                double d2r_dt2[3];

                //v = p/m
                dr_dt[0] = s.p[0] / mass;
                dr_dt[1] = s.p[1] / mass;
                dr_dt[2] = s.p[2] / mass;
                //a = F/m, F = dp/dt
                d2r_dt2[0] = s.dp_dt[0] / mass;
                d2r_dt2[1] = s.dp_dt[1] / mass;
                d2r_dt2[2] = s.dp_dt[2] / mass;

                s.r_t[0] = s.r[0] + dt * (dr_dt[0] + 0.5*d2r_dt2[0]*dt);
                s.r_t[1] = s.r[1] + dt * (dr_dt[1] + 0.5*d2r_dt2[1]*dt);
                s.r_t[2] = s.r[2] + dt * (dr_dt[2] + 0.5*d2r_dt2[2]*dt);
            }

        //Update distances to the ion.
        ion.check_distances(x,y,z, dt!=0 && settings.RECOIL);

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
            double dV_dr = *(forces + i);

            //Scaled by 1/r for converting to cartesian
            dV_dr /= r;

            //Site near us.
            site s = ion.near_sites[i];
            //Select the force array to populate.
            F_at = s.dp_dt;

            //Site location
            ax = s.r[0];
            ay = s.r[1];
            az = s.r[2];

            if(dt!=0 && settings.RECOIL)
            {
                //Get predicted loctations instead.
                ax = s.r_t[0];
                ay = s.r_t[1];
                az = s.r_t[2];
                F_at = s.dp_dt_t;
            }

            //Distances from site to atom
            dx = ax - x;
            dy = ay - y;
            dz = az - z;

            //Convert from magnitude to components of vector
            fx = -2 * dV_dr * dx;
            fy = -2 * dV_dr * dy;
            fz = -2 * dV_dr * dz;

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

        //Superposition of all of the atom-ion interaction
        force[0] = -ftx;
        force[1] = -fty;
        force[2] = -ftz;

        //Add the ion-atom potential to the v_total.
       // if(dt != 0) 
        ion.v_total += Vr_r(ion.near_dists, ion.near_atoms, ion.near);

        //TODO add atom-atom interactions here.
    }
}
