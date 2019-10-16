#include "safio.h"
#include "hameq.h"
#include "potentials.h"
#include <iostream>


void update_site(site &s)
{
    //New Location
    s[0] = s.r_t[0];
    s[1] = s.r_t[1];
    s[2] = s.r_t[2];

    //New Momentum
    s.p[0] = s.p_t[0];
    s.p[1] = s.p_t[1];
    s.p[2] = s.p_t[2];
}

void apply_hameq(ion &ion, lattice &lattice, double dt)
{
    update_site(ion);
    int nearby = ion.near;
    if(settings.RECOIL)
        for(int i = 0; i<nearby; i++)
        {
            site s = ion.near_sites[i];
            update_site(s);
        }
}

void predict_site(site &s, double dt)
{
    //Site near us.
    atom atom = s.atom;
    double mass = atom.mass;

    //v = p/m
    //a = F/m, F = dp/dt
    //r_t = r + vt + 0.5at^2
    s.r_t[0] = s.r[0] + dt * (s.p[0] + 0.5*s.dp_dt[0]*dt) / mass;
    s.r_t[1] = s.r[1] + dt * (s.p[1] + 0.5*s.dp_dt[1]*dt) / mass;
    s.r_t[2] = s.r[2] + dt * (s.p[2] + 0.5*s.dp_dt[2]*dt) / mass;

    //p_t = p + dt * F
    s.p_t[0] = s.p[0] +  dt * s.dp_dt[0];
    s.p_t[1] = s.p[1] +  dt * s.dp_dt[1];
    s.p_t[2] = s.p[2] +  dt * s.dp_dt[2];
}

void run_hameq(ion &ion, lattice &lattice, double *T, double dt)
{
    //Some useful variables.
    double dx, dy, dz, 
           fx, fy, fz, 
           ftx = 0, fty = 0, ftz = 0;

    //Reset V
    ion.V = 0;

    //Force to populate
    double *force;
    //Normally is this.
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
    
    double mass = ion.atom.mass;
    double psq;

    double px = ion.p[0];
    double py = ion.p[1];
    double pz = ion.p[2];

    //if not 0, we are computing for the predicted location.
    if(dt != 0)
    {
        //Use the one for next time step
        force = ion.dp_dt_t;

        predict_site(ion, dt);

        x = ion.r_t[0];
        y = ion.r_t[1];
        z = ion.r_t[2];

        px = ion.p_t[0];
        py = ion.p_t[1];
        pz = ion.p_t[2];
    }

    psq = px*px + py*py + pz*pz;
    *T = 0.5 * psq / mass;

    //Initialize the force array.
    force[0] = 0;
    force[1] = 0;
    force[2] = 0;

    ion.V = 0;
    // debug_file << "Sites Near?: " << ion.near << std::endl;

    //This was set earlier when looking for nearby sites.
    if(ion.near)
    {
        //Check this counter, used for statistics later.
        if(ion.near > ion.max_n)
            ion.max_n = ion.near;

        for(int i = 0; i<ion.near; i++)
        {
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
                predict_site(s, dt);
                //Get predicted loctations instead.
                ax = s.r_t[0];
                ay = s.r_t[1];
                az = s.r_t[2];
                F_at = s.dp_dt_t;
            }

            double r = s.distance(ion, dt!=0);

            //No force if ion is on an atom.
            if(r==0)
            {
                debug_file << "Ion intersected with atom?" <<std::endl;
                continue;
            }

            //Magnitude of force for this location.
            double dV_dr = dVr_dr(r, s.atom.index);

            //Potential for this location.
            ion.V += Vr_r(r, s.atom.index);

            //Scaled by 1/r for converting to cartesian
            dV_dr /= r;

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

        //TODO add atom-atom interactions here.
    }
}
