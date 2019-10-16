#include "safio.h"
#include "hameq.h"
#include "potentials.h"
#include "vec_math.h"
#include <iostream>
#include <cmath>


void update_site(Site &s, double dt)
{
    //Site near us.
    Atom atom = s.atom;
    double mass = atom.mass;

    Vec3d dp;
    Vec3d dr;

    Vec3d F;
    Vec3d F_t;

    F.set(s.dp_dt);
    F_t.set(s.dp_dt_t);

    //Average force * dt
    dp = (F + F_t) * 0.5 * dt;

    //Error in predicted and corrected positions
    dr = (F_t - F) * 0.25 * dt * dt / mass;

    //New Location, adjusted by dr
    s.r[0] = s.r_t[0] - dr[0];
    s.r[1] = s.r_t[1] - dr[1];
    s.r[2] = s.r_t[2] - dr[2];

    //New Momentum, not corrected, just using averages
    s.p[0] += dp[0];
    s.p[1] += dp[1];
    s.p[2] += dp[2];
}

void apply_hameq(Ion &ion, Lattice &lattice, double dt)
{
    update_site(ion, dt);
    if(settings.RECOIL)
    {
        int nearby = ion.near;
        for(int i = 0; i<nearby; i++)
        {
            Site s = ion.near_sites[i];
            Vec3d r;
            Vec3d p;
            r.set(s.r);
            update_site(s, dt);
            p.set(s.p);
            r = r - s.r;
            ion.max_site_displacement = std::max(ion.max_site_displacement, r.norm());
            ion.max_site_momentum = std::max(ion.max_site_momentum, p.norm());
            // std::cout << "Position UD: "<< s.r[0]<<" "<< s.r[1]<<" "<< s.r[2] << std::endl;
            // std::cout << "Momentum UD: "<< s.p[0]<<" "<< s.p[1]<<" "<< s.p[2] << std::endl;
            // std::cout << "Force UD: "<< s.dp_dt[0]<<" "<< s.dp_dt[1]<<" "<< s.dp_dt[2] << std::endl;
        }
    }
}

//Sets the next-time-step values to where the particle would be, given
//the current forces and momenta.
void predict_site_location(Site &s, double dt)
{
    //Site near us.
    Atom atom = s.atom;
    double mass = atom.mass;

    //v = p/m
    //a = F/m, F = dp/dt
    //r_t = r + vdt + 0.5adt^2 = r + dt(p + 0.5*F*dt) / m 
    s.r_t[0] = s.r[0] + dt * (s.p[0] + 0.5*s.dp_dt[0]*dt) / mass;
    s.r_t[1] = s.r[1] + dt * (s.p[1] + 0.5*s.dp_dt[1]*dt) / mass;
    s.r_t[2] = s.r[2] + dt * (s.p[2] + 0.5*s.dp_dt[2]*dt) / mass;
}

void run_hameq(Ion &ion, Lattice &lattice, double dt, bool predicted)
{
    //Some useful variables.
    double dx, dy, dz, 
           fx, fy, fz, 
           ftx = 0, fty = 0, ftz = 0;

    //Reset V
    ion.V = 0;

    //Force to populate
    double *F;

    //Normally is this.
    F = ion.dp_dt;

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
    if(predicted)
    {
        //Use the one for next time step
        F = ion.dp_dt_t;
        //Use predictions
        x = ion.r_t[0];
        y = ion.r_t[1];
        z = ion.r_t[2];
    }

    //Initialize the force array.
    F[0] = 0;
    F[1] = 0;
    F[2] = 0;

    //This was set earlier when looking for nearby sites.
    if(ion.near)
    {
        //Check this counter, used for statistics later.
        if(ion.near > ion.max_n)
            ion.max_n = ion.near;

        for(int i = 0; i<ion.near; i++)
        {
            //Site near us.
            Site s = ion.near_sites[i];
            
            //Select the force array to populate.
            F_at = s.dp_dt;

            //Site location
            ax = s.r[0];
            ay = s.r[1];
            az = s.r[2];

            if(predicted)
            {
                //Set predicted force instead.
                F_at = s.dp_dt_t;
                //Get predicted loctations instead.
                ax = s.r_t[0];
                ay = s.r_t[1];
                az = s.r_t[2];
            }

            //Distances from site to atom
            dx = ax - x;
            dy = ay - y;
            dz = az - z;

            double r = sqrt(dx*dx + dy*dy + dz*dz);

            //No force if ion is on an atom.
            if(r==0)
            {
                debug_file << "Ion intersected with atom?" <<std::endl;
                continue;
            }

            //Magnitude of force for this location.
            double dV_dr = dVr_dr(r, s.atom.index);

            //Potential for this location.
            if(predicted) ion.V_t += Vr_r(r, s.atom.index);
            else ion.V += Vr_r(r, s.atom.index);

            //Scaled by 1/r for converting to cartesian
            dV_dr /= r;

            //Convert from magnitude to components of vector
            fx = -2 * dV_dr * dx;
            fy = -2 * dV_dr * dy;
            fz = -2 * dV_dr * dz;

            //Set the initial momentum changes for the sites here.
            //Later other things might adjust them, but this is where
            //they are initialized.
            F_at[0] = fx;
            F_at[1] = fy;
            F_at[2] = fz;

            //Add force components to net force.
            ftx += fx;
            fty += fy;
            ftz += fz;

            if(settings.RECOIL && !predicted)
            {
                //set the predictions
                predict_site_location(s, dt);
            }
        }

        //Superposition of all of the atom-ion interaction
        F[0] = -ftx;
        F[1] = -fty;
        F[2] = -ftz;

        //TODO add atom-atom interactions here.
    }

    if(!predicted)
    {
        //set the predictions
        predict_site_location(ion, dt);
    }
}
