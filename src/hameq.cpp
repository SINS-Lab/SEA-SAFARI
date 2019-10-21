#include "safio.h"
#include "hameq.h"
#include "potentials.h"
#include "vec_math.h"
#include <iostream>
#include <cmath>

/**
 * This updates the current location/momentum of
 * the particle, to the corrected values based on the
 * forces at the current and predicted locations.
 * 
 * Positions are corrected by the error in the
 * values from the two forces, but momenta are
 * updated with just the average of the two forces.
 * 
 * @param s - the particle to update
 * @param dt - time step for this update.
 */ 
void update_site(Site &s, double dt)
{
    //Site near us.
    Atom &atom = s.atom;
    double mass = atom.mass;

    //We use averaged forces everywhere, so just cut
    //this in half here for 6 less * operations
    dt *= 0.5;

    //New Location, adjusted by error in corrected positions
    //corrected positon is 0.25*dt^2*(F_t - F)/mass
    s.r[0] = s.r_t[0] - dt * dt * (s.dp_dt_t[0] - s.dp_dt[0]) / mass;
    s.r[1] = s.r_t[1] - dt * dt * (s.dp_dt_t[1] - s.dp_dt[1]) / mass;
    s.r[2] = s.r_t[2] - dt * dt * (s.dp_dt_t[2] - s.dp_dt[2]) / mass;

    //New Momentum, not corrected, just using averages
    s.p[0] += dt * (s.dp_dt_t[0] + s.dp_dt[0]);
    s.p[1] += dt * (s.dp_dt_t[1] + s.dp_dt[1]);
    s.p[2] += dt * (s.dp_dt_t[2] + s.dp_dt[2]);
}

/**
 * This predicts the location of the site, assuming the
 * current value of the momentum and force.
 * 
 * @param s - the particle to predict.
 * @param dt - the time step for the prediction
 */ 
void predict_site_location(Site &s, double dt)
{
    //Site near us.
    Atom &atom = s.atom;
    double mass = atom.mass;

    //v = p/m
    //a = F/m, F = dp/dt
    //r_t = r + vdt + 0.5adt^2 = r + dt(p + 0.5*F*dt) / m 
    s.r_t[0] = s.r[0] + dt * (s.p[0] + 0.5*s.dp_dt[0]*dt) / mass;
    s.r_t[1] = s.r[1] + dt * (s.p[1] + 0.5*s.dp_dt[1]*dt) / mass;
    s.r_t[2] = s.r[2] + dt * (s.p[2] + 0.5*s.dp_dt[2]*dt) / mass;
}

void apply_hameq(Ion &ion, Lattice &lattice, double dt)
{
    //Update the ion's location
    update_site(ion, dt);
    if(settings.RECOIL)
    {
        int nearby = ion.near;
        //Update each site.
        for(int i = 0; i<nearby; i++)
        {
            Site &s = *ion.near_sites[i];
            update_site(s, dt);
        }
    }
}

void run_hameq(Ion &ion, Lattice &lattice, double dt, bool predicted)
{
    //Some useful variables.
    double dx = 0, dy = 0, dz = 0, 
           fx = 0, fy = 0, fz = 0, 
           ftx = 0, fty = 0, ftz = 0;

    //Ion coordinates
    double x = 0, y = 0, z = 0;

    //Lattice atom coordintates.
    double ax = 0, ay = 0, az = 0;

    bool recoil = settings.RECOIL;
    bool springs = settings.CORR;
    double atomk = settings.ATOMK; 

    //Initialize ion coordinates
    x = ion[0];
    y = ion[1];
    z = ion[2];

    //Relevant force for this run
    double *F;
    F = ion.dp_dt;

    //Relevant potential for this run
    double *V;
    V = &ion.V;
    
    //if not 0, we are computing for the predicted location.
    if(predicted)
    {
        //Use the one for next time step
        F = ion.dp_dt_t;
        V = &ion.V_t;
        //Use predictions
        x = ion.r_t[0];
        y = ion.r_t[1];
        z = ion.r_t[2];
    }

    //Reset V
    *V = 0;

    //Reset the force array.
    F[0] = 0;
    F[1] = 0;
    F[2] = 0;

    //This was set earlier when looking for nearby sites.
    if(ion.near)
    {
        //Check this counter, used for statistics later.
        if(ion.near > ion.max_n)
            ion.max_n = ion.near;

        //Force on the target.
        double *F_at;

        for(int i = 0; i<ion.near; i++)
        {
            //Site near us.
            Site &s = *ion.near_sites[i];
            
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
                // debug_file << "Ion intersected with atom?" <<std::endl;
                continue;
            }

            ion.r_min = std::min(r, ion.r_min);

            //Magnitude of force for this location.
            double dV_dr = dVr_dr(r, s.atom.index);

            //Potential for this location.
            *V += Vr_r(r, s.atom.index);

            //Scaled by 1/r for converting to cartesian
            dV_dr /= r;

            //Convert from magnitude to components of vector
            //Note, fx = -dV_dr * 2dx, however,
            //we set the force to half of this.
            fx = -dV_dr * dx;
            fy = -dV_dr * dy;
            fz = -dV_dr * dz;

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

            //Apply spring forces.
            if(springs)
            {
                //Use einstein springs.
                if(atomk <= 1e-30)
                {   
                    //F = -kx
                    F_at[0] -= s.atom.spring[0] * (s.r - s.r_0);
                    F_at[1] -= s.atom.spring[1] * (s.r - s.r_0);
                    F_at[2] -= s.atom.spring[2] * (s.r - s.r_0);
                }
                else
                {
                    //Need to loop over all and do some calculation.
                }
                
            }


            if(recoil && !predicted)
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
