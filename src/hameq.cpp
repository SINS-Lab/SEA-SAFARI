#include "hameq.h" 

#include "safio.h"       //For settings
#include "traj.h"        //For nearest check method      
#include "potentials.h"  //For forces/potentials
#include "vec_math.h"    //General maths help
#include "safari.h"      //for exit_fail
#include <cmath>         //sqrt

int num_intersect_max = 1000000;
//This is a counter for number of intersections, if this exceeds max, it terminates safari.
int num_intersect = 0;

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
    Atom* atom = s.atom;
    double mass = atom->mass;

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
    Atom* atom = s.atom;
    double mass = atom->mass;

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
           ftx = 0, fty = 0, ftz = 0,
           rr = 0;

    //Lattice atom coordintates.
    double ax = 0, ay = 0, az = 0;

    bool recoil = settings.RECOIL;
    bool springs = settings.CORR;
    double atomk = settings.ATOMK; 

    //Coordinate of the Ion
    double *r;
    r = ion.r;

    //Relevant force for this run
    double *F;
    F = ion.dp_dt;
    
    //if not 0, we are computing for the predicted location.
    if(predicted)
    {
        //Use the one for next time step
        F = ion.dp_dt_t;
        //Use predictions
        r = ion.r_t;
    }

    //Update "last integration step" for ion
    ion.last_step++;

    //Reset V
    ion.V = 0;

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

            //Reset the sites forces
            if(s.last_step!=ion.last_step)
            {
                F_at[0] = 0;
                F_at[1] = 0;
                F_at[2] = 0;
                s.last_step = ion.last_step;
            }

            //Distances from site to atom
            dx = ax - r[0];
            dy = ay - r[1];
            dz = az - r[2];

            //Lattice/Ion distance no longer need ion's r
            //so we can shadow it here.
            double r = sqrt(dx*dx + dy*dy + dz*dz);

            //No force if ion is on an atom.
            if(r==0)
            {
                debug_file << "Ion intersected with atom?: " << ion.index << " " << s.index << std::endl;
                if(num_intersect++>num_intersect_max) exit_fail("Too many ion-site intersections");
                continue;
            }

            ion.r_min = std::min(r, ion.r_min);

            //Magnitude of force for this location.
            double dV_dr = dVr_dr(r, s.atom->index);

            //Potential for this location.
            if(predicted) ion.V += Vr_r(r, s.atom->index);

            //Scaled by 1/r for converting to cartesian
            dV_dr /= r;

            //Convert from magnitude to components of vector
            //Note, fx = -dV_dr * 2dx, however,
            //we set the force to half of this.
            fx = -dV_dr * dx;
            fy = -dV_dr * dy;
            fz = -dV_dr * dz;

            //Apply momentum change to the site
            F_at[0] += fx;
            F_at[1] += fy;
            F_at[2] += fz;

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
                    double dx_sq = (s.r[0] - s.r_0[0]) * (s.r[0] - s.r_0[0]);
                    double dy_sq = (s.r[1] - s.r_0[1]) * (s.r[1] - s.r_0[1]);
                    double dz_sq = (s.r[2] - s.r_0[2]) * (s.r[2] - s.r_0[2]);
                    //V = 0.5kx^2
                    double V_s_x = 0.5 * s.atom->spring[0] * dx_sq;
                    double V_s_y = 0.5 * s.atom->spring[1] * dy_sq;
                    double V_s_z = 0.5 * s.atom->spring[2] * dz_sq;

                    double V_s = V_s_x + V_s_y + V_s_z;
                    //Check if too much spring energy
                    if(settings.max_spring_V > V_s)
                    {
                        //F = -kx
                        F_at[0] -= s.atom->spring[0] * (s.r[0] - s.r_0[0]);
                        F_at[1] -= s.atom->spring[1] * (s.r[1] - s.r_0[1]);
                        F_at[2] -= s.atom->spring[2] * (s.r[2] - s.r_0[2]);
                    }
                    //These springs are only applied if the site is still "near" the original location.
                }
                else
                {
                    //Force on target.
                    double *F_at2;
                    double *r2;

                    //Use atomk as k for springs between nearest sites.
                    for(int j = 0; j<s.near; j++)
                    {
                        Site &s2 = *s.near_sites[j];

                        if(predicted)
                        {
                            //Set predicted force instead.
                            F_at2 = s2.dp_dt_t;
                            //Get predicted loctations instead.
                            r2 = s2.r_t;
                        }
                        else
                        {
                            F_at2 = s2.dp_dt;
                            r2 = s2.r;
                        }
                        //Reset the sites forces
                        if(s2.last_step!=s.last_step)
                        {
                            F_at2[0] = 0;
                            F_at2[1] = 0;
                            F_at2[2] = 0;
                            s2.last_step = s.last_step;
                        }
                        dx = s.r_0[0] - s2.r_0[0];
                        dy = s.r_0[1] - s2.r_0[1];
                        dz = s.r_0[2] - s2.r_0[2];
                        rr = dx*dx + dy*dy + dz*dz;
                        double l_eq = sqrt(rr);
                        dx = ax - r2[0];
                        dy = ay - r2[1];
                        dz = az - r2[2];
                        rr = dx*dx + dy*dy + dz*dz;
                        double ll_now = rr;
                        double dl = sqrt(ll_now) - l_eq;

                        //Check for spring breaking if dl > 0
                        if(dl > 0)
                        {
                            //V = 0.5*k*x^2
                            double V_s = 0.5 * atomk * dl * dl;

                            //Scale by 1/r^2, this accounts for
                            //differences in bond breaking by distance
                            V_s /= l_eq * l_eq;

                            //Break the spring if too far.
                            if(V_s > settings.max_spring_V) continue;
                        }


                        //F = -kx
                        double f_mag = -atomk * dl;
                        //These are scaled by 1/r^2 for conversion to
                        //the same coordinate system as the fmag was caluclated for
                        double dx_hat = (r2[0] - ax) / ll_now;
                        double dy_hat = (r2[1] - ay) / ll_now;
                        double dz_hat = (r2[2] - az) / ll_now;
                        //Apply to us
                        F_at[0] += dx_hat * f_mag;
                        F_at[1] += dy_hat * f_mag;
                        F_at[2] += dz_hat * f_mag;
                        //Apply to other
                        F_at2[0] -= dx_hat * f_mag;
                        F_at2[1] -= dy_hat * f_mag;
                        F_at2[2] -= dz_hat * f_mag;

                    }
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

    if(settings.use_image)
    {
        F[2] -= dVi_dz(r[2], ion.q);
        if(predicted) ion.V += ion.q * Vi_z(r[2], ion.q);
    }

    if(settings.F_a > 0)
    {
        ion.V += apply_friction(lattice, ion, F, dt);
    }

    if(!predicted)
    {
        //set the predictions
        predict_site_location(ion, dt);
    }
}
