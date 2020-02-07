#include "hameq.h"

#include "safio.h"      //For settings
#include "traj.h"       //For nearest check method
#include "potentials.h" //For forces/potentials
#include "vec_math.h"   //General maths help
#include "safari.h"     //for exit_fail
#include <cmath>        //sqrt

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
    Atom *atom = s.atom;
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
 * This updates the current location/momentum of
 * the particle, to the corrected values based on the
 * forces at the current and predicted locations.
 * 
 * Positions are corrected by the error in the
 * values from the two forces, but momenta are
 * updated with just the average of the two forces.
 * 
 * This also then updates all of the subsites for this site,
 * note that this function does not get optimized as well as
 * the one that only does a single site.
 * 
 * @param s - the particle to update
 * @param 
 * @param dt - time step for this update.
 */
void update_sites(Site &s, int last_step, double dt)
{
    //We have already been updated!
    if (s.last_step != last_step)
    {
        //Flag us as updated, so we don't get this done again.
        s.last_step = last_step;
        //Site near us.
        update_site(s, dt);
    }
    //Update each nearby site as well
    for (int i = 0; i < s.near; i++)
    {
        Site &s2 = *s.near_sites[i];
        update_sites(s2, last_step, dt);
    }
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
    Atom *atom = s.atom;
    double mass = atom->mass;

    //v = p/m
    //a = F/m, F = dp/dt
    //r_t = r + vdt + 0.5adt^2 = r + dt(p + 0.5*F*dt) / m
    s.r_t[0] = s.r[0] + dt * (s.p[0] + 0.5 * s.dp_dt[0] * dt) / mass;
    s.r_t[1] = s.r[1] + dt * (s.p[1] + 0.5 * s.dp_dt[1] * dt) / mass;
    s.r_t[2] = s.r[2] + dt * (s.p[2] + 0.5 * s.dp_dt[2] * dt) / mass;
}

void apply_hameq(Ion &ion, Lattice &lattice, double dt)
{
    if (!settings.useEinsteinSprings)
    {
        //If we are using better correlations,
        //then we need to update more than just
        //the immediate neighbours, so we need to
        //use the slower version of this.
        update_sites(ion, ion.last_step, dt);
    }
    else
    {
        update_site(ion, dt);
        for (int i = 0; i < ion.near; i++)
        {
            Site &s = *ion.near_sites[i];
            update_site(s, dt);
        }
    }
}

double compute_error(Site &site, double dt)
{
    //Error in the associated coordinate.
    double dxp, dyp, dzp;
    //difference in forces between here and destination
    double dpx, dpy, dpz;
    // //Find difference in the forces.
    dpx = site.dp_dt[0] - site.dp_dt_t[0];
    dpy = site.dp_dt[1] - site.dp_dt_t[1];
    dpz = site.dp_dt[2] - site.dp_dt_t[2];

    //Error in positions between the two forces.
    dxp = 0.25 * dt * dt * dpx / site.atom->mass;
    dyp = 0.25 * dt * dt * dpy / site.atom->mass;
    dzp = 0.25 * dt * dt * dpz / site.atom->mass;

    return std::max(fabs(dxp), std::max(fabs(dyp), fabs(dzp)));
}

void apply_ion_lattice(Ion &ion, Site &s, double *F_at, double *r_i, double ax, double ay, double az, double dt, bool predicted, double *F_ion)
{
    double dx = 0, dy = 0, dz = 0;

    //Distances from site to atom
    dx = ax - r_i[0];
    dy = ay - r_i[1];
    dz = az - r_i[2];

    //Lattice/Ion distance
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    if (r < settings.R_MAX && r > 0)
    {
        ion.r_min = std::min(r, ion.r_min);

        //Magnitude of force for this location.
        double dV_dr = dVr_dr(r, s.atom->index);

        //Potential for this location.
        if (predicted)
        {
            ion.V += Vr_r(r, s.atom->index);
            // debug_file << ion.V << " " << r << " " << ax << " " << ay << " " << az << std::endl;
        }
        //Scaled by 1/r for converting to cartesian
        dV_dr /= r;

        //Convert from magnitude to components of vector
        //Note, fx = -dV_dr * 2dx, however,
        //we set the force to half of this.
        F_at[0] = -dV_dr * dx;
        F_at[1] = -dV_dr * dy;
        F_at[2] = -dV_dr * dz;

        //Add force components to net force.
        F_ion[0] += F_at[0];
        F_ion[1] += F_at[1];
        F_ion[2] += F_at[2];
    }
    //No force if ion is on an atom.
    if (r == 0)
    {
        debug_file << "Ion intersected with atom?: "
                   << ion.index
                   << " " << s.index
                   << " " << ion.steps << std::endl;
        if (num_intersect++ > num_intersect_max)
            exit_fail("Too many ion-site intersections");
    }
}

void run_hameq(Ion &ion, Lattice &lattice, double dt, bool predicted, double *dr_max)
{
    //Some useful variables.
    double dx = 0, dy = 0, dz = 0,
           rr = 0, r = 0;

    //Lattice atom coordintates.
    double ax = 0, ay = 0, az = 0;

    bool recoil = settings.RECOIL;
    bool springs = settings.CORR;
    double atomk = settings.ATOMK;

    bool useEinsteinSprings = settings.useEinsteinSprings;
    bool useLennardJones = settings.useLennardJones;
    // bool useAtomSpings = settings.useAtomSpings;//We default to this if not others.

    //Coordinate of the Ion
    double *r_i;
    r_i = ion.r;

    //Relevant force for this run
    double *F;
    F = ion.dp_dt;

    //if not 0, we are computing for the predicted location.
    if (predicted)
    {
        //Use the one for next time step
        F = ion.dp_dt_t;
        //Use predictions
        r_i = ion.r_t;
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
    if (ion.near)
    {
        double F_ion[3];
        F_ion[0] = 0;
        F_ion[1] = 0;
        F_ion[2] = 0;

        //Check this counter, used for statistics later.
        if (ion.near > ion.max_n)
            ion.max_n = ion.near;

        //Force on the target.
        double *F_at;

        for (int i = 0; i < ion.near; i++)
        {
            //Site near us.
            Site &s = *ion.near_sites[i];

            if (predicted)
            {
                //Get predicted force array.
                F_at = s.dp_dt_t;
                //Get predicted loctation array.
                ax = s.r_t[0];
                ay = s.r_t[1];
                az = s.r_t[2];
            }
            else
            {
                //Get current force array.
                F_at = s.dp_dt;
                //Get current location array.
                ax = s.r[0];
                ay = s.r[1];
                az = s.r[2];
            }

            //Apply site forces.
            //Note that this check is here, as when considering
            //Site neighbours, that inner loop will also apply ion
            //forces if needed.
            if (s.last_step != ion.last_step)
            {
                s.last_step = ion.last_step;
                F_at[0] = 0;
                F_at[1] = 0;
                F_at[2] = 0;
                apply_ion_lattice(ion, s, F_at, r_i, ax, ay, az, dt, predicted, F_ion);
            }

            //Apply spring forces.
            if (springs)
            {
                //Use einstein springs.
                if (useEinsteinSprings)
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
                    if (settings.max_spring_V > V_s)
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

                    //Use atomk as k for springs between nearest sites, or lennard jones potentials
                    for (int j = 0; j < s.near; j++)
                    {
                        Site &s2 = *s.near_sites[j];

                        if (predicted)
                        {
                            //Get predicted force array.
                            F_at2 = s2.dp_dt_t;
                            //Get predicted loctations array.
                            r2 = s2.r_t;
                        }
                        else
                        {
                            F_at2 = s2.dp_dt;
                            r2 = s2.r;
                        }
                        //Reset the sites forces
                        if (s2.last_step != ion.last_step)
                        {
                            s2.last_step = ion.last_step;
                            F_at2[0] = 0;
                            F_at2[1] = 0;
                            F_at2[2] = 0;
                            apply_ion_lattice(ion, s2, F_at2, r_i, r2[0], r2[1], r2[2], dt, predicted, F_ion);
                        }
                        //This will be the magnitude of the force.
                        double dV_dr = 0;

                        //Compute lattice site separation here.
                        dx = ax - r2[0];
                        dy = ay - r2[1];
                        dz = az - r2[2];
                        rr = dx * dx + dy * dy + dz * dz;
                        r = sqrt(rr);

                        if (rr < 0.1)
                        {
                            // debug_file << "Somehow lattice site on other?: "
                            //            << s.index << " "
                            //            << s2.index << " "
                            //            << rr << " "
                            //            << s.rr_min_find << " "
                            //            << hameq_tick <<std::endl;
                            continue;
                        }

                        if (useLennardJones)
                        {
                            //In here we use Lennard Jones forces
                            dV_dr = L_J_dV_dr(r, s2.atom->index, s.atom->index);
                            if (fabs(dV_dr) > 1000)
                            {
                                debug_file << "Somehow large lattice force: "
                                           << dV_dr << ", sep: " << r << " "
                                           << ion.last_step << std::endl;
                                dV_dr = 0;
                            }
                            // debug_file << dV_dr << " " << r << " "
                            //            << sqr(F_at2) << " "
                            //            << sqr(F_at) << " "
                            //            << hameq_tick <<std::endl;
                        }
                        else
                        {
                            //In here we use lattice springs.

                            //Compute lattice site separation at rest
                            dx = s.r_0[0] - s2.r_0[0];
                            dy = s.r_0[1] - s2.r_0[1];
                            dz = s.r_0[2] - s2.r_0[2];
                            rr = dx * dx + dy * dy + dz * dz;
                            double l_eq = sqrt(rr);

                            double dl = r - l_eq;

                            //Check for spring breaking if dl > 0
                            if (dl > 0)
                            {
                                //V = 0.5*k*x^2
                                double V_s = 0.5 * atomk * dl * dl;

                                //Scale by 1/r^2, this accounts for
                                //differences in bond breaking by distance
                                V_s /= l_eq * l_eq;

                                //Break the spring if too far.
                                if (V_s > settings.max_spring_V)
                                    continue;
                            }
                            //F = -kx
                            dV_dr = -atomk * dl;
                        }

                        //These are scaled by 1/r for conversion to
                        //the same coordinate system as the fmag was caluclated for
                        double dx_hat = (r2[0] - ax) / r;
                        double dy_hat = (r2[1] - ay) / r;
                        double dz_hat = (r2[2] - az) / r;
                        //Note, fx = -dV_dr * 2dx, however,
                        //we set the force to half of this.
                        //Apply to us
                        F_at[0] += dx_hat * dV_dr;
                        F_at[1] += dy_hat * dV_dr;
                        F_at[2] += dz_hat * dV_dr;
                        //Apply to other
                        F_at2[0] -= dx_hat * dV_dr;
                        F_at2[1] -= dy_hat * dV_dr;
                        F_at2[2] -= dz_hat * dV_dr;

                        if (recoil)
                        {
                            if (!predicted)
                                //set the predictions
                                predict_site_location(s2, dt);
                            else
                                //set account for site errors
                                *dr_max = std::max(compute_error(s2, dt), *dr_max);
                        }
                    }
                }
            }

            if (recoil)
            {
                if (!predicted)
                    //set the predictions
                    predict_site_location(s, dt);
                else
                    //set account for site errors
                    *dr_max = std::max(compute_error(s, dt), *dr_max);
            }
        }
        //We have now exited the lattice considerations.

        //Superposition of all of the atom-ion interaction
        F[0] = -F_ion[0];
        F[1] = -F_ion[1];
        F[2] = -F_ion[2];

        // debug_file << *dr_max << " "
        //             << hameq_tick <<std::endl;
    }

    if (settings.use_image)
    {
        F[2] -= dVi_dz(r_i[2], ion.q);
        if (predicted)
            ion.V += ion.q * Vi_z(r_i[2], ion.q);
    }

    if (settings.F_a > 0)
    {
        ion.V += apply_friction(lattice, ion, F, dt, predicted);
    }

    if (!predicted)
    {
        //set the predictions
        predict_site_location(ion, dt);
    }
    else
    {
        //Compute error
        *dr_max = std::max(compute_error(ion, dt), *dr_max);
    }
}
