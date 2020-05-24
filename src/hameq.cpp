#include "hameq.h"

#include "safio.h"      // for settings
#include "traj.h"       // for nearest check method
#include "potentials.h" // for forces/potentials
#include "vec_math.h"   // General maths help
#include "safari.h"     // for exit_fail
#include <cmath>        // sqrt

void update_site(Site &s, double dt)
{
    Atom *atom = s.atom;
    double mass = atom->mass;

    // We use averaged forces everywhere, so just cut
    // this in half here for 6 less * operations
    dt *= 0.5;

    // New Location, adjusted by error in corrected positions
    // corrected positon is 0.25*dt^2*(F_t - F)/mass
    s.r[0] = s.r_t[0] - dt * dt * (s.dp_dt_t[0] - s.dp_dt[0]) / mass;
    s.r[1] = s.r_t[1] - dt * dt * (s.dp_dt_t[1] - s.dp_dt[1]) / mass;
    s.r[2] = s.r_t[2] - dt * dt * (s.dp_dt_t[2] - s.dp_dt[2]) / mass;

    // New Momentum, not corrected, just using averages
    s.p[0] += dt * (s.dp_dt_t[0] + s.dp_dt[0]);
    s.p[1] += dt * (s.dp_dt_t[1] + s.dp_dt[1]);
    s.p[2] += dt * (s.dp_dt_t[2] + s.dp_dt[2]);

    s.reset_forces();
}

void update_sites(Site &s, int last_update, double dt)
{
    // We have already been updated!
    if (s.last_update != last_update)
    {
        // Flag us as updated, so we don't get this done again.
        s.last_update = last_update;
        // Site near us.
        update_site(s, dt);
        // Update each nearby site as well
        for (int i = 0; i < s.near; i++)
        {
            Site &s2 = *s.near_sites[i];
            update_sites(s2, last_update, dt);
        }
    }
}

void predict_site_location(Site &s, double dt)
{
    Atom *atom = s.atom;
    double mass = atom->mass;

    // v = p/m
    // a = F/m, F = dp/dt
    // r_t = r + vdt + 0.5adt^2 = r + dt(p + 0.5*F*dt) / m
    s.r_t[0] = s.r[0] + dt * (s.p[0] + 0.5 * s.dp_dt[0] * dt) / mass;
    s.r_t[1] = s.r[1] + dt * (s.p[1] + 0.5 * s.dp_dt[1] * dt) / mass;
    s.r_t[2] = s.r[2] + dt * (s.p[2] + 0.5 * s.dp_dt[2] * dt) / mass;
}

bool check_sputter(Ion &ion, Site *s)
{
    // we are already marked!
    if (s->left)
        return false;
    if (settings.cascadeMode)
    {
        bool isHot = (sqr(s->p) * 0.5 / s->atom->mass) > 1;
        double dr = diff_sqr(s->r_0, s->r);
        if (isHot and dr > settings.AX / 2)
        {
            s->left = true;
            return true;
        }
    }
    else if (settings.saveSputter)
    {

        double pz = s->p[2];
        double rz = s->r[2];
        // TODO better conditions for leaving surface
        if (pz > 0 and rz > settings.Z1 / 4)
        {
            s->left = true;
            ion.sputter[ion.sputtered] = s;
            ion.sputtered++;
            return true;
        }
    }
    return false;
}

void apply_hameq(Ion &ion, Lattice *lattice, double dt)
{
    if (!settings.useEinsteinSprings)
    {
        update_site(ion, dt);
        // In this case, we update all of the sites
        if (settings.dynamicNeighbours)
        {
            int num = lattice->active_sites.size();
            for (int i = 0; i < num; i++)
            {
                Site *s = lattice->active_sites[i];
                update_site(*s, dt);
                check_sputter(ion, s);
                s->last_update = ion.steps;
            }
        }
        else
        {
            // Otherwise, just consider ones tracked by the ion
            for (int i = 0; i < ion.near; i++)
            {
                Site *s = ion.near_sites[i];
                update_site(*s, dt);
                s->last_update = ion.steps;
                if (!settings.rigidBounds)
                    for (int j = 0; j < s->near; j++)
                    {
                        Site *s2 = s->near_sites[j];
                        if (s2->last_update != ion.steps)
                        {
                            update_site(*s2, dt);
                            check_sputter(ion, s);
                        }
                        s2->last_update = ion.steps;
                    }
            }
        }
    }
    else
    {
        update_site(ion, dt);
        for (int i = 0; i < ion.near; i++)
        {
            Site *s = ion.near_sites[i];
            update_site(*s, dt);
            check_sputter(ion, s);
        }
    }
}

double compute_error(Site &site, double dt)
{
    // Error in the associated coordinate.
    double dxp, dyp, dzp;
    // difference in forces between here and destination
    double dpx, dpy, dpz;
    // Find difference in the forces.
    dpx = site.dp_dt[0] - site.dp_dt_t[0];
    dpy = site.dp_dt[1] - site.dp_dt_t[1];
    dpz = site.dp_dt[2] - site.dp_dt_t[2];

    // Error in positions between the two forces.
    dxp = 0.25 * dt * dt * dpx / site.atom->mass;
    dyp = 0.25 * dt * dt * dpy / site.atom->mass;
    dzp = 0.25 * dt * dt * dpz / site.atom->mass;

    return std::max(fabs(dxp), std::max(fabs(dyp), fabs(dzp)));
}

void apply_ion_lattice(Ion &ion, Site *s, double *F_at, double *r_i,
                       double ax, double ay, double az,
                       double dt, bool predicted, double *F_ion)
{
    double dx = 0, dy = 0, dz = 0;
    double fx = 0, fy = 0, fz = 0;

    // Only apply this to the valid sites, invalid ones are ions now.
    if (s->valid == ion.index)
        return;

    // Distances from site to atom
    dx = ax - r_i[0];
    dy = ay - r_i[1];
    dz = az - r_i[2];

    // Lattice/Ion distance
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    if (r < settings.R_MAX && r > 0)
    {
        ion.r_min = std::min(r, ion.r_min);
        s->r_min = std::min(r, s->r_min);

        // Magnitude of force for this location.
        double dV_dr = dVr_dr(r, ion.atom->index, s->atom->index);

        // Potential for this location.
        if (predicted)
        {
            ion.V += Vr_r(r, ion.atom->index, s->atom->index);
        }
        // Scaled by 1/r for converting to cartesian
        dV_dr /= r;

        // Convert from magnitude to components of vector
        // Note, fx = -dV_dr * 2dx, however,
        // we set the force to half of this.
        fx = -dV_dr * dx;
        fy = -dV_dr * dy;
        fz = -dV_dr * dz;

        // Apply to object 1
        F_at[0] += fx;
        F_at[1] += fy;
        F_at[2] += fz;

        // Apply to object 2
        F_ion[0] -= fx;
        F_ion[1] -= fy;
        F_ion[2] -= fz;
        // These are -=, as are opposite direction from on atom
    }
    // No force if ion is on an atom.
    if (r == 0)
    {
        debug_file << "Ion intersected with atom?: " << std::endl;
        ion.write_info();
        s->write_info();
        ion.site_site_intersects++;
    }
}

void apply_lattice_lattice(Site *s, Site *s2, Ion &ion, double *F_at,
                           double atomk, double dt, double ax, double ay, double az,
                           bool predicted, bool recoil, bool useLennardJones, bool doubleCount)
{
    // No Self interactions here!
    if (s == s2)
        return;

    // Only apply this to the valid sites, invalid ones are ions now.
    if (s2->valid == ion.index)
        return;

    // Only apply this to the valid sites, invalid ones are ions now.
    if (s->valid == ion.index)
        return;

    // Force on target.
    double *F_at2;
    double *r2;

    double dx, dy, dz;
    double dx1, dy1, dz1;
    double dx_hat, dy_hat, dz_hat;
    double rr, r;

    // This will be the magnitude of the force.
    double dV_dr = 0;

    if (s2->last_ion != ion.index)
    {
        s2->last_ion = ion.index;
        s2->thermal_seed = ion.thermal_seed;
        s2->reset();
    }

    if (predicted)
    {
        // Get predicted force array.
        F_at2 = s2->dp_dt_t;
        // Get predicted loctations array.
        r2 = s2->r_t;
    }
    else
    {
        F_at2 = s2->dp_dt;
        r2 = s2->r;
    }

    // Compute lattice site separation here.
    dx1 = r2[0] - ax;
    dy1 = r2[1] - ay;
    dz1 = r2[2] - az;
    rr = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

    // Skip all of this if we are out of range.
    if (rr > settings.rr_max)
        goto end;

    r = sqrt(rr);
    s->r_min = std::min(r, s->r_min);

    if (r < 0.1)
    {
        debug_file << "Somehow lattice site on other?" << std::endl;
        s->write_info();
        s2->write_info();
        ion.site_site_intersects++;
        goto end;
    }

    if (useLennardJones)
    {
        // In here we use Lennard Jones forces
        dV_dr = L_J_dV_dr(r, s2->atom->index, s->atom->index);
        if (fabs(dV_dr) > 1e20)
        {
            debug_file << "Somehow large lattice force: "
                       << dV_dr << ", sep: " << r << " "
                       << ion.last_step << std::endl;
            dV_dr = 0;
            ion.site_site_intersects++;
            goto end;
        }
    }
    else
    {
        // In here we use lattice springs->

        // Compute lattice site separation at rest
        dx = s->r_0[0] - s2->r_0[0];
        dy = s->r_0[1] - s2->r_0[1];
        dz = s->r_0[2] - s2->r_0[2];
        rr = dx * dx + dy * dy + dz * dz;
        double l_eq = sqrt(rr);

        double dl = r - l_eq;

        // Check for spring breaking if dl > 0
        if (dl > 0)
        {
            // V = 0.5*k*x^2
            double V_s = 0.5 * atomk * dl * dl;

            // Scale by 1/r^2, this accounts for
            // differences in bond breaking by distance
            V_s /= l_eq * l_eq;

            // Break the spring if too far.
            if (V_s > settings.max_spring_V)
                goto end;
        }
        // F = -kx
        dV_dr = -atomk * dl;
    }

    // This is divided by 2, as the other site will also then apply this.
    // This should be removed if the code is updated so this is only
    // run once per site pair.
    if (doubleCount)
        dV_dr /= 2;

    // These are scaled by 1/r for conversion to
    // the same coordinate system as the fmag was caluclated for
    dx_hat = dx1 / r;
    dy_hat = dy1 / r;
    dz_hat = dz1 / r;
    // Note, fx = -dV_dr * 2dx, however,
    // we set the force to half of this
    // Apply to us
    F_at[0] += dx_hat * dV_dr;
    F_at[1] += dy_hat * dV_dr;
    F_at[2] += dz_hat * dV_dr;
    // Apply to other
    F_at2[0] -= dx_hat * dV_dr;
    F_at2[1] -= dy_hat * dV_dr;
    F_at2[2] -= dz_hat * dV_dr;

end:
    if (recoil)
    {
        if (!predicted)
            // set the predictions
            predict_site_location(*s2, dt);
    }
}

void apply_dynamic_lattice(Site *s, Lattice *lattice, int start, Ion &ion,
                           double *F_at, double atomk, double dt,
                           double ax, double ay, double az,
                           bool predicted, bool recoil, bool useLennardJones)
{
    // Use atomk as k for springs between nearest sites, or lennard jones potentials
    for (int j = 0; j < s->near; j++)
    {
        Site *s2 = s->near_sites[j];
        apply_lattice_lattice(s, s2, ion, F_at, atomk, dt, ax, ay, az,
                              predicted, recoil, useLennardJones, true);
    }
}

void run_hameq(Ion &ion, Lattice *lattice, double dt, bool predicted, double *dr_max)
{
    // Some useful variables.
    double dx = 0, dy = 0, dz = 0;

    // Lattice atom coordintates.
    double ax = 0, ay = 0, az = 0;

    bool recoil = settings.RECOIL;
    bool springs = settings.CORR;
    double atomk = settings.ATOMK;

    bool useEinsteinSprings = settings.useEinsteinSprings;
    bool useLennardJones = settings.useLennardJones;

    // Coordinate of the Ion
    double *r_i;
    r_i = ion.r;

    // Relevant force for this run
    double *F;
    F = ion.dp_dt;

    // if not 0, we are computing for the predicted location.
    if (predicted)
    {
        // Use the one for next time step
        F = ion.dp_dt_t;
        // Use predictions
        r_i = ion.r_t;
    }

    // Update "last integration step" for ion
    ion.last_step++;

    // Reset V
    ion.V = 0;

    // This was set earlier when looking for nearby sites.
    if (ion.near)
    {
        // Check this counter, used for statistics later.
        if (ion.near > ion.max_n)
            ion.max_n = ion.near;

        // Force on the target.
        double *F_at;

        for (int i = 0; i < ion.near; i++)
        {
            // Site near us.
            Site *s = ion.near_sites[i];

            if (s->valid == ion.index)
                continue;

            if (predicted)
            {
                // Get predicted force array.
                F_at = s->dp_dt_t;
                // Get predicted loctation array.
                ax = s->r_t[0];
                ay = s->r_t[1];
                az = s->r_t[2];
            }
            else
            {
                // Get current force array.
                F_at = s->dp_dt;
                // Get current location array.
                ax = s->r[0];
                ay = s->r[1];
                az = s->r[2];
            }

            // Apply site forces.
            // Note that this check is here, as when considering
            // Site neighbours, that inner loop will also apply ion
            // forces if needed.
            if (s->last_step != ion.last_step)
            {
                apply_ion_lattice(ion, s, F_at, r_i, ax, ay, az, dt, predicted, F);
            }

            if (s->hameq_tick == ion.hameq_tick)
                continue;
            s->hameq_tick = ion.hameq_tick;

            // Here we consider lattice-lattice correlations.
            // In the dynamic neighbours case, they will apply to
            // every site in the found list, and is done later
            // otherwise it considers the sites found during the
            // initial setup of the lattice.
            if (springs && !settings.dynamicNeighbours)
            {
                // Displacement from rest
                // This is used in einstein springs
                // but also for determining if
                // we actually want lattice correlations
                dx = s->r[0] - s->r_0[0];
                dy = s->r[1] - s->r_0[1];
                dz = s->r[2] - s->r_0[2];

                // Use einstein springs.
                if (useEinsteinSprings)
                {
                    double dx_sq = dx * dx;
                    double dy_sq = dy * dy;
                    double dz_sq = dz * dz;
                    // V = 0.5kx^2
                    double V_s_x = 0.5 * s->atom->spring[0] * dx_sq;
                    double V_s_y = 0.5 * s->atom->spring[1] * dy_sq;
                    double V_s_z = 0.5 * s->atom->spring[2] * dz_sq;

                    double V_s = V_s_x + V_s_y + V_s_z;
                    // Check if too much spring energy
                    if (settings.max_spring_V > V_s)
                    {
                        // F = -kx
                        F_at[0] -= s->atom->spring[0] * dx;
                        F_at[1] -= s->atom->spring[1] * dy;
                        F_at[2] -= s->atom->spring[2] * dz;
                    }
                    // These springs are only applied if the site is still "near" the original location.
                }
                else
                {
                    // Use atomk as k for springs between nearest sites, or lennard jones potentials
                    for (int j = 0; j < s->near; j++)
                    {
                        Site *s2 = s->near_sites[j];
                        apply_lattice_lattice(s, s2, ion, F_at, atomk, dt, ax, ay, az,
                                              predicted, recoil, useLennardJones, true);
                    }
                }
            }

            if (recoil)
            {
                if (!predicted)
                    // set the predictions
                    predict_site_location(*s, dt);
                else
                    //set account for site errors
                    *dr_max = std::max(compute_error(*s, dt), *dr_max);
            }
        }
        if (springs && settings.dynamicNeighbours)
        {
            for (auto s : lattice->active_sites)
            {
                if (predicted)
                {
                    // Get predicted force array.
                    F_at = s->dp_dt_t;
                    // Get predicted loctation array.
                    ax = s->r_t[0];
                    ay = s->r_t[1];
                    az = s->r_t[2];
                }
                else
                {
                    // Get current force array.
                    F_at = s->dp_dt;
                    // Get current location array.
                    ax = s->r[0];
                    ay = s->r[1];
                    az = s->r[2];
                }
                apply_dynamic_lattice(s, lattice, 0, ion, F_at, atomk, dt, ax, ay, az,
                                      predicted, recoil, useLennardJones);
            }
        }

        // We have now exited the lattice considerations.
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
        // set the predictions
        predict_site_location(ion, dt);
    }
    else
    {
        // Compute error
        *dr_max = std::max(compute_error(ion, dt), *dr_max);
    }
}
