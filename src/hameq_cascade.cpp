#include "hameq.h"

#include "safio.h"      // for settings
#include "traj.h"       // for nearest check method
#include "potentials.h" // for forces/potentials
#include "vec_math.h"   // General maths help
#include "safari.h"     // for exit_fail
#include <cmath>        // sqrt

Ion *makeIonFromSite(Site *s, int index)
{
    Ion *ion = new Ion();

    std::copy(s->r_0, s->r_0 + 3, ion->r_0);
    std::copy(s->r, s->r + 3, ion->r);
    std::copy(s->p, s->p + 3, ion->p);

    s->invalidate(index);
    ion->index = index;
    ion->atom = s->atom;

    return ion;
}

void apply_hameq(std::vector<Ion *> &ions, Lattice *lattice, double dt)
{
    int num = ions.size();
    Ion *orig = ions[0];
    int steps = orig->steps;
    for (int i = 0; i < num; i++)
    {
        Ion &ion = *ions[i];

        if (ion.done)
            continue;

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
                    if (s->valid == ion.index)
                        continue;
                    if (s->last_update == steps)
                        continue;
                    s->last_update = steps;
                    update_site(*s, dt);
                    bool newIon = check_sputter(ion, s);
                    if (newIon)
                    {
                        ions.push_back(makeIonFromSite(s, orig->index));
                    }
                }
            }
            else
            {
                // Otherwise, just consider ones tracked by the ion
                for (int i = 0; i < ion.near; i++)
                {
                    Site *s = ion.near_sites[i];
                    if (s->valid == ion.index)
                        continue;
                    if (s->last_update == steps)
                        continue;
                    s->last_update = steps;
                    update_site(*s, dt);
                    bool newIon = check_sputter(ion, s);
                    if (newIon)
                    {
                        ions.push_back(makeIonFromSite(s, orig->index));
                    }
                    if (!settings.rigidBounds)
                        for (int j = 0; j < s->near; j++)
                        {
                            Site *s2 = s->near_sites[j];
                            if (s2->valid == ion.index)
                                continue;
                            if (s2->last_update != steps)
                            {
                                s2->last_update = steps;
                                update_site(*s2, dt);
                                bool newIon = check_sputter(ion, s2);
                                if (newIon)
                                {
                                    ions.push_back(makeIonFromSite(s2, orig->index));
                                }
                            }
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
                if (s->valid == ion.index)
                    continue;
                if (s->last_update == steps)
                    continue;
                s->last_update = steps;
                update_site(*s, dt);
                bool newIon = check_sputter(ion, s);
                if (newIon)
                {
                    ions.push_back(makeIonFromSite(s, orig->index));
                }
            }
        }
    }
}

void apply_ion_ion(Ion &ion, Ion *s, double *F_at, double *r_i,
                   double ax, double ay, double az,
                   double dt, bool predicted, double *F_ion)
{
    double dx = 0, dy = 0, dz = 0;
    double fx = 0, fy = 0, fz = 0;

    if (s->done || ion.done)
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
    if (r < 0.1)
    {
        debug_file << "Ion intersected with atom?: " << std::endl;
        ion.write_info();
        s->write_info();
        ion.site_site_intersects++;
    }
}

void run_hameq(std::vector<Ion *> &ions, Lattice *lattice, double dt, bool predicted, double *dr_max)
{
    Ion *orig = ions[0];
    int num = ions.size();
    int tick = orig->hameq_tick;
    int hameqstep = orig->last_step;
    for (int i = 0; i < num; i++)
    {
        Ion &ion = *ions[i];
        ion.steps = orig->steps;

        if (ion.done)
            continue;
        ion.hameq_tick = tick;
        // Find the forces at the current location
        ion.hameq_tick++;

        // Run for the ion-ion interaction
        // Force on the target.
        double *F_at;

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

        // Lattice atom coordintates.
        double ax = 0, ay = 0, az = 0;
        for (int j = i + 1; j < num; j++)
        {
            // Site near us.
            Ion *s = ions[j];

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
            apply_ion_ion(ion, s, F_at, r_i, ax, ay, az, dt, predicted, F);
        }
        ion.last_step = hameqstep++;
        run_hameq(ion, lattice, dt, predicted, dr_max);
    }
    orig->last_step = hameqstep;
}
