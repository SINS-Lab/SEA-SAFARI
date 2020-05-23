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

void run_hameq(std::vector<Ion *> &ions, Lattice *lattice, double dt, bool predicted, double *dr_max)
{
    Ion *orig = ions[0];
    int num = ions.size();
    int tick = orig->hameq_tick;
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
            apply_ion_lattice(ion, s, F_at, r_i, ax, ay, az, dt, predicted, F);
        }
        run_hameq(ion, lattice, dt, predicted, dr_max);
    }
}
