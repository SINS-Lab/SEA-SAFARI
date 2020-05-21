#include "traj.h"       // header for this
#include "hameq.h"      // the integrator
#include "potentials.h" // the potential used
#include "safio.h"      // general settings
#include <cmath>        // trig functions and sqrt
#include <algorithm>    // std::sort

void update_dynamic_neighbours(std::vector<Ion *> &ions, Ion *ion_ptr, Site *site, Lattice *lattice, int radius, int target_num, double max_rr, bool re_sort, bool updateCells)
{
}

int fill_nearest(std::vector<Ion *> &ions, Ion *ion_ptr, Site *site, Lattice *lattice, int radius, int target_num, double max_rr, bool re_sort, bool updateCells)
{
    return 0;
}

void traj(std::vector<Ion *> &ions, Lattice *lattice, bool &log, bool &xyz, Detector &detector)
{
    if (!settings.useEinsteinSprings)
    {
        // We need to call this before each run, as the run can
        // shuffle the sites between cells.
        lattice->init_springs(settings.neighbour_count);
    }

    // Get some constants for the loop

    // Radius of atom search
    const int d_search = settings.DIST_SEARCH;
    // Max atoms to find in the search
    const int n_parts = settings.NPART;
    // Lowest allowed time step
    const double dt_low = settings.DELLOW;
    // Highest allowed time step
    const double dt_high = settings.DELT0;

    // Time step, initialized at 0.1
    double dt = 0.1;

    // Parameters for checking if things need re-calculating
    // Location of last nearby update
    std::vector<Vec3d *> rs;
    Vec3d r;
    // change in location since last nearby update
    Vec3d dr;
    // Distancesq we travelled since last nearby update
    double diff = 0;
    // Maximum error on displacment during integration.
    double dr_max = 0;

    // Distance in angstroms to consider far enough moved.
    // After it moves this far, it will re-calculate nearest.
    double r_reset = d_search / 10.0;
    // This is set true whenever r_reset's condition is met
    bool sort = true;
    bool reindex = true;

    // Multiplier on timestep.
    double change;

    Ion &orig = *ions[0];

    // Reset some values before the integration loop
    r.set(orig.r);
    orig.last_step = 0;
    orig.max_active = 0;
    orig.site_site_intersects = 0;

    // Only bother to call clear if we will use this
    if (settings.saveSputter)
    {
        if (orig.sputter != NULL)
            delete orig.sputter;
        orig.sputter = new Site *[1024];
        orig.sputtered = 0;
    }
    int ion_count = 1;
    double E1 = 0;

start:

    // Reset this to 0.
    dr_max = 0;

    // This is set back true later if needed.
    sort = false;
    reindex = false;
    // Increment counter for how many steps we have taken
    orig.steps++;

    // Verify time step is in range.
    dt = std::min(std::max(dt, dt_low), dt_high);

    int num = ions.size();
    for (int i = 0; i < num; i++)
    {
        Ion &ion = *ions[i];
        if (ion.done)
            continue;

        // Find nearby lattice atoms
        if (settings.dynamicNeighbours)
        {
            update_dynamic_neighbours(&ion, &ion, lattice, d_search, n_parts,
                                      settings.rr_max, sort or ion.resort, 
                                      false, reindex or ion.reindex);
        }
        else
        {
            fill_nearest(&ion, &ion, lattice, d_search, n_parts,
                         settings.rr_max, sort or ion.resort, 
                         false, reindex or ion.reindex);
        }
        // check if we are still in a good state to run.
        validate(ion, E1);

        // These are our standard exit conditions
        if (ion.froze || ion.buried || ion.stuck || ion.left || ion.off_edge)
        {
            ion.done = true;
            continue;
        }
    }

    // Find the forces at the current location
    run_hameq(ions, lattice, dt, false, &dr_max);
    // Find the forces at the next location
    run_hameq(ions, lattice, dt, true, &dr_max);

    // Now we do some checks to see if the timestep needs to be adjusted
    if (dr_max != 0)
    {
        // check if energy has changed too much.
        change = pow(settings.error_scale / dr_max, settings.error_exponent);
        // don't allow speedups of more than 2x, as those can
        // cause major discontinuities later
        if (change >= 2)
            change = 2;

        // Large change in energy, try to reduce timestep.
        if (change < .2 && dt > dt_low)
        {
            dt *= change;
            // Make sure it is still at least dt_low.
            if (dt < dt_low)
                dt = dt_low;
            if (log)
            {
                debug_file << "change_down: " << change << "\n";
            }
            // Reset the last index the ion saw,
            // this forces a re-check of nearby atoms
            reindex = true;
            sort = true;
            // Return back to start of traj run.
            goto start;
        }
    }
    else
    {
        // By default, try to increase the time step.
        // This ensures things speed up again after leaving.
        change = 2;
    }

    apply_hameq(ions, lattice, dt);

    // Apply differences to timestep for the next loop
    dt *= change;

    ion_count = 0;

    orig.time += dt;

    num = ions.size();
    for (int i = 0; i < num; i++)
    {
        Ion &ion = *ions[i];
        if (ion.done)
            continue;

        Vec3d *r = NULL;
        int rsc = rs.size();
        if (rsc < i)
            r = rs[i];
        else
        {
            r = new Vec3d();
            rs.push_back(r);
        }
        dr = *r - ion.r;
        
        diff = sqr(dr.v);

        // Reset at either 1/10 of the search distance, or 1/3 the distance to nearest, whichever is larger.
        r_reset = std::max(0.1 * d_search, ion.rr_min_find / 9.0);
        if (diff > r_reset)
        {
            // Update ion location for last check
            r->set(ion.r);
            // We need to re-sort the lists.
            ion.resort = true;
        }

        ion_count++;
        // Update some parameters for saving.
        ion.log_z = fmin(ion.r[2], ion.log_z);

        // Increment the time step for the ion
        ion.time = orig.time;

        orig.site_site_intersects = orig.site_site_intersects or ion.site_site_intersects;
    }

    if (orig.site_site_intersects)
    {
        goto end;
    }

    if (ion_count == 0)
        goto end;

    // Back up we go.
    goto start;

end:

    Ion &ion = *ions[0];

    // Log the sputtered data first, the ion's index is set as a flag in the logging process
    // We do not want to save if it was discontinous, or froze!
    if (settings.saveSputter && !ion.discont && !ion.froze)
    {
        double E = 0;
        if (ion.stuck)
        {
            E = -100;
        }
        else if (ion.buried)
        {
            E = -200;
        }
        else if (ion.froze)
        {
            E = -300;
        }
        else if (ion.off_edge)
        {
            E = -400;
        }
        else if (ion.discont)
        {
            E = -500;
        }
        else if (ion.site_site_intersects)
        {
            E = -900;
        }

        double x_0 = ion.r_0[0];
        double y_0 = ion.r_0[1];

        double x_1 = ion.r[0];
        double y_1 = ion.r[1];

        double x_min, x_max, y_min, y_max;

        double dr = settings.R_MAX / 2;

        x_min = std::min(x_0, x_1) - dr;
        y_min = std::min(y_0, y_1) - dr;
        x_max = std::max(x_0, x_1) + dr;
        y_max = std::max(y_0, y_1) + dr;

        for (int i = 0; i < ion.sputtered; i++)
        {
            Site *s = ion.sputter[i];
            // Cache this value
            int oldIndex = s->index;

            x_0 = s->r_0[0];
            y_0 = s->r_0[1];

            // If we are too close to the edge of the detection bounds,
            // Just skip this particle, as is probably not valid anyway.
            if (x_0 < x_min || x_0 > x_max)
                continue;
            if (y_0 < y_min || y_0 > y_max)
                continue;

            // This means we are actually going downwards, not valid sputter!
            if (s->p[2] < 0 or s->r[2] < settings.Z1 / 4)
                continue;

            // Copy some values over from the ion
            s->steps = ion.steps;
            // Recorde the ion's error flag as the wieght instead
            s->weight = E;
            s->index = ion.index;
            s->Eerr_max = ion.Eerr_max;
            s->time = ion.time;
            detector.log(sptr_file, *s, lattice, false, false, false, false, false, true);
            // Revert the change to index.
            s->index = oldIndex;
        }
    }
    int index = ion.index;
    num = ions.size();
    for (int i = 0; i < num; i++)
    {
        ion = *ions[i];
        ion.index = index;
        lattice->sum_active += ion.max_active;
        lattice->count_active++;
        // Output data
        detector.log(out_file, ion, lattice,
                     ion.stuck, ion.buried, ion.froze,
                     ion.off_edge, ion.discont, false);
    }

    return;
}
