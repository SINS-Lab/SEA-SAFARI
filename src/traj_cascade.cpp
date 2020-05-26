#include "traj.h"       // header for this
#include "hameq.h"      // the integrator
#include "potentials.h" // the potential used
#include "safio.h"      // general settings
#include <cmath>        // trig functions and sqrt
#include <algorithm>    // std::sort

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

    Ion *orig = ions[0];

    // Reset some values before the integration loop
    std::copy(orig->r, orig->r + 3, orig->r_check);
    orig->last_step = 0;
    orig->steps = 0;
    orig->max_active = 0;
    orig->site_site_intersects = 0;

    // Only bother to call clear if we will use this
    if (settings.saveSputter)
    {
        if (orig->sputter != NULL)
            delete orig->sputter;
        orig->sputter = new Site *[1024];
        orig->sputtered = 0;
    }
    int ion_count = 1;
    double E1 = 0;
    int num = 0;
    int new_num = 0;

    int total_steps = 0;
    int total_threshold = 20 * settings.MAX_STEPS;

start:

    if (total_steps++ >= total_threshold)
    {
        goto end;
    }

    // Reset this to 0.
    dr_max = 0;

    // This is set back true later if needed.
    sort = false;
    reindex = false;

    // Verify time step is in range.
    dt = std::min(std::max(dt, dt_low), dt_high);

    num = ions.size();
    for (int i = 0; i < num; i++)
    {
        Ion *ion = ions[i];
        if (ion->done)
            continue;

        // Increment counter for how many steps we have taken
        ion->steps++;
        ion->reset_forces();

        // Find nearby lattice atoms
        if (settings.dynamicNeighbours)
        {
            update_dynamic_neighbours(ion, ion, lattice, d_search, n_parts,
                                      settings.rr_max, sort or ion->resort,
                                      false, reindex or ion->reindex);
        }
        else
        {
            fill_nearest(ion, ion, lattice, d_search, n_parts,
                         settings.rr_max, sort or ion->resort,
                         false, reindex or ion->reindex);
        }
        // check if we are still in a good state to run.
        ion->done = !validate(ion, E1);
    }

    // Find the forces at the current location
    run_hameq(ions, lattice, dt, false, &dr_max);
    // Find the forces at the next location
    run_hameq(ions, lattice, dt, true, &dr_max);

    new_num = ions.size();
    if (new_num != num)
    {
        reindex = true;
        sort = true;
        // Return back to start of traj run.
        goto start;
    }

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

    orig->time += dt;

    num = ions.size();
    for (int i = 0; i < num; i++)
    {
        Ion *ion = ions[i];
        if (ion->done)
            continue;

        r.set(ion->r_check);
        dr = r - ion->r;

        diff = sqr(dr.v);

        // Reset at either 1/10 of the search distance, or 1/3 the distance to nearest, whichever is larger.
        r_reset = std::max(0.1 * d_search, ion->rr_min_find / 9.0);
        if (diff > r_reset)
        {
            // Update ion location for last check
            std::copy(ion->r, ion->r + 3, ion->r_check);
            // We need to re-sort the lists.
            ion->resort = true;
        }

        ion_count++;
        // Update some parameters for saving.
        ion->log_z = fmin(ion->r[2], ion->log_z);

        // Increment the time step for the ion
        ion->time = orig->time;

        orig->site_site_intersects = orig->site_site_intersects or ion->site_site_intersects;
    }

    if (orig->site_site_intersects)
    {
        goto end;
    }

    if (ion_count == 0)
        goto end;

    // Back up we go.
    goto start;

end:

    int size = ions.size();
    lattice->max_active = std::max(lattice->max_active, size);
    lattice->sum_active += size;
    lattice->count_active++;

    int index = orig->index;
    num = ions.size();
    for (int i = 0; i < num; i++)
    {
        Ion *ion = ions[i];
        ion->index = index;
        lattice->total_hits++;

        if (total_steps >= total_threshold)
        {
            ion->froze = true;
            ion->steps = total_steps;
        }

        // Output data
        if (ion->atom->index == settings.ion.index)
        {
            detector.log(out_file, *ion, lattice,
                         ion->stuck, ion->buried, ion->froze,
                         ion->off_edge, ion->discont, false);
        }
        else
        {
            detector.log(sptr_file, *ion, lattice,
                         ion->stuck, ion->buried, ion->froze,
                         ion->off_edge, ion->discont, false);
        }
        if (i > 0)
            delete ion;
    }
    ions.clear();
}
