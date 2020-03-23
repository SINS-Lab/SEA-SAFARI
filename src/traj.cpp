
#include "ion.h"
#include "traj.h"
#include "space_math.h"
#include "hameq.h"
#include "potentials.h"
#include "safio.h"
#include <cmath>
#include <algorithm> // std::sort

int fill_nearest(Ion *ion_ptr, Site &site, Lattice &lattice, int radius, int target_num, double max_rr, bool re_sort, bool updateCells)
{
    //Initial locations are where ion is.
    double cell_x = site.r[0];
    double cell_y = site.r[1];
    double cell_z = site.r[2];

    int pos_hash = to_hash(cell_x, cell_y, cell_z);

    if (updateCells && pos_hash != site.cell_number)
    {
        Cell* from = lattice.get_cell(site.cell_number);
        Cell* to = lattice.make_cell(pos_hash);
        moveSite(site, from, to);
    }

    //New site, so we need to re-calculate things
    if (pos_hash != site.last_index)
    {
        //We need to resort this if we find new index.
        re_sort = true;

        //Update last index checked.
        site.last_index = pos_hash;

        //Number in current cell being checked.
        int num;

        //Reset total nearby, as we are now looking for them.
        site.total_near = 0;

        //Initialize this big.
        site.rr_min_find = 1e6;

        //Location of mask
        Vec3d loc;

        //Initialize this to max allowed distance.
        double rr_min = max_rr;

        //volume of the mask to check.
        int nmax = pow(2 * radius + 1, 3);

        int last_check = -1;

        int check_index = ion_ptr == NULL ? -site.index : site.index;
        int check_step = ion_ptr == NULL ? -2 : ion_ptr->steps;

        //Loop over the mask, this is a radial loop.
        for (int n = 0; n < nmax; n++)
        {
            //Gets x,y,z for the centered cube.
            index_to_loc(n, loc);
            //Translates to ion coordinate
            double x = loc[0] * CELL_SIZE + cell_x;
            double y = loc[1] * CELL_SIZE + cell_y;
            double z = loc[2] * CELL_SIZE + cell_z;

            pos_hash = to_hash(x, y, z);
            //Due to mask and cell sizes differing, this
            //can happen easily, so lets skip if it did.
            if (last_check == pos_hash)
                continue;
            last_check = pos_hash;

            Cell *cell = lattice.get_cell(x, y, z);
            //No cell here? skip.
            if (cell == NULL)
                continue;
            //Reset check stamp if new ion.
            if (cell->ion_stamp != check_index)
                cell->check_stamp = -1;
            //Already seen this cell this tick? skip.
            if (cell->check_stamp == check_step)
                continue;
            //Stamp cell so it gets skipped if seen again.
            cell->check_stamp = check_step;
            cell->ion_stamp = check_index;
            //Get number of sites from the cell
            num = cell->num;
            //Check each site in this cell.
            for (int i = 0; i < num; i++)
            {
                Site *s = cell->sites[i];
                //If some other ion has seen the site, reset it here.
                //This reset puts it back to where it should be.
                if (ion_ptr != NULL && s->last_ion != ion_ptr->index)
                {
                    s->last_ion = ion_ptr->index;
                    s->thermal_seed = ion_ptr->thermal_seed;
                    s->reset();
                }
                //In this case, we want to make sure we are not including self.
                if (ion_ptr == NULL)
                {
                    if (s->index == site.index)
                        continue;
                }
                //Check if site is close enough
                double rr = diff_sqr(site.r, s->r);
                if (rr > max_rr)
                    continue;
                rr_min = std::min(rr_min, rr);
                site.rr_min_find = std::min(rr, site.rr_min_find);
                //Add the site to our tracked sites.
                site.near_sites[site.total_near] = s;
                site.total_near++;
                //If we have enough, goto end.
                if (site.total_near > 255)
                    goto end;
            }
        }
    }
end:
    if (re_sort)
    {
        //Reset number nearby.
        site.near = 0;
        //Sort the sites, and update near to min of near or target_num
        std::sort(site.near_sites, site.near_sites + site.total_near, [&site](Site *a, Site *b) { return site.compare_sites(a, b); });
        //Update nearby number, taking min of these two.
        site.near = std::min(site.total_near, target_num);
    }

    if (ion_ptr != NULL && settings.useLennardJones)
    {
        //Update each nearby site as well
        for (int i = 0; i < site.near; i++)
        {
            Site &s2 = *site.near_sites[i];
            fill_nearest(NULL, s2, lattice, radius, target_num, max_rr, true, true);
        }
    }
    //Sets this to 0, so that the max check later is fine.
    if (site.rr_min_find == 1e6)
        site.rr_min_find = 0;
    return site.near;
}

bool validate(Ion &ion, bool *buried, bool *off_edge, bool *stuck,
              bool *froze, bool *left, double &E)
{
    //left crystal
    if (ion.r[2] > settings.Z1)
    {
        *left = true;
        return false;
    }

    //Ranges of the crystal
    double x_max = settings.AX * settings.RAX;
    double y_max = settings.AY * settings.RAY;

    //Fell of the edge
    if (ion.r[0] > x_max || ion.r[0] < -x_max ||
        ion.r[1] > y_max || ion.r[1] < -y_max)
    {
        *off_edge = true;
        return false;
    }
    //Buried
    if (ion.r[2] < -settings.BDIST)
    {
        *buried = true;
        return false;
    }
    //Too many steps
    if (ion.steps > settings.MAX_STEPS)
    {
        *froze = true;
        return false;
    }
    //Stuck in surface
    if (E < settings.SENRGY)
    {
        *stuck = true;
        return false;
    }
    return true;
}

void log_xyz(Ion &ion, Lattice &lattice, int &lattice_num, char *buffer)
{
    //Update the lattice site momenta...
    for (int i = 0; i < ion.near; i++)
    {
        Site &s = *ion.near_sites[i];
        lattice.sites[s.index] = &s;
        //Flag as nearby
        s.near_check = true;
    }

    //First line for xyz file is number involved.
    xyz_file << lattice_num << "\n";

    //Next line is a "comment", we will stuff the time here.
    xyz_file << ion.time << "\n";

    //Next stuff the symbol, position, momentum and mass for the ion itself
    //The ion is given index 0
    sprintf(buffer, "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
            ion.atom->symbol.c_str(), ion.r[0], ion.r[1], ion.r[2], // Sym, x, y, z,
            ion.p[0], ion.p[1], ion.p[2],                           //     px,py,pz,
            ion.atom->mass, 0, 1);                                  //mass, index, near
    xyz_file << buffer;

    int num = lattice.sites.size();
    //Then stuff in the entire lattice, why not...
    for (int i = 0; i < num; i++)
    {
        Site &s = *lattice.sites[i];

        //Check if we want to ignore the lattice site
        if (settings.SCAT_TYPE)
        {
            //check out of bounds, from r_0
            if (s.r_0[0] < lattice.xyz_bounds[0]     //x
                || s.r_0[0] > lattice.xyz_bounds[3]  //x
                || s.r_0[1] < lattice.xyz_bounds[1]  //y
                || s.r_0[1] > lattice.xyz_bounds[4]  //y
                || s.r_0[2] < lattice.xyz_bounds[2]  //z
                || s.r_0[2] > lattice.xyz_bounds[5]) //z
                continue;
        }

        //Note that this is same format as the ion, index has 1 added to it, as ion is 0
        sprintf(buffer, "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
                s.atom->symbol.c_str(), s.r[0], s.r[1], s.r[2], // Sym, x, y, z,
                s.p[0], s.p[1], s.p[2],                         //     px,py,pz,
                s.atom->mass, s.index + 1,                      //mass, index,
                s.near_check);                                  //Near
        xyz_file << buffer;
        //Reset this flag for next run
        s.near_check = false;
    }
}

void traj(Ion &ion, Lattice &lattice, bool &log, bool &xyz,
          Detector &detector)
{
    //Get some constants for the loop

    //Radius of atom search
    const int d_search = settings.DIST_SEARCH;
    //Max atoms to find in the search
    const int n_parts = settings.NPART;
    //Energy difference to consider a failure
    const int de_fail = settings.FAILED_DE;
    //Lowest allowed time step
    const double dt_low = settings.DELLOW;
    //Highest allowed time step
    const double dt_high = settings.DELT0;
    //Ion mass
    const double mass = ion.atom->mass;

    //Time step, initialized at 0.1
    double dt = 0.1;
    //Magnitude squared of the ion momentum
    double psq;

    //exit conditions

    //Below -BDIST
    bool buried = false;
    //Took longer than MAX_STEPS
    bool froze = false;
    //Energy lower than SENRGY
    bool stuck = false;
    //Left RAX/RAY
    bool off_edge = false;
    //Too much change in energy over last 2 ticks
    bool discont = false;
    //Above Z1
    bool left = false;

    //Kinetic Energy
    double T;

    //Parameters for checking if things need re-calculating
    //Location of last nearby update
    Vec3d r;
    //change in location since last nearby update
    Vec3d dr;
    //Distancesq we travelled since last nearby update
    double diff = 0;
    //Maximum error on displacment during integration.
    double dr_max = 0;

    //Energy now and up to 2 steps previously
    double E1, E2, E3;

    //Distance in angstroms to consider far enough moved.
    //After it moves this far, it will re-calculate nearest.
    double r_reset = d_search / 10.0;
    //This is set true whenever r_reset's condition is met
    bool sort = true;

    //Used for xyz output
    int lattice_num = 1;

    //Multiplier on timestep.
    double change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;
    E1 = E2 = E3 = T;

    if (log)
    {
        debug_file << "\n\nStarting Ion Trajectory Output\n\n";
        debug_file << "Recoil: " << (settings.RECOIL ? "T" : "F") << "\n";
        debug_file << "\nIon:\n";
        ion.write_info();
        traj_file << "x\ty\tz\tpx\tpy\tpz\tt\tn\tT\tV\tE\tnear\tdt\tdr_max\n";
        //Log initial state
        sprintf(buffer, "%f\t%f\t%f\t%.3f\t%.3f\t%.3f\t%f\t%d\t%.3f\t%.3f\t%.3f\t%d\t%f\t%.3f\n",
                ion.r[0], ion.r[1], ion.r[2], ion.p[0], ion.p[1], ion.p[2],
                ion.time, ion.steps, T, ion.V, (T + ion.V), ion.near, dt, dr_max);
        traj_file << buffer;
    }

    if (xyz)
    {
        int num = lattice.sites.size();
        //Here we initialize lattice_num, it was set to 1 earlier for the ion.
        //Now it gets incremented for each valid site
        for (int i = 0; i < num; i++)
        {
            Site &s = *lattice.sites[i];
            //Check if we want to ignore the lattice site
            if (settings.SCAT_TYPE)
            {
                //check out of bounds, from r_0
                if (s.r_0[0] < lattice.xyz_bounds[0]     //x
                    || s.r_0[0] > lattice.xyz_bounds[3]  //x
                    || s.r_0[1] < lattice.xyz_bounds[1]  //y
                    || s.r_0[1] > lattice.xyz_bounds[4]  //y
                    || s.r_0[2] < lattice.xyz_bounds[2]  //z
                    || s.r_0[2] > lattice.xyz_bounds[5]) //z
                    continue;
            }
            lattice_num++;
        }
        //Log initial state
        log_xyz(ion, lattice, lattice_num, buffer);
    }

    r.set(ion.r);
    ion.last_step = 0;

start:
    //Reset this to 0.
    dr_max = 0;
    //Find nearby lattice atoms
    fill_nearest(&ion, ion, lattice, d_search, n_parts, settings.rr_max, sort, false);
    //This is set back true later if needed.
    sort = false;
    //Increment counter for how many steps we have taken
    ion.steps++;

    //Verify time step is in range.
    dt = std::min(std::max(dt, dt_low), dt_high);

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, E1);

    //These are our standard exit conditions
    if (froze || buried || stuck || left || off_edge)
    {
        if (log)
            debug_file << "Exited " << froze
                       << " " << buried
                       << " " << stuck
                       << " " << left
                       << " " << off_edge
                       << std::endl;
        goto end;
    }

    //Check our energy trackers for discontinuities
    //This essentially gives the value of the second derivative,
    //Showing any major jumps in energy of the particle.
    //This value should average around 0.01, so if larger than 50,
    //Then we have a major jump.
    ion.Eerr_max = std::max(ion.Eerr_max, fabs(E3 - 2 * E2 + E1));
    if (ion.Eerr_max > de_fail)
    {
        discont = true;
        goto end;
    }

    //Find the forces at the current location
    run_hameq(ion, lattice, dt, false, &dr_max);
    //Find the forces at the next location
    run_hameq(ion, lattice, dt, true, &dr_max);

    //Now we do some checks to see if the timestep needs to be adjusted
    if (dr_max != 0)
    {
        //check if energy has changed too much.
        change = pow(settings.error_scale / dr_max, settings.error_exponent);
        //don't allow speedups of more than 2x, as those can
        //cause major discontinuities later
        if (change >= 2)
            change = 2;
        //If timestep tries to speed up between 1 and 2 times,
        //just leave it where it is.
        if (change > 1 && change < 2)
            change = 1;

        //Large change in energy, try to reduce timestep.
        if (change < .2 && dt > dt_low)
        {
            dt *= change;
            //Make sure it is still at least dt_low.
            if (dt < dt_low)
                dt = dt_low;
            if (log)
            {
                debug_file << "change_down: " << change << std::endl;
            }
            //Reset the last index the ion saw,
            //this forces a re-check of nearby atoms
            ion.last_index = -1;
            sort = true;
            //Return back to start of traj run.
            goto start;
        }
    }
    else
    {
        //By default, try to increase the time step.
        //This ensures things speed up again after leaving.
        change = 2;
    }

    //Apply changes, this updates lattice loctations.
    apply_hameq(ion, lattice, dt);

    //check if we have gone too far, and need to re-calculate nearby atoms
    dr = r - ion.r;
    diff = sqr(dr.v);
    //Reset at either 1/10 of the search distance, or 1/3 the distance to nearest, whichever is larger.
    r_reset = std::max(0.1 * d_search, ion.rr_min_find / 9.0);
    if (diff > r_reset)
    {
        //Update ion location for last check
        r.set(ion.r);
        //We need to re-sort the lists.
        sort = true;
    }

    //Update some parameters for saving.
    ion.r_0[2] = fmin(ion.r[2], ion.r_0[2]);

    //Increment the time step for the ion
    ion.time += dt;

    //Update kinetic energy.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;

    //Update energy trackers, propogate them backwards.
    E3 = E2;
    E2 = E1;
    E1 = T + ion.V;

    //Log things if needed
    if (log)
    {
        //This log file is just the ion itself, not the full xyz including the lattice
        sprintf(buffer, "%f\t%f\t%f\t%.3f\t%.3f\t%.3f\t%f\t%d\t%.3f\t%.3f\t%.3f\t%d\t%f\t%.3f\n",
                ion.r[0], ion.r[1], ion.r[2], ion.p[0], ion.p[1], ion.p[2],
                ion.time, ion.steps, T, ion.V, (T + ion.V), ion.near, dt, dr_max);
        traj_file << buffer;
    }
    if (xyz)
    {
        log_xyz(ion, lattice, lattice_num, buffer);
    }

    //Apply differences to timestep for the next loop
    dt *= change;

    //Back up we go.
    goto start;

end:

    if (log)
    {
        debug_file << "\n\nFinished Traj Run\n\n";
        debug_file << "\nIon:\n";
        ion.write_info();
    }

    double E = T;
    //z-momentum squared, exit theta, exit phi
    double pzz = 0, theta, phi;
    /**
     * Checks failure flags, and sets energy accordingly
     * 
     * In all failure cases, theta = 0, phi = 90
     * 
     * Energy failure flags as follows:
     * 
     * -5   : Ion exited at phi out of absolute max detector bounds
     * -10  : ion got stuck due to image charge.
     * -100 : got stuck, ie E too low
     * -200 : got buried, ie z below -BDIST
     * -300 : froze, ie took too many steps
     * -400 : off edge, ie left crystal via x or y
     * -500 : discont, had a discontinuity in E, so was dropped.
     * 
     * Note that -5 condition is only applied if settings.detector_type is greater than 0
     */
    if (stuck)
    {
        theta = 0;
        phi = 90;
        E = -100;
        lattice.stuck_num++;
    }
    else if (buried)
    {
        theta = 0;
        phi = 90;
        E = -200;
        lattice.buried_num++;
    }
    else if (froze)
    {
        theta = 0;
        phi = 90;
        E = -300;
        lattice.froze_num++;
    }
    else if (off_edge)
    {
        theta = 0;
        phi = 90;
        E = -400;
        lattice.left_num++;
    }
    else if (discont)
    {
        theta = 0;
        phi = 90;
        E = -500;
        lattice.err_num++;
    }
    else
    {
        double px = ion.p[0];
        double py = ion.p[1];
        double pz = ion.p[2];
        psq = sqr(ion.p);
        //Find the momentum at infinity
        if (settings.use_image)
        {
            //Image charge would pull it towards surface, this accounts
            //for that effect.
            pzz = (pz * pz) + (2 * mass * Vi_z(settings.Z1, ion.q));
            pzz = pzz < 0 ? -sqrt(-pzz) : sqrt(pzz);
            //Recalulate this, as pz has changed
            //We are fine with pzz being -ve, as that case
            //will be dropped due to ion not escaping.
            psq = pzz * pzz + px * px + py * py;
        }
        else
        {
            //No image, so this is the same as it was.
            pzz = pz;
        }
        double p = sqrt(psq);
        //Ion is not escaping.
        if (pzz <= 0 || p <= 0)
        {
            theta = 0;
            phi = 90;
            E = -10;
            lattice.trapped_num++;
        }
        else
        {
            //Recalculate E, incase image affected it
            E = 0.5 * psq / mass;

            //calculate theta, depends on pz
            theta = acos(pzz / p) * 180 / M_PI;

            if (px == 0 && py == 0)
            {
                //phi isn't well defined here,
                //so just set it as same as incoming
                phi = settings.PHI0;
            }
            else
            {
                //Calculate phi, depends on px, py
                phi = atan2(py, px) * 180 / M_PI;
            }
        }
    }

    //Output data
    detector.log(ion, lattice, E, theta, phi);
    return;
}