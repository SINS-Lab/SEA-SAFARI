
#include "ion.h"
#include "traj.h"
#include "space_math.h"
#include "hameq.h"
#include "potentials.h"
#include "safio.h"
#include <cmath>

bool validate(Ion &ion, bool *buried, bool *off_edge, bool *stuck,
                        bool *froze, bool *left, double E)
{
    //left crystal
    if(ion.r[2] > settings.Z1)
    {
        *left = true;
        return false;
    }
    
    //Ranges of the crystal
    double x_max = settings.AX * settings.RAX;
    double y_max = settings.AY * settings.RAY;
    
    //Fell of the edge
    if(ion.r[0] > x_max || ion.r[0] < -x_max ||
       ion.r[1] > y_max || ion.r[1] < -y_max)
    {
        *off_edge = true;
        return false;
    }
    //Buried
    if(ion.r[2] < -settings.BDIST)
    {
        *buried = true;
        return false;
    }
    //Too many steps
    if(ion.steps > settings.MAX_STEPS)
    {
        *froze = true;
        return false;
    }
    //Stuck in surface
    if(E < settings.SENRGY)
    {
        *stuck = true;
        return false;
    }
    return true;
}

//File output cache related things
int save_index = 0;                //Current index for cache
const int cache_size = 100;        //Size of cache
std::string save_cache[cache_size];//The cache itself

void save(char* buffer)
{
    //If NULL, flush the cache
    if(buffer==NULL)
    {
        //Save all to outfile
        for(int n = 0; n<save_index; n++)
        {
            out_file << save_cache[n];
        }
        //Reset index
        save_index = 0;
    }
    else
    {
        //If not enough in cache, add to cache
        if(save_index < cache_size)
        {
            save_cache[save_index] = buffer;
            save_index++;
        }
        else
        {
            //Save all to outfile
            for(int n = 0; n<save_index; n++)
            {
                out_file << save_cache[n];
            }
            //Output the latest one as well
            out_file << buffer;
            //Reset cache index
            save_index = 0;
        }
    }
}

void traj(Ion &ion, Lattice &lattice, bool log, bool xyz)
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
    //Error in the associated coordinate.
    double dxp, dyp, dzp;
    //difference in forces between here and destination
    double dpx, dpy, dpz;
    //Maximum error on displacment for the ion.
    double dr_max;

    //Energy now and up to 2 steps previously
    double E1, E2, E3;

    //Distance in angstroms to consider far enough moved.
    //After it moves this far, it will re-calculate nearest.
    double r_reset = d_search/10.0;
    //This is set true whenever r_reset's condition is met
    bool sort = true;

    //Multiplier on timestep.
    double change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;
    E1 = E2 = E3 = T;

    //Initialize these to 0
    ion.steps = 0;
    ion.V = 0;
    ion.V_t = 0;

    if(log)
    {
        debug_file << "\n\nStarting Ion Trajectory Output\n\n";
        debug_file << "Recoil: " << (settings.RECOIL?"T":"F") << "\n";
        debug_file << "\nIon:\n";
        ion.write_info();
        traj_file << "x\ty\tz\tpx\tpy\tpz\tt\tn\tT\tV\tE\tnear\tdt\tdr_max\n";
    }
    r.set(ion.r);

start:
    //Find nearby lattice atoms
    ion.fill_nearest(lattice, d_search, n_parts, sort);
    //This is set back true later if needed.
    sort = false;
    //Increment counter for how many steps we have taken
    ion.steps++;
    
    //Verify time step is in range.
    dt = std::min(std::max(dt, dt_low), dt_high);

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, E1);

    //These are our standard exit conditions
    if(froze || buried || stuck || left || off_edge)
    {
        if(log) debug_file << "Exited " << froze
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
    if(fabs(E3 - 2*E2 + E1) > de_fail)
    {
        discont = true;
        goto end;
    }

    //Find the forces at the current location
    run_hameq(ion, lattice, dt, false);
    //Find the forces at the next location
    run_hameq(ion, lattice, dt, true);
    
    //Find difference in the forces.
    dpx = ion.dp_dt[0] - ion.dp_dt_t[0];
    dpy = ion.dp_dt[1] - ion.dp_dt_t[1];
    dpz = ion.dp_dt[2] - ion.dp_dt_t[2];

    //Error in positions between the two forces.
    dxp = 0.25 * dt * dt * dpx / mass;
    dyp = 0.25 * dt * dt * dpy / mass;
    dzp = 0.25 * dt * dt * dpz / mass;

    //TODO decide on whether to include lattice here.

    dr_max = std::max(fabs(dxp), std::max(fabs(dyp), fabs(dzp)));

    //Now we do some checks to see if the timestep needs to be adjusted
    if(dr_max != 0)
    {
        //check if energy has changed too much.
        change = pow(settings.error_scale/dr_max, settings.error_exponent);
        //don't allow speedups of more than 2x, as those can
        //cause major discontinuities later
        if(change >= 2)
            change = 2;
        //If timestep tries to speed up between 1 and 2 times,
        //just leave it where it is.
        if(change > 1 && change < 2)
            change = 1;

        //Large change in energy, try to reduce timestep.
        if(change < .2 && dt > dt_low)
        {
            dt *= change;
            //Make sure it is still at least dt_low.
            if(dt < dt_low)
                dt = dt_low;
            if(log)
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
    if(diff > r_reset)
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
    if(log)
    {
        //This log file is just the ion itself, not the full xyz including the lattice
        sprintf(buffer, "%f\t%f\t%f\t%.3f\t%.3f\t%.3f\t%f\t%d\t%.3f\t%.3f\t%.3f\t%d\t%f\t%.3f\n",
                ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],
                ion.time,ion.steps,T,ion.V,(T+ion.V),ion.near,dt,dr_max);
        traj_file << buffer;
    }
    if(xyz)
    {
        //Update the lattice site momenta...
        for(int i = 0; i<ion.near; i++)
        {
            Site &s = *ion.near_sites[i];
            lattice.sites[s.index] = &s;
            //Flag as nearby
            s.near_check = true;
        }

        //First line for xyz file is number involved.
        int num = lattice.sites.size() + 1;
        xyz_file << num << "\n";

        //Next line is a "comment", we will stuff the time here.
        sprintf(buffer, "%f\t%d\t%.3f\t%.3f\t%d\t%f\t%.3f\n",
                ion.time,ion.steps,T,ion.V,ion.near,dt,dr_max);
        xyz_file << buffer;

        //Next stuff the symbol, position, momentum and mass for the ion itself
        //The ion is given index 0
        sprintf(buffer, "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
                ion.atom->symbol.c_str(),ion.r[0],ion.r[1],ion.r[2], // Sym, x, y, z,
                                         ion.p[0],ion.p[1],ion.p[2], //     px,py,pz,
                                         ion.atom->mass,0,1);        //mass, index, near
        xyz_file << buffer;

        //Then stuff in the entire lattice, why not...
        for(int i = 0; i<num-1; i++)
        {
            Site &s = *lattice.sites[i];
            //Note that this is same format as the ion, index has 1 added to it, as ion is 0
            sprintf(buffer, "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\n",
                    s.atom->symbol.c_str(),s.r[0],s.r[1],s.r[2],     // Sym, x, y, z,
                                           s.p[0],s.p[1],s.p[2],     //     px,py,pz,
                                           s.atom->mass,s.index + 1, //mass, index,
                                           s.near_check);            //Near
            xyz_file << buffer;
            //Reset this flag for next run
            s.near_check = false;
        }
    }

    //Apply differences to timestep for the next loop
    dt *= change;

    //Back up we go.
    goto start;

end:

    if(log)
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
     * -10  : ion got stuck due to image charge.
     * -100 : got stuck, ie E too low
     * -200 : got buried, ie z below -BDIST
     * -300 : froze, ie took too many steps
     * -400 : off edge, ie left crystal via x or y
     * -500 : discont, had a discontinuity in E, so was dropped.
     */ 
    if(stuck)
    {
        theta = 0;
        phi = 90;
        E = -100;
    }
    else if(buried)
    {
        theta = 0;
        phi = 90;
        E = -200;
    }
    else if(froze)
    {
        theta = 0;
        phi = 90;
        E = -300;
    }
    else if(off_edge)
    {
        theta = 0;
        phi = 90;
        E = -400;
    }
    else if(discont)
    {
        theta = 0;
        phi = 90;
        E = -500;
    }
    else
    {
        double px = ion.p[0];
        double py = ion.p[1];
        double pz = ion.p[2];
        psq = sqr(ion.p);
        //Find the momentum at infinity
        if(settings.use_image)
        {
            //Image charge would pull it towards surface, this accounts
            //for that effect.
            pzz = (pz*pz) - (2*mass*Vi_z(settings.Z1, ion.q));
            pzz = pzz < 0 ? -sqrt(-pzz) : sqrt(pzz);
            //Recalulate this, as pz has changed
            //We are fine with pzz being -ve, as that case
            //will be dropped due to ion not escaping.
            psq = pzz + px*px + py*py;
        }
        else
        {
            //No image, so this is the same as it was.
            pzz = pz;
        }
        //Ion is not escaping.
        if(pzz <= 0 || psq <= 0)
        {
            theta = 0;
            phi = 90;
            E = -10;
        }
        else
        {
            //Recalculate E, incase image affected it
            E = 0.5 * psq / mass;

            //calculate theta, depends on pz
            theta = acos(pzz/sqrt(psq)) * 180 / M_PI;

            if(px == 0 && py==0)
            {
                //phi isn't well defined here,
                //so just set it as 90
                phi = 90;
            }
            else
            {
                //Calculate phi, depends on px, py
                phi = atan2(py, px) * 180 / M_PI;
            }
        }
    }

    //Output data, first stuff it in the buffer
    sprintf(buffer, "%f\t%f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\n",
            ion.r_0[0],ion.r_0[1],ion.r_0[2],
            E,theta,phi,
            ion.index,1.0, //TODO make the second 1 be an area based on gridscat depth?
            ion.max_n, ion.r_min, ion.steps, ion.time);
    //Then save it
    save(buffer);

    return;
}