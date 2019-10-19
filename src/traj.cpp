
#include "ion.h"
#include "traj.h"
#include "space_math.h"
#include "hameq.h"
#include "potentials.h"
#include <cmath>

bool validate(Ion &ion, bool *buried, bool *off_edge, bool *stuck,
                        bool *froze, bool *left, double E)
{
    //left crystal
    if(ion[2] > settings.Z1)
    {
        *left = true;
        return false;
    }
    
    //Ranges of the crystal
    double x_max = settings.AX * settings.RAX;
    double y_max = settings.AY * settings.RAY;
    
    //Fell of the edge
    if(ion[0] > x_max || ion[0] < -x_max ||
       ion[1] > y_max || ion[1] < -y_max)
    {
        *off_edge = true;
        return false;
    }
    //Buried
    if(ion[2] < -settings.BDIST)
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
int save_index = 0;
const int cache_size = 100;
std::string save_cache[cache_size];

//This saves the current site, in batches of cache_size
//If buffer is NULL, this saves all remaining trajectories.
void save(char* buffer)
{
    if(buffer==NULL)
    {
        for(int n = 0; n<save_index; n++)
        {
            out_file << save_cache[n];
        }
        save_index = 0;
    }
    else
    {
        if(save_index < cache_size)
        {
            save_cache[save_index] = buffer;
            save_index++;
        }
        else
        {
            for(int n = 0; n<save_index; n++)
            {
                out_file << save_cache[n];
            }
            save_index = 0;
        }
    }
}

void traj(Ion &ion, Lattice &lattice, bool log, bool xyz)
{
    //Time step, initialized at 0.1
    double dt = 0.1;
    //Magnitude squared of the ion momentum
    double psq;
    //Ion mass
    double mass = ion.atom.mass;

    //exit conditions
    bool buried = false;
    bool froze = false;
    bool stuck = false;
    bool off_edge = false;
    bool discont = false;
    bool left = false;

    //Kinetic Energy
    double T;
    
    //Parameters for checking if things need re-calculating
    //Location of last nearby update
    Vec3d r;
    //change in location since last nearby update
    Vec3d dr;
    //Distancesq we travelled since last nearby update
    double diff;
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
    double r_reset = 1;

    //Multiplier on timestep.
    double change;

    //Used for printing output.
    char buffer[200];

    //Some initial conditions.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;
    E1 = E2 = E3 = T;

    ion.steps = 0;
    ion.V = 0;
    ion.V_t = 0;

    if(log)
    {
        debug_file << "\n\nStarting Ion Trajectory Output\n\n";
        debug_file << "Recoil: " << (settings.RECOIL?"T":"F") << "\n";
        debug_file << "\nIon:\n";
        ion.write_info();
        traj_file << "x\ty\tz\tpx\tpy\tpz\tt\tn\tT\tV\tE\tnear\tdt\tdr_max\tdV\n";
    }
    r.set(ion.r);

start:
    //Find nearby lattice atoms
    ion.fill_nearest(lattice, 4, settings.NPART);
    //Increment counter for how many steps we have taken
    ion.steps++;
    
    //Verify time step is in range.
    dt = std::min(std::max(dt, settings.DELLOW), settings.DELT0);

    //check if we are still in a good state to run.
    validate(ion, &buried, &off_edge, &stuck, &froze, &left, E1);

    //These are our standard exit conditions
    if(froze || buried || stuck || left || off_edge)
    {
        if(log) debug_file << "Exited" << std::endl;
        goto end;
    }

    //Check our energy trackers for discontinuities
    //This essentially gives the value of the second derivative,
    //Showing any major jumps in energy of the particle.
    //This value should average around 0.01, so if larger than 10,
    //Then we have a major jump.
    if(fabs(E3 - 2*E2 + E1) > 10)
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
        change = pow(settings.ABSERR/dr_max, settings.DEMAX);
        if(change >= 2)
            change = 2;
        if(change > 1 && change < 2)
            change = 1;

        //Large change in energy, try to reduce timestep.
        if(change < .2 && dt > settings.DELLOW)
        {
            dt *= change;
            if(dt < settings.DELLOW)
                dt = settings.DELLOW;
            if(log)
            {
                debug_file << "change_down: " << change << std::endl;
            }
            ion.last_index = -1;
            goto start;
        }
    }
    else
    {
        change = 2;
    }

    //Apply changes, this updates lattice loctations.
    apply_hameq(ion, lattice, dt);
    
    //check if we have gone too far, and need to re-calculate nearby atoms
    dr = r - ion.r;
    diff = sqr(dr.v);
    //Reset at either 1A, or 1/3 the distance to nearest, whichever is larger.
    r_reset = std::max(1.0, ion.rr_min_find / 9.0);
    if(diff > r_reset)
    {
        //Update ion location for last check
        r.set(ion.r);
        //Reset the index considered as "last found",
        //This forces it to re-calculate nearest atoms
        ion.last_index = -1;
    }

    //Update some parameters for saving.
    ion.r_0[2] = fmin(ion.r[2], ion.r_0[2]);
    ion.time += dt;

    //Update kinetic energy.
    psq = sqr(ion.p);
    T = psq * 0.5 / mass;

    //Update energy trackers, propogate the energy backwards.
    E3 = E2;
    E2 = E1;
    E1 = T + (ion.V + ion.V_t) / 2;

    //Log things if needed
    if(log)
    {
        //Average potential energy
        double V = (ion.V + ion.V_t) / 2;
        double dV = ion.V_t - ion.V;
        //This log file is just the ion itself, not the full xyz including the lattice
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\n",
                ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],
                ion.time,ion.steps,T,V,(T+V),ion.near,dt,dr_max,dV);
        traj_file << buffer;
    }
    if(xyz)
    {
        double V = (ion.V + ion.V_t) / 2;

        //Update the lattice site momenta...
        for(int i = 0; i<ion.near; i++)
        {
            Site *s = ion.near_sites[i];
            lattice.sites[s->index] = ion.near_sites[i];
        }

        //First line for xyz file is number involved.
        int num = lattice.sites.size() + 1;
        xyz_file << num << "\n";

        //Next line is a "comment", we will stuff the time here.
        sprintf(buffer, "%f\t%d\t%f\t%f\t%d\t%f\t%f\n",
                ion.time,ion.steps,T,V,ion.near,dt,dr_max);
        xyz_file << buffer;

        //Next stuff the symbol, position, momentum and mass for the ion itself
        sprintf(buffer, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                ion.atom.symbol.c_str(),ion.r[0],ion.r[1],ion.r[2],ion.p[0],ion.p[1],ion.p[2],ion.atom.mass);
        xyz_file << buffer;

        //Then stuff in the entire lattice, why not...
        for(int i = 0; i<num-1; i++)
        {
            Site s = *lattice.sites[i];
            //Note that this is same format as the ion, except also contains the index.
            sprintf(buffer, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",
                    s.atom.symbol.c_str(),s.r[0],s.r[1],s.r[2],s.p[0],s.p[1],s.p[2],s.atom.mass,s.index);
            xyz_file << buffer;
        }
    }

    //Apply differences to timestep;
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
    double pp = 0, pzz = 0, theta, phi;
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
        if(settings.IMAGE)
        {
            pp = (pz*pz) - (2*mass*Vi_z(settings.Z1, ion.q));
            pzz = pp < 0 ? -sqrt(-pp) : sqrt(pp);
        }
        else
        {
            pp = pz * pz;
            pzz = pz;
        }
        //Ion is not escaping.
        if(pzz < 0)
        {
            theta = 0;
            phi = 90;
            E = -10;
        }
        else
        {
            //Recalculate E, incase image affected it
            E = 0.5 * psq / mass;
            theta = acos(pzz/sqrt(psq)) * 180 / M_PI;
            if(px == 0 && py==0)
            {
                phi = 90;
            }
            else
            {
                phi = atan2(py, px) * 180 / M_PI;
            }
        }
    }

    //Output data.
    sprintf(buffer, "%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%d\t%f\t%d\t%f\n",
            ion.r_0[0],ion.r_0[1],ion.r_0[2],
            E,theta,phi,
            1,1.0/settings.NUMCHA,
            ion.max_n, ion.r_min, ion.steps, ion.time);

    save(buffer);

    return;
}