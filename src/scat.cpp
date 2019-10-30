#include "scat.h"
#include "safio.h"
#include "traj.h"
#include "temps.h"

/**
 * Fires the given ion at the given lattice.
 * This also initializes the ions kinetic energy, and sets the
 * impact parameters to the given coordinates. 
 * The ion is assigned an index as well.
 * 
 * @param lattice - The lattice to shoot at
 * @param ion - The ion to shoot at the lattice
 * @param x - the x-impact parameter for the ion
 * @param y - the y-impact parameter for the ion
 * @param index - the index for this ion
 * @param log - whether to log the trajectory
 * @param xyz - whether to log xyz file version of the trajectory.
 * 
 */ 
void fire(Lattice &lattice, Ion &ion, double x, double y, int index, bool log, bool xyz)
{
    ion.E0 = settings.E0;
    if(settings.ESIZE > 0)
    {
        thermaize_ion(ion);
    }
    ion.set_KE(settings.THETA0, settings.PHI0, x, y);
    ion.index = index;
    traj(ion, lattice, log, xyz);
}

void montecarloscat(Lattice &lattice, int *num)
{
    //Make a new RNG instance, and then set the seed to what it should be
    std::default_random_engine rng;
    debug_file << "Initializing RNG, seed: " << settings.SEED << '\n';
    //Initialize the RNG
    rng.seed(make_seed(settings.SEED));

    //Find the range of the region we are covering
    double x_size = (settings.XSTOP - settings.XSTART);
    double y_size = (settings.YSTOP - settings.YSTART);

    for(int n = 0; n < settings.NUMCHA; n++)
    {
        //Find a new random location, these rands are 0-1
        double rx = frand(rng);
        double ry = frand(rng);
        //Use the random numbers to find the impact parameter
        double x = settings.XSTART + x_size * rx;
        double y = settings.YSTART + y_size * ry;

        //Fire the ion at the random spot.
        Ion ion;
        fire(lattice, ion, x, y, n, false, false);
    }
    *num = settings.NUMCHA;
    return;
}

void gridscat(Lattice &lattice, int *num)
{
    int n = 0;
    if(settings.NUMCHA == 1)
    {
        //If we are in log mode, we only run this for the first coord
        Ion ion;
        //We also use the settings ion_index for single shot mode, to
        //guarentee that any thermal effects can be repeated.
        fire(lattice, ion, settings.XSTART, settings.YSTART, settings.ion_index, true, true);
        n++;
    }
    else
    {
        for(double x = settings.XSTART; x <= settings.XSTOP; x+=settings.XSTEP)
            for(double y = settings.YSTART; y <= settings.YSTOP; y+=settings.YSTEP)
            {
                Ion ion;
                fire(lattice, ion, x, y, n, false, false);
                n++;
            }
    }
    
    *num = n;
    return;
}

void chainscat(Lattice &lattice, int *num)
{
    double dx = settings.XSTART - settings.XSTOP;
    double dy = settings.YSTART - settings.YSTOP;

    //Computes the step sizes in x and y
    dx /= settings.NUMCHA;
    dy /= settings.NUMCHA;

    double x, y;
    int n = 0;
    for(n = 0; n < settings.NUMCHA; n++)
    {
        //Step in the line defined by dx, dy.
        x = settings.XSTART + dx * n;
        y = settings.YSTART + dy * n;
        Ion ion;
        //Fire at next point in the line.
        fire(lattice, ion, x, y, n, false, false);
        n++;
    }
    *num = n;
    return;
}
