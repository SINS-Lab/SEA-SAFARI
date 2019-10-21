#include "scat.h"
#include "safio.h"
#include "traj.h"

/**
 * Returns a random, floating point
 * number from 0 to 1;
 */ 
double frand()
{
    double max = rng.max();
    return rng()/max;
}

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
    ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x, y);
    ion.index = index;
    traj(ion, lattice, log, xyz);
}

void montecarloscat(Lattice &lattice, int *num)
{
    //Before montecarloscat is called, the RNG has its
    //Seed set, this will result in each seperate call
    //being identical, unless the seed is changed.

    //Find the range of the region we are covering
    double x_size = (settings.XSTOP - settings.XSTART);
    double y_size = (settings.YSTOP - settings.YSTART);

    for(int n = 0; n < settings.NUMCHA; n++)
    {
        //Find a new random location, these rands are 0-1
        double rx = frand();
        double ry = frand();
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
    bool log = settings.NUMCHA == 1;

    int n = 0;
    if(log)
    {
        //If we are in log mode, we only run this for the first coord
        Ion ion;
        fire(lattice, ion, settings.XSTART, settings.YSTART, n, log, log);
        n++;
    }
    else
    {
        for(double x = settings.XSTART; x <= settings.XSTOP; x+=settings.XSTEP)
            for(double y = settings.YSTART; y <= settings.YSTOP; y+=settings.YSTEP)
            {
                Ion ion;
                fire(lattice, ion, x, y, n, log, log);
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
