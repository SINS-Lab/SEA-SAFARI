#include "scat.h"
#include "safio.h"
#include "traj.h"
#include "temps.h"

/**
 * Fires the given ion at the given lattice->
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
bool fire(Lattice *lattice, Ion &ion, double x, double y, int index, bool log, bool xyz)
{
    if (!lattice->mask.inside(x, y))
    {
        lattice->out_of_mask++;
        return false;
    }
    ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x, y);
    ion.index = index;
    ion.thermal_seed = index;
    traj(ion, lattice, log, xyz, default_detector);
    return true;
}

void montecarloscat(Lattice *lattice, int ionStart, int numcha, double seed)
{
    //Make a new RNG instance, and then set the seed to what it should be
    std::default_random_engine rng;
    debug_file << "Initializing RNG, seed: " << seed << '\n';
    //Initialize the RNG
    rng.seed(make_seed(seed));

    //Find the range of the region we are covering
    double x_size = (settings.XSTOP - settings.XSTART);
    double y_size = (settings.YSTOP - settings.YSTART);

    for (int n = 0; n < numcha; n++)
    {
        //Find a new random location, these rands are 0-1
        double rx = frand(rng);
        double ry = frand(rng);
        //Use the random numbers to find the impact parameter
        double x = settings.XSTART + x_size * rx;
        double y = settings.YSTART + y_size * ry;

        //Fire the ion at the random spot.
        Ion ion;
        if (!fire(lattice, ion, x, y, n + ionStart, false, false))
            n--;
    }
}

void gridscat(Lattice *lattice, int *num)
{
    int n = 0;
    if (settings.NUMCHA == 1)
    {
        //If we are in log mode, we only run this for the first coord
        Ion ion;
        //We also use the settings ion_index for single shot mode, to
        //guarentee that any thermal effects can be repeated.
        if (settings.SCAT_TYPE)
        {
            Lattice *initial = new Lattice(*lattice);
            std::cout << "Single Shot Partial Lattice Mode\n";
            //Dry run to find bounds of operation
            fire(initial, ion, settings.XSTART, settings.YSTART, settings.ion_index, true, false);
            delete initial;

            //Set min bounds from ion positions
            lattice->xyz_bounds[0] = std::min(ion.r_0[0], ion.r[0]) - settings.AX;
            lattice->xyz_bounds[1] = std::min(ion.r_0[1], ion.r[1]) - settings.AY;
            lattice->xyz_bounds[2] = std::min(ion.r_0[2], ion.r[2]) - 2 * settings.AZ;

            //Set max bounds from ion positions
            lattice->xyz_bounds[3] = std::max(ion.r_0[0], ion.r[0]) + settings.AX;
            lattice->xyz_bounds[4] = std::max(ion.r_0[1], ion.r[1]) + settings.AY;
            lattice->xyz_bounds[5] = std::max(ion.r_0[2], ion.r[2]);

            //Lets make this a square slab, so find midpoint of x,y
            //and then make it equal sized.
            double dx = lattice->xyz_bounds[3] - lattice->xyz_bounds[0];
            double dy = lattice->xyz_bounds[4] - lattice->xyz_bounds[1];
            double dr = std::max(dx, dy) / 2 + 5; //Add extra 5 for some padding
            //Set x/y to avg +- dr
            lattice->xyz_bounds[0] = (lattice->xyz_bounds[3] + lattice->xyz_bounds[0]) / 2 - dr;
            lattice->xyz_bounds[3] = (lattice->xyz_bounds[3] + lattice->xyz_bounds[0]) / 2 + dr;
            lattice->xyz_bounds[1] = (lattice->xyz_bounds[4] + lattice->xyz_bounds[1]) / 2 - dr;
            lattice->xyz_bounds[4] = (lattice->xyz_bounds[4] + lattice->xyz_bounds[1]) / 2 + dr;

            std::cout << "Single Shot Pass 2\n";

            //Run to output the xyz file
            fire(lattice, ion, settings.XSTART, settings.YSTART, settings.ion_index, false, true);
        }
        else
        {
            std::cout << "Single Shot Full Lattice Mode\n";
            //Fire, log entire lattice, and also output xyz at once
            fire(lattice, ion, settings.XSTART, settings.YSTART, settings.ion_index, true, true);
        }
        if (ion.r[2] < settings.Z1)
        {
            std::cout << "Ion buried!\n";
        }
        else
            std::cout << "Ion Escaped!\n";
        n++;
    }
    else
    {
        for (double y = settings.YSTART; y <= settings.YSTOP; y += settings.YSTEP)
            for (double x = settings.XSTART; x <= settings.XSTOP; x += settings.XSTEP)
            {
                Ion ion;
                if (fire(lattice, ion, x, y, n, false, false))
                    n++;
            }
    }

    *num = n;
    return;
}

void chainscat(Lattice *lattice, int *num)
{
    double dx = settings.XSTART - settings.XSTOP;
    double dy = settings.YSTART - settings.YSTOP;

    //Computes the step sizes in x and y
    dx /= settings.NUMCHA;
    dy /= settings.NUMCHA;

    double x, y;
    int n = 0;
    int steps = settings.NUMCHA / 100;
    steps = std::max(1, steps);
    for (n = 0; n < settings.NUMCHA; n++)
    {
        //Step in the line defined by dx, dy.
        x = settings.XSTART + dx * n;
        y = settings.YSTART + dy * n;
        Ion ion;
        //Fire at next point in the line.
        fire(lattice, ion, x, y, n, false, false);
        if ((n % steps) == steps - 1)
        {
            std::cout << "\x1B[32mx\x1B[0m" << std::flush;
        }
    }
    std::cout << "\n";
    *num = n;
    return;
}

void save_adaptive(Ion &ion, Lattice *lattice, double E, double theta, double phi)
{
    if (settings.save_errored || E > 0)
    {
        char buffer[200];
        //first stuff it in the buffer
        sprintf(buffer, "%f\t%f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\n",
                ion.r_0[0], ion.r_0[1], ion.r_0[2],
                E, theta, phi,
                ion.index, 1.0,
                ion.max_n, ion.r_min, ion.steps, ion.time);
        //Then save it
        out_file << buffer;
    }
}

void adaptivegridscat(double xstart, double xstep, double xstop,
                      double ystart, double ystep, double ystop,
                      Lattice *lattice, Detector &detector,
                      int max_depth, int current_depth, int *num, int *index, int iter)
{
    if (current_depth > max_depth)
        return;
    bool print = iter == 0 and current_depth == 0;
    bool log = false;
    bool xyz = false;
    int d = current_depth + 1;
    for (double y = ystart; y <= ystop; y += ystep)
        for (double x = xstart; x <= xstop; x += xstep)
        {
            if (lattice->mask.inside(x, y))
            {
                Ion ion;
                ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x, y);
                ion.index = (*index) + 1;
                *index = ion.index;
                ion.weight = current_depth;
                ion.thermal_seed = iter;
                traj(ion, lattice, log, xyz, detector);
                *num = *num + 1;
                if (ion.index)
                {
                    if (print)
                    {
                        if ((x + xstep) > xstop)
                            std::cout << "\x1B[32mx\x1B[0m\n"
                                      << std::flush;
                        else
                            std::cout << "\x1B[32mx\x1B[0m" << std::flush;
                    }
                    //Do a higher resolution scan around the ion.
                    double dx = xstep / 2;
                    double dy = ystep / 2;
                    double x_min = x - dx / 2;
                    double y_min = y - dy / 2;
                    double x_max = x + dx;
                    double y_max = y + dy;
                    adaptivegridscat(x_min, dx, x_max,
                                     y_min, dy, y_max,
                                     lattice, detector,
                                     max_depth, d, num, index, iter);
                }
                else if (print)
                {
                    if ((x + xstep) > xstop)
                        std::cout << "\x1B[31mo\x1B[0m\n"
                                  << std::flush;
                    else
                        std::cout << "\x1B[31mo\x1B[0m" << std::flush;
                }
            }
            else
            {
                lattice->out_of_mask++;
                if (print)
                {
                    if ((x + xstep) > xstop)
                        std::cout << "\x1B[30m.\x1B[0m\n"
                                  << std::flush;
                    else
                        std::cout << "\x1B[30m.\x1B[0m" << std::flush;
                }
            }
        }
}