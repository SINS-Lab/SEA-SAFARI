#include "scat.h"
#include "safio.h"
#include "traj.h"
#include "temps.h"
#include "safari.h"

void MonteCarlo::run(Lattice *lattice, int *num)
{
    debug_file << "Running Montecarlo\n";
    std::cout << "Running Montecarlo\n"<< std::flush;

    double seeds[THREADCOUNT];
    std::default_random_engine rng;
    //Initialize the RNG
    rng.seed(make_seed(settings.SEED));
    for (int i = 0; i < THREADCOUNT; i++)
    {
        seeds[i] = frand(rng);
    }
    int ions_per_thread = settings.NUMCHA / THREADCOUNT;

    setupDetector();

    #pragma omp parallel for num_threads(THREADCOUNT)
    for (int i = 0; i < THREADCOUNT; i++)
    {
        int start = i * ions_per_thread;
        //Copy the lattice
        Lattice *toUse = new Lattice(*lattice);
        mutx.lock();
        std::cout << "Starting Thread " << i << "\n"
                    << std::flush;
        mutx.unlock();
        toUse->clear_stats();
        montecarloscat(toUse, start, ions_per_thread, seeds[i]);
        lattice->add_stats(toUse);
        mutx.lock();
        std::cout << "Finished Thread " << i << "\n"
                    << std::flush;
        mutx.unlock();
        delete toUse;
    }
    *num = ions_per_thread * THREADCOUNT;

    cleanUpDetector();
}

void MonteCarlo::montecarloscat(Lattice *lattice, int ionStart, int numcha, double seed)
{
    //Make a new RNG instance, and then set the seed to what it should be
    mutx.lock();
    std::default_random_engine rng;
    debug_file << "Initializing RNG, seed: " << seed << '\n';
    mutx.unlock();
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

void GridScat::run(Lattice *lattice, int *num)
{
    debug_file << "Running Grid Scat\n";
    std::cout << "Running Grid Scat\n" << std::flush;

    setupDetector();
    
    int n = 0;
	for (double y = settings.YSTART; y <= settings.YSTOP; y += settings.YSTEP)
		for (double x = settings.XSTART; x <= settings.XSTOP; x += settings.XSTEP)
		{
			Ion ion;
			if (fire(lattice, ion, x, y, n, false, false))
				n++;
		}
    *num = n;

    cleanUpDetector();
}

void SingleShot::run(Lattice *lattice, int *num)
{
    debug_file << "Running Single Shot\n";
    std::cout << "Running Single Shot\n" << std::flush;

    setupDetector();

    // We want to print out a copy of the lattice with the given thermalisation, if T > 0.
    if (settings.TEMP > 0)
    {
        std::ofstream crys_xyz_file;
        std::ofstream debug_file;
        crys_xyz_file.open(settings.output_name + "_thermalised.crys.xyz");
        crystal_file.open(settings.output_name + "_thermalised.crys");
        debug_file << "Printing Lattice\n";
        std::cout << "Printing Lattice\n"
                << std::flush;
        int num = lattice->sites.size();
        crys_xyz_file << num << "\n\n";
        char buffer[200];
        for (int i = 0; i < num; i++)
        {
            Site &s = *(lattice->sites[i]);
            s.thermal_seed = settings.ion_index;
            s.reset();
            Atom *a = s.atom;
            sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",
                    s.r[0], s.r[1], s.r[2], a->charge, a->mass);
            crystal_file << buffer;
            sprintf(buffer, "%s\t%f\t%f\t%f\n",
                    a->symbol.c_str(), s.r[0], s.r[1], s.r[2]);
            crys_xyz_file << buffer;
        }
        crys_xyz_file.close();
        crystal_file.close();
    }

    //If we are in log mode, we only run this for the first coord
    Ion ion;
    //We also use the settings ion_index for single shot mode, to
    //guarentee that any thermal effects can be repeated.
    if (settings.SCAT_TYPE)
    {
        Lattice *initial = new Lattice(*lattice);
        std::cout << "Single Shot Partial Lattice Mode\n" << std::flush;
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

        std::cout << "Single Shot Pass 2\n" << std::flush;

        //Run to output the xyz file
        fire(lattice, ion, settings.XSTART, settings.YSTART, settings.ion_index, false, true);
        *num = 2;
    }
    else
    {
        std::cout << "Single Shot Full Lattice Mode\n" << std::flush;
        //Fire, log entire lattice, and also output xyz at once
        fire(lattice, ion, settings.XSTART, settings.YSTART, settings.ion_index, true, true);
        *num = 1;
    }
    if (ion.r[2] < settings.Z1)
    {
        std::cout << "Ion buried!\n" << std::flush;
    }
    else
        std::cout << "Ion Escaped!\n" << std::flush;

    cleanUpDetector();
}

void ChainScat::run(Lattice *lattice, int *num)
{            
    debug_file << "Running Chainscat\n";
    std::cout << "Running Chainscat\n" << std::flush;

    setupDetector();

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

    cleanUpDetector();
}

void AdaptiveGrid::run(Lattice *lattice, int *num)
{
    // Otherwise this is an adaptive scat, with arguments of
    // the maximum depth to persue.
    debug_file << "Running Adaptive Grid\n";
    std::cout << "Running Adaptive Grid\n" << std::flush;

    // We use numcha for the number of iterations of adaptive grid.
    // If temperature is 0, we only need 1 run.
    int runs = settings.TEMP > 0.0001 ? settings.NUMCHA : 1;
    int ind = 0;

    setupDetector();

    if (runs == 1)
    {
        adaptivegridscat(settings.XSTART, settings.XSTEP, settings.XSTOP,
                            settings.YSTART, settings.YSTEP, settings.YSTOP,
                            lattice, settings.SCAT_TYPE, 0, num, &ind, 0);
    }
    else
    {
        #pragma omp parallel for num_threads(THREADCOUNT)
        for (int i = 0; i < runs; i++)
        {
            //Copy the lattice
            Lattice *toUse = new Lattice(*lattice);
            toUse->clear_stats();
            adaptivegridscat(settings.XSTART, settings.XSTEP, settings.XSTOP,
                                settings.YSTART, settings.YSTEP, settings.YSTOP,
                                toUse, settings.SCAT_TYPE, 0, num, &ind, i);
            lattice->add_stats(toUse);
            delete toUse;
        }
    }

    cleanUpDetector();
}

void AdaptiveGrid::adaptivegridscat(double xstart, double xstep, double xstop,
                      double ystart, double ystep, double ystop,
                      Lattice *lattice,
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
                traj(ion, lattice, log, xyz, spot_detector, area_detector);
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
                                     lattice,
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
