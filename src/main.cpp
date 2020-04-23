#include <stdlib.h> /* srand, rand */
#include <cstdio>   /* the output file streams*/
#include <time.h>   /* for runtime tracker */
#include <iomanip>  /* precision setting on output numbers */

#include "safio.h"        /* settings */
#include "vec_math.h"     /* General math */
#include "space_math.h"   /* fast lookup stuff */
#include "string_utils.h" /* Used for argument parsing */
#include "scat.h"         /* The actual scattering routines */
#include "potentials.h"   /* We initialize these here */
#include "tests.h"        /* Debugging stuff */
#include "temps.h"        /* These are also initialized */
#include "safari.h"       /* This includes the exit_fail function*/

// Initialize the global variables.
Safio settings;
std::ofstream out_file;
std::ofstream sptr_file;
std::ofstream debug_file;
std::ofstream traj_file;
std::ofstream xyz_file;
std::ofstream crystal_file;
Detector default_detector;

double space_mask[3375][3];

int main(int argc, char *argv[])
{
    // Starts total runtime timer.
    double load = clock();

    std::map<std::string, ArgValue> args = get_arguments(argc, argv);
    std::string safio_file = args["-i"].as_string();

    args["-f"] = "t";

    // We had no input file specified via command line
    if (!args["-i"])
    {
        std::ifstream input;
        input.open("safari.input");
        if (input.is_open())
        {
            while (getline(input, safio_file))
            {
                break;
            }
            input.close();
        }
        else
        {
            std::cout << "Error opening Safari.input\n";
            exit(EXIT_FAILURE);
        }
        // Add the safio file to the arguments.
        args["-i"] = safio_file;
    }

    std::cout << "Loading Info From: " << safio_file + ".input\n";
    // Load the input file
    settings.load(args);
    debug_file << "Loaded Settings, Initializing Potentials and Temperatures\n";

    // Initialize potentials
    init_potentials();

    // Initialize temperatures
    init_temps();

    debug_file << "Initialized Potentials, Building Lattice\n";

    char buffer[200];

    Lattice lattice;

    if (settings.load_crystal)
    {
        std::ifstream input;
        std::string crys_in_file = safio_file + ".crys_in";
        debug_file << "Loading Lattice from " << crys_in_file << '\n';
        input.open(crys_in_file);
        lattice.load_lattice(input);
        input.close();
        debug_file << "Lattice loaded\n";
        std::cout << "Lattice loaded\n";
    }
    else
    {
        lattice.build_lattice();
    }

    // Initialize springs if not using einstein
    if (!settings.useEinsteinSprings)
    {
        lattice.init_springs(settings.neighbour_count);
    }

    // Only update the mask if it has more than 2 points, and has same for x and y.
    if (settings.n_x_mask == settings.n_y_mask && settings.n_y_mask > 2)
    {
        debug_file << "Initializing surface mask of size ";
        std::cout << "Initializing surface mask of size ";

        for (int i = 0; i < settings.n_y_mask; i++)
        {
            Point point;
            point.x = settings.x_mask_points[i];
            point.y = settings.y_mask_points[i];
            lattice.mask.points[i] = point;
        }
        debug_file << settings.n_y_mask << "\n";
        std::cout << settings.n_y_mask << "\n";
        lattice.mask.num = settings.n_y_mask;
    }

    std::ofstream crys_xyz_file;
    crys_xyz_file.open(settings.output_name + ".crys.xyz");
    debug_file << "Printing Lattice\n";
    std::cout << "Printing Lattice\n"
              << std::flush;
    int num = lattice.sites.size();
    crys_xyz_file << num << "\n\n";
    for (int i = 0; i < num; i++)
    {
        Site &s = *lattice.sites[i];
        Atom *a = s.atom;
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",
                s.r_0[0], s.r_0[1], s.r_0[2], a->charge, a->mass);
        crystal_file << buffer;
        sprintf(buffer, "%s\t%f\t%f\t%f\n",
                a->symbol.c_str(), s.r_0[0], s.r_0[1], s.r_0[2]);
        crys_xyz_file << buffer;
    }
    crys_xyz_file.close();
    crystal_file.close();

    // Initialize the space_math's lookup table
    init_lookup();

    // Starts trajectory timer.
    double start = clock();
    int n = 0;

    default_detector.phi = settings.PHI0;
    default_detector.theta = 45;
    default_detector.e_min = settings.EMIN;

    if (settings.SCAT_FLAG == 666)
    {
        settings.scat_started = true;
        out_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tion index\tweight\tmax_n\tmin_r\tsteps\tMax Error\ttotal time" << std::endl;
        if (settings.saveSputter)
            sptr_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tion index\tweight\tmax_n\tmin_r\tsteps\tMax Error\ttotal time" << std::endl;

        if (settings.singleshot)
        {
            debug_file << "Running Single Shot"
                       << "\n";
            std::cout << "Running Single Shot"
                      << "\n"
                      << std::flush;
            singleshot(&lattice, &n);
        }
        else if (settings.adaptivegrid)
        {
            // Otherwise this is an adaptive scat, with arguments of
            // the maximum depth to persue.
            debug_file << "Running Adaptive Grid\n";
            std::cout << "Running Adaptive Grid\n"
                      << std::flush;

            default_detector.theta = settings.detect_parameters[0];
            default_detector.dtheta = settings.detect_parameters[1];
            default_detector.dphi = settings.detect_parameters[2];

            // We use numcha for the number of iterations of adaptive grid.
            // If temperature is 0, we only need 1 run.
            int runs = settings.TEMP > 0.0001 ? settings.NUMCHA : 1;
            int ind = 0;

            if (runs == 1)
            {
                adaptivegridscat(settings.XSTART, settings.XSTEP, settings.XSTOP,
                                 settings.YSTART, settings.YSTEP, settings.YSTOP,
                                 &lattice, default_detector, settings.SCAT_TYPE, 0, &n, &ind, 0);
            }
            else
            {
                #pragma omp parallel for num_threads(THREADCOUNT)
                for (int i = 0; i < runs; i++)
                {
                    //Copy the lattice
                    Lattice *toUse = new Lattice(lattice);
                    toUse->clear_stats();
                    adaptivegridscat(settings.XSTART, settings.XSTEP, settings.XSTOP,
                                     settings.YSTART, settings.YSTEP, settings.YSTOP,
                                     toUse, default_detector, settings.SCAT_TYPE, 0, &n, &ind, i);
                    lattice.add_stats(toUse);
                    delete toUse;
                }
            }
        }
        else if (settings.gridscat)
        {
            debug_file << "Running Grid Scat"
                       << "\n";
            std::cout << "Running Grid Scat "
                      << "\n"
                      << std::flush;
            gridscat(&lattice, &n);
        }
        else if (settings.montecarlo)
        {
            debug_file << "Running Montecarlo "
                       << "\n";
            std::cout << "Running Montecarlo "
                      << "\n"
                      << std::flush;

            double seeds[THREADCOUNT];
            std::default_random_engine rng;
            //Initialize the RNG
            rng.seed(make_seed(settings.SEED));
            for (int i = 0; i < THREADCOUNT; i++)
            {
                seeds[i] = frand(rng);
            }
            int ions_per_thread = settings.NUMCHA / THREADCOUNT;

            #pragma omp parallel for num_threads(THREADCOUNT)
            for (int i = 0; i < THREADCOUNT; i++)
            {
                int start = i * ions_per_thread;
                //Copy the lattice
                Lattice *toUse = new Lattice(lattice);
                std::cout << "Starting Thread " << i << "\n"
                          << std::flush;
                toUse->clear_stats();
                montecarloscat(toUse, start, ions_per_thread, seeds[i]);
                lattice.add_stats(toUse);
                std::cout << "Finished Thread " << i << "\n"
                          << std::flush;
                delete toUse;
            }
            n = ions_per_thread * THREADCOUNT;
        }
        else if (settings.chainscat)
        {
            debug_file << "Running Chainscat"
                       << "\n";
            std::cout << "Running Chainscat"
                      << "\n"
                      << std::flush;
            chainscat(&lattice, &n);
        }

        // Log some debug info from the lattice
        debug_file << "Total Out of Phi  (-5): " << lattice.undetectable_num << "\n";
        debug_file << "Total Trapped     (-10): " << lattice.trapped_num << "\n";
        debug_file << "Total Stuck       (-100): " << lattice.stuck_num << "\n";
        debug_file << "Total Buried      (-200): " << lattice.buried_num << "\n";
        debug_file << "Total Froze       (-300): " << lattice.froze_num << "\n";
        debug_file << "Total OOB         (-400): " << lattice.left_num << "\n";
        debug_file << "Total Errored     (-500): " << lattice.err_num << "\n";
        debug_file << "Total Intersected (-600): " << lattice.intersections << "\n";
        debug_file << "Total Out of Mask: " << lattice.out_of_mask << "\n";

        if (settings.dynamicNeighbours)
        {
            debug_file << "\nMax active sites: " << lattice.max_active << "\n";
            debug_file << "Mean active sites: " << (lattice.sum_active / lattice.count_active) << "\n\n";
        }

        // Compute time per trajectory.
        double dt = (clock() - start) / CLOCKS_PER_SEC;
        dt /= n;
        // Convert to ms;
        dt *= 1000;

        std::cout << "\nFinished Running " << safio_file << "\n"
                  << "\n";
        std::cout << "Time per particle: " << std::setprecision(4) << dt << "ms\n";
        debug_file << "\nTotal number particles: " << n << "\n";
        debug_file << "Time per particle: " << std::setprecision(4) << dt << "ms\n";
        // End final timer.
        dt = (clock() - load) / CLOCKS_PER_SEC;
        std::cout << "Total Runtime: " << std::setprecision(4) << dt << "s\n";
        debug_file << "\nTotal Runtime: " << std::setprecision(4) << dt << "s\n";
    }
    else
    {
        n = 1;
        if (settings.SCAT_FLAG == 555)
        {
            test_cache();
        }
        else if (settings.SCAT_FLAG == 777)
        {
            switch (settings.SCAT_TYPE)
            {
            case 666:
                test_lattice_copy(lattice);
                break;
            case 777:
                test_lattice_springs(lattice);
                break;
            default:
                break;
            }
        }
        else if (settings.SCAT_FLAG == 888)
        {
            test_rngs();
        }
        else if (settings.SCAT_FLAG == 999)
        {
            test_mask(lattice);
        }
        // Compute time per trajectory.
        double dt = (clock() - start) / CLOCKS_PER_SEC;
        std::cout << "Finished Running Tests, Time taken: " << dt << "s\n";
        debug_file << "Finished Running Tests, Time taken: " << dt << "s\n";
    }

    // Close files.
    out_file.close();
    debug_file.close();
    traj_file.close();
    sptr_file.close();
    return 0;
}

void exit_fail(std::string reason)
{
    std::cout << "Exiting Early, for reason:\n";
    std::cout << reason << "\n";
    debug_file << "Exiting Early, for reason:\n";
    debug_file << reason << "\n";

    // Close files.
    out_file.close();
    debug_file.close();
    traj_file.close();
    sptr_file.close();

    // Exit
    exit(EXIT_FAILURE);
}
