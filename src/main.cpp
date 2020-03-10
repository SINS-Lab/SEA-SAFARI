#include <stdlib.h> /* srand, rand */
#include <cstdio>
#include <time.h> /* For runtime tracker */
#include <iomanip>

#include "safio.h"        /* settings */
#include "vec_math.h"     /* General math */
#include "space_math.h"   /* fast lookup stuff */
#include "string_utils.h" /* Used for argument parsing */
#include "scat.h"         /* The actual scattering routines */
#include "potentials.h"   /* We initialize these here */
#include "tests.h"        /* Debugging stuff */
#include "temps.h"        /* These are also initialized */
#include "safari.h"       /* This includes the exit_fail function*/

#define THREADCOUNT 10

// Initialize the global variables.
Safio settings;
std::ofstream out_file;
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
            std::cout << "Error opening Safari.input" << '\n';
            exit(EXIT_FAILURE);
        }
        // Add the safio file to the arguments.
        args["-i"] = safio_file;
    }

    std::cout << "Loading Info From: " << safio_file + ".input" << '\n';
    // Load the input file
    settings.load(args);
    debug_file << "Loaded Settings, Initializing Potentials and Temperatures" << '\n';

    // Initialize potentials
    init_potentials();

    // Initialize temperatures
    init_temps();

    debug_file << "Initialized Potentials, Building Lattice" << '\n';

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
        debug_file << "Lattice loaded" << '\n';
        std::cout << "Lattice loaded" << '\n';
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
        debug_file << settings.n_y_mask << std::endl;
        std::cout << settings.n_y_mask << std::endl;
        lattice.mask.num = settings.n_y_mask;
    }

    std::ofstream crys_xyz_file;
    crys_xyz_file.open(settings.output_name + ".crys.xyz");
    debug_file << "Printing Lattice" << '\n';
    std::cout << "Printing Lattice" << '\n';
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
        out_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tion index\tweight\tmax_n\tmin_r\tsteps\ttotal time\n";

        // 100 bifurcations is currently computationally infeasable, so
        // this will never be a valid choice for this, at least until
        // computers get many, many orders of magnitude better, when
        // that happens, here is where this needs to be changed!
        if (settings.SCAT_TYPE < 100 && settings.SCAT_TYPE)
        {
            // Otherwise this is an adaptive scat, with arguments of
            // the maximum depth to persue.
            debug_file << "Running Adaptive Grid " << std::endl;
            std::cout << "Running Adaptive Grid " << std::endl;

            default_detector.theta = settings.detect_parameters[0];
            default_detector.dtheta = settings.detect_parameters[1];
            default_detector.dphi = settings.detect_parameters[2];

            // We use numcha for the number of iterations of adaptive grid.
            // If temperature is 0, we only need 1 run.
            int runs = settings.TEMP > 0.0001 ? settings.NUMCHA : 1;

            if (runs == 1)
            {
                adaptivegridscat(settings.XSTART, settings.XSTEP, settings.XSTOP,
                                 settings.YSTART, settings.YSTEP, settings.YSTOP,
                                 lattice, default_detector, settings.SCAT_TYPE, 0, &n, 0);
            }
            else
            {
                #pragma omp parallel for num_threads(THREADCOUNT)
                for (int i = 0; i < runs; i++)
                {
                    //Cop the lattice
                    Lattice toUse = lattice;
                    toUse.clear_stats();
                    adaptivegridscat(settings.XSTART, settings.XSTEP, settings.XSTOP,
                                     settings.YSTART, settings.YSTEP, settings.YSTOP,
                                     toUse, default_detector, settings.SCAT_TYPE, 0, &n, i);
                    lattice.add_stats(toUse);
                }
            }
        }
        else if (settings.NUMCHA == 1 || settings.SCAT_TYPE == 777)
        {
            if (settings.NUMCHA == 1)
            {
                debug_file << "Running Single Shot" << std::endl;
                std::cout << "Running Single Shot" << std::endl;
            }
            else
            {
                debug_file << "Running Grid Scat" << std::endl;
                std::cout << "Running Grid Scat " << std::endl;
            }
            gridscat(lattice, &n);
        }
        else if (settings.SCAT_TYPE == 666)
        {
            debug_file << "Running Montecarlo " << std::endl;
            std::cout << "Running Montecarlo " << std::endl;
            montecarloscat(lattice, &n);
        }
        else if (settings.SCAT_TYPE == 888)
        {
            chainscat(lattice, &n);
        }

        // Log some debug info from the lattice
        debug_file << "Total Out of Phi(-5): " << lattice.undetectable_num << std::endl;
        debug_file << "Total Trapped  (-10): " << lattice.trapped_num << std::endl;
        debug_file << "Total Stuck   (-100): " << lattice.stuck_num << std::endl;
        debug_file << "Total Buried  (-200): " << lattice.buried_num << std::endl;
        debug_file << "Total Froze   (-300): " << lattice.froze_num << std::endl;
        debug_file << "Total OOB     (-400): " << lattice.left_num << std::endl;
        debug_file << "Total Errored (-500): " << lattice.err_num << std::endl;
        debug_file << "Total Out of Mask: " << lattice.out_of_mask << std::endl;

        // Compute time per trajectory.
        double dt = (clock() - start) / CLOCKS_PER_SEC;
        dt /= n;
        // Convert to ms;
        dt *= 1000;

        std::cout << "\nFinished Running " << safio_file << "\n"
                  << std::endl;
        std::cout << "Time per particle: " << std::setprecision(4) << dt << "ms" << std::endl;
        debug_file << "\nTotal number particles: " << n << std::endl;
        debug_file << "Time per particle: " << std::setprecision(4) << dt << "ms" << std::endl;
        // End final timer.
        dt = (clock() - load) / CLOCKS_PER_SEC;
        std::cout << "Total Runtime: " << std::setprecision(4) << dt << "s" << std::endl;
        debug_file << "\nTotal Runtime: " << std::setprecision(4) << dt << "s" << std::endl;
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
        std::cout << "Finished Running Tests, Time taken: " << dt << "s" << std::endl;
        debug_file << "Finished Running Tests, Time taken: " << dt << "s" << std::endl;
    }

    // Close files.
    out_file.close();
    debug_file.close();
    traj_file.close();
    return 0;
}

void exit_fail(std::string reason)
{
    std::cout << "Exiting Early, for reason:" << std::endl;
    std::cout << reason << std::endl;
    debug_file << "Exiting Early, for reason:" << std::endl;
    debug_file << reason << std::endl;

    // Close files.
    out_file.close();
    debug_file.close();
    traj_file.close();

    // Exit
    exit(EXIT_FAILURE);
}