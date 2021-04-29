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
std::mutex mutx;

double space_mask_cube[N_CUBE_MASK][3];

int main(int argc, char *argv[])
{
    std::cout << "Starting SAFARI!\n" << std::flush;
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
    std::cout << "Loaded Settings, Initializing Potentials and Temperatures\n";

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

    if (settings.SCAT_FLAG == 666)
    {
        settings.scat_started = true;
        if(out_file.is_open())
            out_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tion index\tweight\tmax_n\tmin_r\tsteps\tMax Error\ttotal time" << std::endl;
        if (settings.saveSputter)
            sptr_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tion index\tion err flag\tmax_n\tmin_r\tsteps\tMax Error\ttotal time" << std::endl;

        ScatRoutine* toRun = NULL;

        if (settings.singleshot)
        {
            toRun = new SingleShot();
        }
        else if (settings.adaptivegrid)
        {
            toRun = new AdaptiveGrid();
        }
        else if (settings.gridscat)
        {
            toRun = new GridScat();
        }
        else if (settings.montecarlo)
        {
            toRun = new MonteCarlo();
        }
        else if (settings.chainscat)
        {
            toRun = new ChainScat();
        }

        if(toRun != NULL)
        {
            toRun->run(&lattice, &n);
        }
        else
        {
            exit_fail("Unkown Scat Mode??\n");
        }

        // Log some debug info from the lattice
        const char* debugInfo = "\n"
        "Total Out of Phi  (-5):   %d\n"
        "Total Trapped     (-10):  %d\n"
        "Total Stuck       (-100): %d\n"
        "Total Buried      (-200): %d\n"
        "Total Froze       (-300): %d\n"
        "Total OOB         (-400): %d\n"
        "Total Errored     (-500): %d\n"
        "Total Intersected (-600): %d\n"
        "Total Out of Mask: %d\n";

        char buffer[1024];
        sprintf(buffer, debugInfo, lattice.undetectable_num, lattice.trapped_num, 
               lattice.stuck_num, lattice.buried_num, lattice.froze_num, lattice.left_num, 
               lattice.err_num, lattice.intersections, lattice.out_of_mask);

        debug_file << buffer;

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

        std::cout << "\nFinished Running " << safio_file << "\n";
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
    if(out_file.is_open()) out_file.close();
    if(debug_file.is_open()) debug_file.close();
    if(traj_file.is_open()) traj_file.close();
    if(sptr_file.is_open()) sptr_file.close();
    return 0;
}

void exit_fail(std::string reason)
{
    std::cout << "Exiting Early, for reason:\n";
    std::cout << reason << "\n";
    debug_file << "Exiting Early, for reason:\n";
    debug_file << reason << "\n";

    // Close files.
    if(out_file.is_open()) out_file.close();
    if(debug_file.is_open()) debug_file.close();
    if(traj_file.is_open()) traj_file.close();
    if(sptr_file.is_open()) sptr_file.close();

    // Exit
    exit(EXIT_FAILURE);
}