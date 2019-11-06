#include <stdlib.h>     /* srand, rand */
#include <cstdio>
#include <time.h>       /* For runtime tracker */
#include <iomanip>

#include "safio.h"       /* settings */
#include "vec_math.h"    /* General math */
#include "space_math.h"  /* fast lookup stuff */
#include "string_utils.h"/* Used for argument parsing */
#include "scat.h"        /* The actual scattering routines */
#include "potentials.h"  /* We initialize these here */
#include "tests.h"       /* Debugging stuff */
#include "temps.h"       /* These are also initialized */

//Initialize the global variables.
Safio settings;
std::ofstream out_file;
std::ofstream debug_file;
std::ofstream traj_file;
std::ofstream xyz_file;
std::ofstream crystal_file;

double space_mask[3375][3];

int main(int argc, char* argv[])
{
    //Starts total runtime timer.
    double load = clock();

    std::map<std::string, ArgValue> args = get_arguments(argc, argv);
    std::string safio_file = args["-i"].as_string();

    args["-f"] = "t";

    //We had no input file specified via command line
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
        //Add the safio file to the arguments.
        args["-i"] = safio_file;
    }

    std::cout << "Loading Info From: " << safio_file + ".input" << '\n';
    //Load the input file
    settings.load(args);
    debug_file << "Loaded Settings, Initializing Potentials and Temperatures" << '\n';

    //Initialize potentials
    init_potentials();

    //Initialize temperatures
    init_temps();

    debug_file << "Initialized Potentials, Building Lattice" << '\n';

    char buffer[200];

    Lattice lattice;

    if(settings.load_crystal)
    {
        std::ifstream input;
        std::string crys_in_file = safio_file +".crys_in";
        debug_file << "Loading Lattice from "<< crys_in_file << '\n';
        input.open(crys_in_file);
        lattice.load_lattice(input);
        input.close();
        debug_file << "Lattice loaded"<<'\n';
        std::cout << "Lattice loaded" << '\n';
    }
    else
    {
        lattice.build_lattice();
    }
    
    std::ofstream crys_xyz_file;
    crys_xyz_file.open(settings.output_name+".crys.xyz");
    debug_file << "Printing Lattice" << '\n';
    std::cout << "Printing Lattice" << '\n';
    int num = lattice.sites.size();
    crys_xyz_file << num << "\n\n";
    for (int i = 0; i < num; i++)
    {
        Site &s = *lattice.sites[i];
        Atom* a = s.atom;
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",
                         s.r_0[0], s.r_0[1], s.r_0[2], a->charge, a->mass);
        crystal_file << buffer;
        sprintf(buffer, "%s\t%f\t%f\t%f\n",
                         a->symbol.c_str(), s.r_0[0], s.r_0[1], s.r_0[2]);
        crys_xyz_file << buffer;
    }
    crys_xyz_file.close();
    crystal_file.close();

    //Initialize the space_math's lookup table
    init_lookup();

    //Starts trajectory timer.
    double start = clock();
    int n = 0;

    if (settings.SCAT_FLAG == 666)
    {
        out_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tion index\tweight\tmax_n\tmin_r\tsteps\ttotal time\n";
        if (settings.NUMCHA == 1 || settings.SCAT_TYPE == 777)
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
        //Compute time per trajectory.
        double dt = (clock() - start) / CLOCKS_PER_SEC;
        dt /= n;
        //Convert to ms;
        dt *= 1000;

        std::cout << "\nFinished Running\n" << std::endl;
        std::cout << "Time per particle: " << std::setprecision(4) << dt << "ms" << std::endl;
        debug_file << "\nTotal number particles: " << n << std::endl;
        debug_file << "Time per particle: " << std::setprecision(4) << dt << "ms" << std::endl;
        //End final timer.
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
            test_lattice_copy(lattice);
        }
        else if (settings.SCAT_FLAG == 888)
        {
            test_rngs();
        }
        else if (settings.SCAT_FLAG == 999)
        {
            test_mask();
        }
        //Compute time per trajectory.
        double dt = (clock() - start) / CLOCKS_PER_SEC;
        std::cout << "Finished Running Tests, Time taken: " << dt << "s" << std::endl;
        debug_file << "Finished Running Tests, Time taken: " << dt << "s" << std::endl;
    }

    //Close files.
    out_file.close();
    debug_file.close();
    traj_file.close();
    return 0;
}
