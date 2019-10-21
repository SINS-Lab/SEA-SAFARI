#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <cstdio>
#include "safio.h"
#include "vec_math.h"
#include "space_math.h"
#include "traj.h"
#include <time.h>
#include <iomanip>
#include "scat.h"
#include "potentials.h"

//Initialize the global variables.
Safio settings;
std::ofstream out_file;
std::ofstream debug_file;
std::ofstream traj_file;
std::ofstream xyz_file;
std::ofstream crystal_file;

std::default_random_engine rng;
double space_mask[3375][3];

void test_cache()
{
    double v_tot = 0;
    double v_r_tot = 0;
    double r;
    double r_min = settings.DR_MIN_TAB*8.0/1e8;
    double dt_cache;
    double dt_compu;
    int n = 0;
    //Testing stuff goes here.

    clock_t timer = clock();
    //Testing Vr_r calls.
    for(int i = 1; i<=1e8; i++)
    {
        r = 1 + i * r_min;
        v_tot += Vr_r(r, 1);
        v_r_tot += dVr_dr(r, 1);
        n++;
    }
    dt_cache = ((double)clock() - timer) / CLOCKS_PER_SEC;    
    dt_cache /= n;
    //Convert to ms;
    dt_cache *= 1000;
    std::cout << v_tot << " " << v_r_tot << " ms per: " << dt_cache <<std::endl;


    v_tot = 0;
    v_r_tot = 0;
    n = 0;
    timer = clock();
    //Testing Vr_r_init calls.
    for(int i = 1; i<=1e8; i++)
    {
        r = 1 + i * r_min;
        v_tot += Vr_r_init(r, 1);
        v_r_tot += dVr_dr_init(r, 1);
        n++;
    }
    dt_compu = ((double)clock() - timer) / CLOCKS_PER_SEC;    
    dt_compu /= n;
    //Convert to ms;
    dt_compu *= 1000;
    std::cout << v_tot << " " << v_r_tot << " ms per: " << dt_compu <<std::endl;
    std::cout << "compu/cache: " << (dt_compu/dt_cache) <<std::endl;
}

int main(int argc,char* argv[])
{
    //Starts total runtime timer.
    clock_t load = clock();

    std::string safio_file = "safari.input";
    if(argc < 2)
    {
        std::ifstream input;
        input.open(safio_file.c_str());
        if (input.is_open())
        {
            while(getline(input,safio_file))
            {
                std::cout << "Loading Info From: " << safio_file + ".input" << '\n';
                break;
            }
            input.close();
        }
        else
        {
            std::cout << "Error opening Safari.input" << '\n';
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        safio_file = argv[1];
    }
    
    //Load the input file
    settings.load(safio_file);

    //Initialize potentials
    init_potentials();

    char buffer[200];

    Lattice lattice;
    lattice.build_lattice();

    int num = lattice.sites.size();

    for(int i = 0; i<num; i++)
    {
        Site s = *lattice.sites[i];
        Atom &a = s.atom;
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",s[0],s[1],s[2],a.charge,a.mass);
        crystal_file << buffer;
    }

    crystal_file.close();

    //Initialize the RNG
    rng.seed(settings.SEED);

    //Initialize the space_math's lookup table
    init_lookup();

    //Starts trajectory timer.
    clock_t start = clock();
    int n = 0;
    
    out_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tlevel\tweight\tmax_n\tmin_r\tsteps\ttotal time\n";

    if(settings.SCAT_FLAG == 666)
    {
        if(settings.NUMCHA==1 || settings.SCAT_TYPE == 777)
        {
            if(settings.NUMCHA==1)
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
        else if(settings.SCAT_TYPE == 666)
        {
            debug_file << "Running Montecarlo " << std::endl;
            std::cout << "Running Montecarlo " << std::endl;
            montecarloscat(lattice, &n);
        }
        else
        {
            chainscat(lattice, &n);
        }
        //Forces output of remaining trajectories.
        save(NULL);
    }
    else
    {
        test_cache();
    }
    
    
    //Compute time per trajectory.
    double dt = ((double)clock() - start) / CLOCKS_PER_SEC;
    dt /= n;
    //Convert to ms;
    dt *= 1000;

    std::cout << "\nFinished Running\n"<<std::endl;
    std::cout << "Time per particle: " << std::setprecision(2) << dt <<"ms"<<std::endl;
    debug_file << "\nTotal number particles: " << n <<std::endl;
    debug_file << "Time per particle: " << std::setprecision(2) << dt <<"ms"<<std::endl;
    //End final timer.
    dt = ((double)clock() - load) / CLOCKS_PER_SEC;
    std::cout << "Total Runtime: " << std::setprecision(2) << dt <<"s"<<std::endl;
    debug_file << "\nTotal Runtime: " << std::setprecision(2) << dt <<"s"<<std::endl;
    
    //Close files.
    out_file.close();
    debug_file.close();
    traj_file.close();
    return 0;
}
