#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <cstdio>
#include "safio.h"
#include "vec_math.h"
#include "space_math.h"
#include "ion.h"
#include <time.h>
#include "scat.h"

//Initialize the global variables.
Safio settings;
std::ofstream out_file;
std::ofstream debug_file;
std::ofstream traj_file;
std::default_random_engine rng;
std::ofstream crystal_file;

int main()
{
    //Load the input file
    settings.load();

    crystal_file.open("crystal.input");

    char buffer[200];

    Lattice lattice;
    lattice.build_lattice();

    int num = lattice.sites.size();

    for(int i = 0; i<num; i++)
    {
        Site s = lattice.sites[i];
        Atom a = s.atom;
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",s[0],s[1],s[2],a.charge,a.mass);
        crystal_file << buffer;
    }

    crystal_file.close();

    //Initialize the RNG
    rng.seed(settings.SEED);

    clock_t start = clock();
    int n = 0;

    
    out_file << "X0\tY0\tZm\tE\tTHETA\tPHI\tt\tsteps\tmax_n\tr_min\n";

    if(settings.NWRITX == 666)
    {
        if(settings.NUMCHA==1 || settings.NWRITY == 777)
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
        else if(settings.NWRITY == 666)
        {
            debug_file << "Running Montecarlo " << std::endl;
            std::cout << "Running Montecarlo " << std::endl;
            montecarloscat(lattice, &n);
        }
        else
        {
            chainscat(lattice, &n);
        }
    }
    else
    {
        n = 1;
        //Testing stuff
        Vec3d loc;
        int r = 2;
        int var = pow(2*r + 1, 3);
        for(int n = 0; n < var; n++)
        {
            index_to_loc(n, loc);
            // debug_file << loc[0] << " " << loc[1] << " " << loc[2] << std::endl;
        }
    }
    

    double dt = ( (double)clock() - start ) / CLOCKS_PER_SEC;
    dt /= n;
    //Convert to ms;
    dt *= 1000;

    printf ( "%f\n", dt );

    debug_file << "Time per run: "<< dt <<"ms"<<std::endl;
    out_file.close();
    debug_file.close();
    traj_file.close();
    return 0;
}
