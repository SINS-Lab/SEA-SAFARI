#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <cstdio>
#include "safio.h"
#include "vec_math.h"
#include "ion.h"
#include <time.h>
#include "scat.h"

safio settings;
std::ofstream out_file;
std::ofstream debug_file;

int main()
{
    //Load the input file
    settings.load();

    char buffer[100];

    lattice lattice;
    lattice.build_lattice();
    std::ofstream myfile;
    myfile.open ("crystal.input");

    int num = lattice.sites.size();

    for(int i = 0; i<num; i++)
    {
        site s = lattice.sites[i];
        atom a = lattice.atoms[i];
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",s[0],s[1],s[2],a.charge,a.mass);
        myfile << buffer;
    }
    myfile.close();


    clock_t start = clock();
    int n = 0;
    srand(settings.SEED);

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

    double dt = ( (double)clock() - start ) / CLOCKS_PER_SEC;
    dt /= n;
    printf ( "%f\n", dt );

    debug_file << "Time per run: "<<dt<<"ms"<<std::endl;
    out_file.close();
    debug_file.close();
    return 0;
}
