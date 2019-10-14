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

safio settings;
std::ofstream out_file;
std::ofstream debug_file;

int main()
{
    //Load the input file
    settings.load();

    char buffer[100];

    lattice l;
    l.build_lattice();
    std::ofstream myfile;
    myfile.open ("crystal.input");

    int num = l.sites.size();

    for(int i = 0; i<num; i++)
    {
        site s = l.sites[i];
        atom a = l.atoms[i];
        sprintf(buffer, "%f\t%f\t%f\t%f\t%f\n",s[0],s[1],s[2],a.charge,a.mass);
        myfile << buffer;
    }
    myfile.close();


    clock_t start = clock();

    int n = 0;
    srand(1);
    bool log = settings.NUMCHA == 1;
    for(double x = settings.XSTART; x <= settings.XSTOP; x+=settings.XSTEP)
        for(double y = settings.YSTART; y <= settings.YSTOP; y+=settings.YSTEP)
        {
            double x_0 = rand()%20;
            double y_0 = rand()%20;

            x_0 = 0;
            y_0 = 16;

            ion ion;
            debug_file << "ion: " << n << " x: " << x_0 << " y: "<< y_0 << std::endl;
            ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x_0, y_0);
            ion.index = n;
            debug_file << "fire: "  << ion[0] << " " << ion[1] << " " << ion[2] << std::endl;
            traj(ion, l, log);
            n++;
            if(log || n>1000)
                goto out;
        }

out:
    double dt = ( (double)clock() - start ) / CLOCKS_PER_SEC;
    dt /= n;
    printf ( "%f\n", dt );

    debug_file << "Time per run: "<<dt<<"ms"<<std::endl;
    out_file.close();
    debug_file.close();
    return 0;
}
