#include "scat.h"
#include "safio.h"
#include "ion.h"
#include <stdlib.h>     /* srand, rand */

void montecarloscat(lattice &lattice, int *num)
{
    double x_size = (settings.XSTART - settings.XSTOP)/1e4;
    double y_size = (settings.YSTART - settings.YSTOP)/1e4;

    for(int i = 0; i < settings.NUMCHA; i++)
    {
        int rx = rand()%10000;
        int ry = rand()%10000;
        double x = settings.XSTART + x_size * rx;
        double y = settings.YSTART + y_size * ry;
        ion ion;
        ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, y, x);
        *num = i;
        ion.index = *num;
        traj(ion, lattice, false);
    }
}

void gridscat(lattice &lattice, int *num)
{
    bool log = settings.NUMCHA == 1;
    for(double x = settings.XSTART; x <= settings.XSTOP; x+=settings.XSTEP)
        for(double y = settings.YSTART; y <= settings.YSTOP; y+=settings.YSTEP)
        {
            ion ion;
            ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, y, x);
            ion.index = *num;
            *num = *num + 1;
            traj(ion, lattice, log);
            if(log) goto out;
        }
out:
    return;
}

void chainscat(lattice &lattice, int *num)
{

}
