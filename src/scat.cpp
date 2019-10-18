#include "scat.h"
#include "safio.h"
#include "ion.h"

double frand()
{
    double max = rng.max();
    return rng()/max;
}

void montecarloscat(Lattice &lattice, int *num)
{
    double x_size = (settings.XSTOP - settings.XSTART);
    double y_size = (settings.YSTOP - settings.YSTART);

    for(int i = 0; i < settings.NUMCHA; i++)
    {
        double rx = frand();
        double ry = frand();
        double x = settings.XSTART + x_size * rx;
        double y = settings.YSTART + y_size * ry;
        Ion ion;
        ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x, y);
        ion.index = i;
        traj(ion, lattice, false, false);
    }
    *num = settings.NUMCHA;
}

void gridscat(Lattice &lattice, int *num)
{
    bool log = settings.NUMCHA == 1;
    int n = 0;
    for(double x = settings.XSTART; x <= settings.XSTOP; x+=settings.XSTEP)
        for(double y = settings.YSTART; y <= settings.YSTOP; y+=settings.YSTEP)
        {
            Ion ion;
            ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x, y);
            ion.index = n;
            n++;
            traj(ion, lattice, log, log);
            if(log) goto out;
        }
out:
    *num = n;
    return;
}

void chainscat(Lattice &lattice, int *num)
{

}
