#include "potentials.h"
#include "safio.h"
#include <math.h>
#include <iostream>


double Vr_r(double r[], int n[], int num)
{
    double out = 0;
    if(settings.IPOT==1)
    {
        double a,b,c,d,x;
        for(int i = 0; i < num; i++)
        {
            //The potpars are in groups of 4,
            //in order of the listed basis atoms
            //The first atom is listed as 1.
            int index = (n[i]-1)*4;
            a = settings.POTPAR[index];
            b = -settings.POTPAR[index + 1];
            c = settings.POTPAR[index + 2];
            d = -settings.POTPAR[index + 3];
            x = r[i];
            out += a*exp(b*x) + c*exp(d*x);
        }
    }
    else
    {
        debug_file << "ERROR WITH Vr_r" << std::endl;
        exit(EXIT_FAILURE);
    }
    return out;
}

double * dVr_dr(double r[], int n[], int num)
{
    double arr[num];

    if(settings.IPOT==1)
    {
        double a,b,c,d,x;
        for(int i = 0; i < num; i++)
        {
            //The potpars are in groups of 4,
            //in order of the listed basis atoms
            //The first atom is listed as 1.
            int index = (n[i]-1)*4;
            a = settings.POTPAR[index];
            b = -settings.POTPAR[index + 1];
            c = settings.POTPAR[index + 2];
            d = -settings.POTPAR[index + 3];
            x = r[i];
            arr[i] = b*a*exp(b*x) + d*c*exp(d*x);
        }
    }
    else
    {
        debug_file << "ERROR WITH dVr_dr" << std::endl;
        exit(EXIT_FAILURE);
    }
    double *out = arr;
    return out;
}

double Vi_z(double z, int q)
{
    double z_min = settings.PIMPAR[1];
    double v_min = settings.PIMPAR[2];
    if(settings.IIMPOT == 1)
    {
        if(z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25*eqsr*q;
            double eq_v = eq/v_min;
            return eq/sqrt(dz*dz + eq_v*eq_v);
        }
    }
    else
    {
        debug_file << "ERROR WITH Vi_z" << std::endl;
        exit(EXIT_FAILURE);
    }
    return -v_min;
}

double dVi_dz(double z, int q)
{
    double z_min = settings.PIMPAR[1];
    double v_min = settings.PIMPAR[2];
    if(settings.IIMPOT == 1)
    {
        if(z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25*eqsr*q;
            double eq_v = eq/v_min;
            return eq*dz/pow(dz*dz + eq_v*eq_v, 1.5);
        }
    }
    else
    {
        debug_file << "ERROR WITH dVi_dz" << std::endl;
        exit(EXIT_FAILURE);
    }
    return 0;
}
