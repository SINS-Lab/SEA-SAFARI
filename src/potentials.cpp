#include "potentials.h"
#include "safio.h"
#include <math.h>
#include <iostream>

double ** Vr_r_cache;
double ** dVr_dr_cache;
double r_max;
double dr_min;
double z_max;
double z_min;
int n_rmax;

double dVr_dr_init(double r, int n)
{
    if(settings.IPOT==1)
    {
        double a,b,c,d;
        //The potpars are in groups of 4,
        //in order of the listed basis atoms
        //The first atom is listed as 1.
        int index = (n-1)*4;
        a = settings.POTPAR[index];
        b = settings.POTPAR[index + 1];
        c = settings.POTPAR[index + 2];
        d = settings.POTPAR[index + 3];
        return -b*a*exp(-b*r) - d*c*exp(-d*r);
    }
    else
    {
        debug_file << "ERROR WITH dVr_dr" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double Vr_r_init(double r, int n)
{
    if(settings.IPOT==1)
    {
        double a,b,c,d;
        //The potpars are in groups of 4,
        //in order of the listed basis atoms
        //The first atom is listed as 1.
        int index = (n-1)*4;
        a = settings.POTPAR[index];
        b = settings.POTPAR[index + 1];
        c = settings.POTPAR[index + 2];
        d = settings.POTPAR[index + 3];
        return a*exp(-b*r) + c*exp(-d*r);
    }
    else
    {
        debug_file << "ERROR WITH Vr_r" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void init_potentials()
{
    double r;
    dr_min = settings.DR_MIN_TAB;
    r_max = settings.R_MAX;
    if(dr_min == 0) return;
    n_rmax = r_max / dr_min;

    Vr_r_cache = new double*[settings.NTYPES];
    dVr_dr_cache = new double*[settings.NTYPES];
    
    for(int n = 0; n<settings.NTYPES; n++)
    {
        Vr_r_cache[n] = new double[n_rmax];
        Vr_r_cache[n][0] = 0;
        dVr_dr_cache[n] = new double[n_rmax];
        dVr_dr_cache[n][0] = 0;
    }

    //Start at 1, as these guys are not defined for r=0 anyway.
    for(int i = 1; i<n_rmax; i++)
    {
        r = dr_min * i;
        for(int n = 0; n<settings.NTYPES; n++)
        {
            dVr_dr_cache[n][i] = dVr_dr_init(r, n+1);
            Vr_r_cache[n][i] = Vr_r_init(r, n+1);
        }
    }
}

double interp_r(double r, double *arr)
{
    //Index for before r
    int i_bef = (int)(r/dr_min);
    //Index for after r
    int i_aft = i_bef+1;
    //Anything out of bounds of table is 0.
    if(i_aft >= n_rmax) return 0;
    //Fraction of way between before and after
    double frac = (r - dr_min*i_bef)/dr_min;
    //Value before
    double bef = arr[i_bef];
    //Value after
    double aft = arr[i_aft];
    //Linearly interpolate between the two.
    return (1.0 - frac)*bef + frac * aft;
}

double Vr_r(double r, int n)
{
    if(dr_min == 0) return Vr_r_init(r, n);
    return interp_r(r, Vr_r_cache[n-1]);
}

double dVr_dr(double r, int n)
{
    if(dr_min == 0) return dVr_dr_init(r, n);
    return interp_r(r, dVr_dr_cache[n-1]);
}

double Vi_z(double z, int q)
{
    double z_min = settings.PIMPAR[0];
    double v_min = settings.PIMPAR[1];
    if(settings.IIMPOT == 1)
    {
        if(z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25*eqsr*q;
            double eq_v = eq/v_min;
            return -eq/sqrt(dz*dz + eq_v*eq_v);
        }
        else 
        {
            return -v_min;
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
    double z_min = settings.PIMPAR[0];
    double v_min = settings.PIMPAR[1];
    if(settings.IIMPOT == 1)
    {
        if(z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25*eqsr*q;
            double eq_v = eq/v_min;
            return eq*dz/pow(dz*dz + eq_v*eq_v, 1.5);
        }
        return 0;
    }
    else
    {
        debug_file << "ERROR WITH dVi_dz" << std::endl;
        exit(EXIT_FAILURE);
    }
    return 0;
}

void apply_friction(Ion &ion, double* F)
{
    double z = ion.r[2];
    //TODO some better "above surfaceness"
    if(z > 0) return;

    double vx = ion.p[0]/ion.atom->mass;
    double vy = ion.p[1]/ion.atom->mass;
    double vz = ion.p[2]/ion.atom->mass;

    double v_sq = vx*vx + vy*vy + vz*vz;
    double v = sqrt(v_sq);
    //No velocity, no friction.
    if(v == 0) return;

    //Components of the velocity direction.
    vx/=v;
    vy/=v;
    vz/=v;

    //assume magnitude of friction is:
    //Av + Bv^2, where v is magnitude of velocity.
    double fric = (settings.F_a*v + settings.F_b*v_sq);
    //Direction of friction is opposite of velocity.
    F[0] -= vx*fric;
    F[1] -= vy*fric;
    F[2] -= vz*fric;
}
