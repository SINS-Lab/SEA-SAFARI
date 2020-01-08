#include "potentials.h"

#include "safio.h"   //settings
#include "safari.h"  //exit_fail

#include <math.h>    //exp, sqrt, etc

double ** Vr_r_cache;
double ** dVr_dr_cache;
double r_max;
double dr_min;
double z_max;
double z_min;
int n_rmax;

double dVr_dr_init(double r, int n)
{
    if(settings.binary_potential_type==1)
    {
        double a,b,c,d;
        //The potpars are in groups of 4,
        //in order of the listed basis atoms
        //The first atom is listed as 1.
        int index = (n-1)*4;
        a = settings.binary_potential_parameters[index];
        b = settings.binary_potential_parameters[index + 1];
        c = settings.binary_potential_parameters[index + 2];
        d = settings.binary_potential_parameters[index + 3];
        return -b*a*exp(-b*r) - d*c*exp(-d*r);
    }
    else
    {
        exit_fail("ERROR WITH dVr_dr");
    }
    return 0;
}

double Vr_r_init(double r, int n)
{
    if(settings.binary_potential_type==1)
    {
        double a,b,c,d;
        //The potpars are in groups of 4,
        //in order of the listed basis atoms
        //The first atom is listed as 1.
        int index = (n-1)*4;
        a = settings.binary_potential_parameters[index];
        b = settings.binary_potential_parameters[index + 1];
        c = settings.binary_potential_parameters[index + 2];
        d = settings.binary_potential_parameters[index + 3];
        return a*exp(-b*r) + c*exp(-d*r);
    }
    else
    {
        exit_fail("ERROR WITH Vr_r");
    }
    return 0;
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
    if(settings.image_potential_type == 0)
    {
        //This shouldn't be used for below the surface.
        if(z < settings.Z1) return 0;
        
        //Type 0: Only apply to the entry and exit trajectories.
        //Does not actually do anything during the numerical integration
        double z_min = settings.image_parameters[0];
        double v_min = settings.image_parameters[1];
        if(z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25*eqsr*q;
            double eq_v = eq/v_min;
            return -eq/sqrt(dz*dz + eq_v*eq_v);
        }
        else 
        {
            return -q * v_min;
        }
    }
    else if(settings.image_potential_type == 1)
    {
        //Type 1: Saturated image potential
        //Saturates based on the given z_min and v_min
        double z_min = settings.image_parameters[0];
        double v_min = settings.image_parameters[1];
        if(z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25*eqsr*q;
            double eq_v = eq/v_min;
            return -eq/sqrt(dz*dz + eq_v*eq_v);
        }
        else 
        {
            return -q * v_min;
        }
    }
    else
    {
        exit_fail("ERROR WITH Vi_z");
    }
    return 0;
}

double dVi_dz(double z, int q)
{
    if(settings.image_potential_type == 0)
    {
        //Type 0: Only apply to the entry and exit trajectories.
        //Does not actually do anything during the numerical integration
        return 0;
    }
    if(settings.image_potential_type == 1)
    {
        double z_min = settings.image_parameters[0];
        double v_min = settings.image_parameters[1];
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
        exit_fail("ERROR WITH dVi_dz");
    }
    return 0;
}

double electron_density(Lattice &lattice, Ion &ion)
{
    //TODO see https://github.com/SINS-Lab/SAFARI/blob/10305e6f9ee597e89a6df7acfede3554371f42bc/src/hameqinel.f
    // it has calculations for electron density.
    return 0;
}

double apply_friction(Lattice &lattice, Ion &ion, double* F, double dt)
{
    double vx = ion.p[0]/ion.atom->mass;
    double vy = ion.p[1]/ion.atom->mass;
    double vz = ion.p[2]/ion.atom->mass;

    double v_sq = vx*vx + vy*vy + vz*vz;
    double v = sqrt(v_sq);
    //No velocity, no friction.
    if(v == 0) return 0;

    //Components of the velocity direction.
    vx/=v;
    vy/=v;
    vz/=v;

    double density = electron_density(lattice, ion);
    //No electron density, no friction
    if(density == 0) return 0;

    //assume magnitude of friction is:
    //Av + Bv^2, where v is magnitude of velocity.
    //This is then scaled by electron density for the location
    double fric = (settings.F_a*v + settings.F_b*v_sq)
                      * density;
    //Direction of friction is opposite of velocity.
    F[0] -= vx*fric;
    F[1] -= vy*fric;
    F[2] -= vz*fric;

    //TODO return effective "potential" from this interaction.
    return 0;
}
