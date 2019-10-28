#include "temps.h"
#include <cmath>
#include "safio.h"

#define boltz 8.617e-5 /* Units of eV/K */

std::default_random_engine temperature_rng;

double frand(std::default_random_engine &rng)
{
    double max = rng.max();
    return rng()/max;
}

void init_temps()
{
    init_temp_seed();
    //Convert from temperature to energy.
    double energy = boltz * settings.TEMP;
    for(Atom &a : settings.ATOMS)
    {
        if(settings.TEMP > 0)
        {
            for(int i = 0; i<3; i++)
            {
                //dx = sqrt(2E/k)
                a.dev_r[i] = sqrt(2.0*energy/a.spring[i]);
                //dp = sqrt(2mE)
                a.dev_p[i] = sqrt(2.0*a.mass*energy);
            }
        }
        else
        {
            for(int i = 0; i<3; i++)
            {
                //Both 0 at 0K
                a.dev_p[i] = 0;
                a.dev_r[i] = 0;
            }
        }
    }
}

void init_temp_seed()
{
    temperature_rng.seed(settings.SEED * UINT_FAST32_MAX);
}

void thermaize(Site &site)
{
    if(settings.TEMP > 0)
    {
        //2 random numbers
        double r1 = 0, r2 = 0;
        for(int i = 0; i<3; i++)
        {
            r1 = frand(temperature_rng);
            //Ensure this isn't 0
            while(r1==0) r1 = frand(temperature_rng);
            //Needed to not be 0 so we could do this.
            r1 = sqrt(-log(r1));
            //Other RNG
            r2 = frand(temperature_rng);
            //Adjusted position
            site.r[i] = site.r_0[i] + site.atom->dev_r[i]*r1*cos(2*M_PI*r2);
            site.p[i] = site.p_0[i] + site.atom->dev_p[i]*r1*sin(2*M_PI*r2);
        }
    }
}

void thermaize_ion(Ion &ion)
{

}