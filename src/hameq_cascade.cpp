#include "hameq.h"

#include "safio.h"      // for settings
#include "traj.h"       // for nearest check method
#include "potentials.h" // for forces/potentials
#include "vec_math.h"   // General maths help
#include "safari.h"     // for exit_fail
#include <cmath>        // sqrt

void apply_hameq(std::vector<Ion*> &ions, Lattice *lattice, double dt)
{
    int num = ions.size();
    for(int i = 0; i<num; i++)
    {
        Ion &ion = *ions[i];
        if (!settings.useEinsteinSprings)
        {
            update_site(ion, dt);
            // In this case, we update all of the sites
            if (settings.dynamicNeighbours)
            {
                int num = lattice->active_sites.size();
                for (int i = 0; i < num; i++)
                {
                    Site *s = lattice->active_sites[i];
                    update_site(*s, dt);
                    check_sputter(ion, s);
                    s->last_update = ion.steps;
                }
            }
            else
            {
                // Otherwise, just consider ones tracked by the ion
                for (int i = 0; i < ion.near; i++)
                {
                    Site *s = ion.near_sites[i];
                    update_site(*s, dt);
                    s->last_update = ion.steps;
                    if (!settings.rigidBounds)
                        for (int j = 0; j < s->near; j++)
                        {
                            Site *s2 = s->near_sites[j];
                            if (s2->last_update != ion.steps)
                            {
                                update_site(*s2, dt);
                                check_sputter(ion, s);
                            }
                            s2->last_update = ion.steps;
                        }
                }
            }
        }
        else
        {
            update_site(ion, dt);
            for (int i = 0; i < ion.near; i++)
            {
                Site *s = ion.near_sites[i];
                update_site(*s, dt);
                check_sputter(ion, s);
            }
        }
    }
}

void run_hameq(std::vector<Ion *> &ions, Lattice *lattice, double dt, bool predicted, double *dr_max)
{

}
