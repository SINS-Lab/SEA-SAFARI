#include "ion.h"

#include "potentials.h" // needed for image effect on initial KE
#include "safio.h"      // settings
#include "temps.h"      // thermalize_ion, thermalize

#include <algorithm>  // std::transform
#include <cmath>      // sqrt, cos, sin, tan, etc
#include <functional> // std::minus

double zeros[3] = {0, 0, 0};

void Ion::set_KE(double E0, double theta0, double phi0, double x, double y)
{
    //Start by resetting the ion to its initial state.
    reset();
    //Set the energy of the ion
    this->E0 = E0;
    //Attempt to thermalize the ion
    if (settings.ESIZE > 0)
    {
        thermaize_ion(*this);
    }
    this->Theta0 = theta0;
    this->Phi0 = phi0;

    atom = &settings.ion;
    //TODO lookup table for atomic symbols...
    //TODO charge config somewhere.
    q = 1;

    //Convert incoming angle to radians.
    theta0 = theta0 * M_PI / 180;
    phi0 = phi0 * M_PI / 180;

    double p0 = sqrt(atom->two_mass * E0);
    double p_trans = p0 * sin(theta0);

    //This is the initial momentum, before surface effects.
    double p_x0 = p_trans * cos(phi0);
    double p_y0 = p_trans * sin(phi0);
    double p_z0 = -p0 * cos(theta0);

    //Impact parameters offsets, this aims at the impact point.
    r[0] = r_t[0] = -settings.Z1 * tan(theta0) * cos(phi0) + x;
    r[1] = r_t[1] = -settings.Z1 * tan(theta0) * sin(phi0) + y;
    r[2] = r_t[2] = settings.Z1;

    //Reset the forces, this cleans up some of the debug output.
    std::copy(zeros, zeros + 3, dp_dt);
    std::copy(zeros, zeros + 3, dp_dt_t);

    //Set the "initial" location to the targetted impact point.
    r_0[0] = x;
    r_0[1] = y;
    r_0[2] = settings.Z1;

    //If we have image effect, account for that here.
    if (settings.use_image)
    {
        //The image part is subtracted here, as we have already assigned
        //the sign outside the sqrt function, the output of Vi_z is negative.
        p_z0 = -sqrt((p_z0 * p_z0) - (atom->two_mass * Vi_z(settings.Z1, q)));
    }

    //Set the initial momentum of the ion
    p[0] = p_x0;
    p[1] = p_y0;
    p[2] = p_z0;
}

void Ion::reset()
{
    E0 = 0;
    Theta0 = 0;
    Phi0 = 0;
    steps = 0;
    time = 0;
    max_n = 0;
    near = 0;
    r_min = 1e3;
    rr_min_find = 1e3;
    last_index = -1;
    total_near = 0;
    V = 0;
}

void Site::reset()
{
    //Reset positions and momenta
    std::copy(r_0, r_0 + 3, r);
    std::copy(r_0, r_0 + 3, r_t);
    std::copy(r_0, r_0 + 3, r_u);
    std::copy(p_0, p_0 + 3, p);
    last_step = -1;
    left_origin = false;
    unbound = false;
    sputter_tick = -1;

    reset_forces();

    hameq_tick = -1;

    //Reset our tracked sites last steps too.
    for (int i = 0; i < near; i++)
    {
        if (near_sites[i]->last_ion != last_ion)
            near_sites[i]->last_step = -1;
    }
    log_z = r_0[2];

    validate();

    //Thermalize the site
    thermaize(this);
}

void Site::write_info()
{
    debug_file << "Atom: " << atom->symbol << std::endl;
    debug_file << "r_0: " << r_0[0] << " " << r_0[1] << " " << r_0[2] << std::endl;
    debug_file << "r  : " << r[0] << " " << r[1] << " " << r[2] << std::endl;
    debug_file << "p  : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    debug_file << "r_t: " << r_t[0] << " " << r_t[1] << " " << r_t[2] << std::endl;
    debug_file << "F  : " << dp_dt[0] << " " << dp_dt[1] << " " << dp_dt[2] << std::endl;
    debug_file << "F_t: " << dp_dt_t[0] << " " << dp_dt_t[1] << " " << dp_dt_t[2] << std::endl;
    debug_file << "pnt: " << r << " " << p << " " << r_0 << " " << last_ion << " " << last_step << " " << last_update << std::endl;
}