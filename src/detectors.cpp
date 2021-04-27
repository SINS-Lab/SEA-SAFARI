
#include "lattice.h"
#include "potentials.h" // the potential used
#include "detector.h"  // Our base header
#include "safio.h"     // For settings


void Detector::log(std::ofstream &out_file, Site &ion, Lattice *lattice,
                   bool stuck, bool buried, bool froze, bool off_edge, bool discont,
                   bool ignore_bounds)
{
    double psq = sqr(ion.p);
    double mx2 = ion.atom->two_mass;
    double E = psq / mx2;
    // z-momentum squared, exit theta, exit phi
    double pzz = 0, theta, phi;
    /**
     * Checks failure flags, and sets energy accordingly
     * 
     * In all failure cases, theta = 0, phi = 90
     * 
     * Energy failure flags as follows:
     * 
     * -5   : Ion exited at phi out of absolute max detector bounds
     * -10  : ion got stuck due to image charge.
     * -100 : got stuck, ie E too low
     * -200 : got buried, ie z below -BDIST
     * -300 : froze, ie took too many steps
     * -400 : off edge, ie left crystal via x or y
     * -500 : discont, had a discontinuity in E, so was dropped.
     * 
     * Note that -5 condition is only applied if settings.detector_type is greater than 0
     */
    if (stuck)
    {
        theta = 0;
        phi = 90;

        if (ion.r[2] > 0)
        {
            E = -10;
            lattice->trapped_num++;
        }
        else
        {
            E = -100;
            lattice->stuck_num++;
        }
    }
    else if (buried)
    {
        theta = 0;
        phi = 90;
        E = -200;
        lattice->buried_num++;
    }
    else if (froze)
    {
        theta = 0;
        phi = 90;
        E = -300;
        lattice->froze_num++;
    }
    else if (off_edge)
    {
        theta = 0;
        phi = 90;
        E = -400;
        lattice->left_num++;
    }
    else if (discont)
    {
        theta = 0;
        phi = 90;
        E = -500;
        lattice->err_num++;
    }
    else if (ion.site_site_intersects)
    {
        theta = 0;
        phi = 90;
        E = -900;
        lattice->intersections++;
    }
    else
    {
        double px = ion.p[0];
        double py = ion.p[1];
        double pz = ion.p[2];
        // Find the momentum at infinity
        if (settings.use_image and ion.q != 0)
        {
            // Image charge would pull it towards surface, this accounts
            // for that effect.
            pzz = (pz * pz) + (mx2 * Vi_z(settings.Z1, ion.q));
            pzz = pzz < 0 ? -sqrt(-pzz) : sqrt(pzz);
            // Recalulate this, as pz has changed
            // We are fine with pzz being -ve, as that case
            // will be dropped due to ion not escaping.
            psq = pzz * pzz + px * px + py * py;
        }
        else
        {
            // No image, so this is the same as it was.
            pzz = pz;
        }
        double p = sqrt(psq);
        // Ion is not escaping.
        if (pzz <= 0 || p <= 0)
        {
            theta = 0;
            phi = 90;
            E = -10;
            lattice->trapped_num++;
        }
        else
        {
            // Recalculate E, incase image affected it
            E = psq / mx2;

            // calculate theta, depends on pz
            theta = acos(pzz / p) * RAD2DEG;

            if (px == 0 && py == 0)
            {
                // phi isn't well defined here,
                // so just set it as same as incoming
                phi = settings.PHI0;
            }
            else
            {
                // Calculate phi, depends on px, py
                phi = atan2(py, px) * RAD2DEG;
            }
        }
    }

    bool did_hit = E > -10;

    if (did_hit && !hit(E, theta, phi) && !ignore_bounds)
    {
        theta = 0;
        phi = 90;
        E = -5;
        lattice->undetectable_num++;
        did_hit = false;
    }
    if (did_hit || settings.save_errored)
    {
        /**
             * This uses the default saving behaviour
             */
        char buffer[200];
        // first stuff it in the buffer
        sprintf(buffer, "%f\t%f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\t%.3f\n",
                ion.r_0[0], ion.r_0[1], ion.r_0[2],
                E, theta, phi,
                ion.thermal_seed, ion.weight,
                ion.max_n, ion.r_min, ion.steps,
                ion.Eerr_max, ion.time);
        // Then save it
        mutx.lock();
        out_file << buffer << std::flush;
        mutx.unlock();
    }
    // We use this as a after saving.
    ion.index = did_hit;
}

void SpectrumDetector::log(std::ofstream &out_file, Site &ion, Lattice *lattice,
                bool stuck, bool buried, bool froze, bool off_edge, bool discont,
                bool ignore_bounds)
{
    double psq = sqr(ion.p);
    double mx2 = ion.atom->two_mass;
    double E = psq / mx2;
    // z-momentum squared, exit theta, exit phi
    double pzz = 0, theta, phi;
    /**
    * Checks failure flags, and sets energy accordingly
    * 
    * In all failure cases, theta = 0, phi = 90
    * 
    * Energy failure flags as follows:
    * 
    * -5   : Ion exited at phi out of absolute max detector bounds
    * -10  : ion got stuck due to image charge.
    * -100 : got stuck, ie E too low
    * -200 : got buried, ie z below -BDIST
    * -300 : froze, ie took too many steps
    * -400 : off edge, ie left crystal via x or y
    * -500 : discont, had a discontinuity in E, so was dropped.
    * 
    * Note that -5 condition is only applied if settings.detector_type is greater than 0
    */
    if (stuck)
    {
        theta = 0;
        phi = 90;

        if (ion.r[2] > 0)
        {
            E = -10;
            lattice->trapped_num++;
        }
        else
        {
            E = -100;
            lattice->stuck_num++;
        }
    }
    else if (buried)
    {
        theta = 0;
        phi = 90;
        E = -200;
        lattice->buried_num++;
    }
    else if (froze)
    {
        theta = 0;
        phi = 90;
        E = -300;
        lattice->froze_num++;
    }
    else if (off_edge)
    {
        theta = 0;
        phi = 90;
        E = -400;
        lattice->left_num++;
    }
    else if (discont)
    {
        theta = 0;
        phi = 90;
        E = -500;
        lattice->err_num++;
    }
    else if (ion.site_site_intersects)
    {
        theta = 0;
        phi = 90;
        E = -900;
        lattice->intersections++;
    }
    else
    {
        double px = ion.p[0];
        double py = ion.p[1];
        double pz = ion.p[2];
        // Find the momentum at infinity
        if (settings.use_image and ion.q != 0)
        {
            // Image charge would pull it towards surface, this accounts
            // for that effect.
            pzz = (pz * pz) + (mx2 * Vi_z(settings.Z1, ion.q));
            pzz = pzz < 0 ? -sqrt(-pzz) : sqrt(pzz);
            // Recalulate this, as pz has changed
            // We are fine with pzz being -ve, as that case
            // will be dropped due to ion not escaping.
            psq = pzz * pzz + px * px + py * py;
        }
        else
        {
            // No image, so this is the same as it was.
            pzz = pz;
        }
        double p = sqrt(psq);
        // Ion is not escaping.
        if (pzz <= 0 || p <= 0)
        {
            theta = 0;
            phi = 90;
            E = -10;
            lattice->trapped_num++;
        }
        else
        {
            // Recalculate E, incase image affected it
            E = psq / mx2;

            // calculate theta, depends on pz
            theta = acos(pzz / p) * RAD2DEG;

            if (px == 0 && py == 0)
            {
                // phi isn't well defined here,
                // so just set it as same as incoming
                phi = settings.PHI0;
            }
            else
            {
                // Calculate phi, depends on px, py
                phi = atan2(py, px) * RAD2DEG;
            }
        }
    }
    bool did_hit = E > 0;
    if (did_hit)
    {
        float dE = settings.E0 / ERES;
        int E_bin = E / dE;

        if(theta - dtheta < 0)
        {
            modulo = 180;
        }
        else
        {
            modulo = 360;
        }
        double phi_ = phi;
        while(phi_ < 0)
            phi_ += modulo;

        int theta_bin = theta * (TRES / 180.0);
        int phi_bin = phi_ * (PRES / 360.0);

        if(E_bin < ERES && E_bin >= 0 && 
            theta_bin >= 0 && theta_bin < TRES && 
            phi_bin >= 0 && phi_bin < PRES)
        {
            logNum++;
            int num = counts[E_bin][theta_bin][phi_bin] + 1;
            counts[E_bin][theta_bin][phi_bin] = num;
            //if(num > 3) std::cout << "Another one! "<<num<<"\n" << std::flush;
        }
        else
        {
            std::cout << "Error finding bin? "<< E_bin<< " "<< theta_bin << " " << phi_bin << " " << phi_ << '\n' << std::flush;
        }
        if(logNum > saveNum)
        {
            save(out_file);
            logNum = 0;
        }
    }
    // We use this as a after saving.
    ion.index = did_hit;
}


void Detector::finish(std::ofstream &out_file)
{
}

void SpectrumDetector::save(std::ofstream &out_file)
{
    // Then save it
    mutx.lock();

    out_file.close();
    std::string filename = settings.output_name + ".spec";
    out_file.open(filename, std::ofstream::trunc);

    float dE = settings.E0 / ERES;
    float dT = (180 / TRES);
    // We will do this the long way for now, can speed it up later I guess?
    for(int i = 0; i < ERES; i++)
    {
        out_file << dE * i;
        for(int j = 0; j < TRES;j++)
        {
            out_file << '\n' << '\t' << dT * j;
            for(int k = 0; k < PRES; k++)
            {
                out_file << '\t';
                out_file << counts[i][j][k];
            }
        }
        out_file << '\n';
    }
    out_file << std::flush;
    mutx.unlock();
}
