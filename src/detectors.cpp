
#include "lattice.h"
#include "potentials.h"// the potential used
#include "detector.h"  // Our base header
#include "safio.h"     // For settings

void LogCheck::checkConditionsForLog(Site &ion, Lattice *lattice,
                   bool stuck, bool buried, bool froze, bool off_edge, bool discont,
                   bool ignore_bounds)
{
    double psq = sqr(ion.p);
    double mx2 = ion.atom->two_mass;
    E = psq / mx2;
    // z-momentum squared
    double pzz = 0;
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
     */
    if (stuck)
    {
        theta = 0;
        phi = 90;

        if (ion.r[2] > 0)
        {
            E = -10;
            if(lattice != NULL) lattice->trapped_num++;
        }
        else
        {
            E = -100;
            if(lattice != NULL) lattice->stuck_num++;
        }
    }
    else if (buried)
    {
        theta = 0;
        phi = 90;
        E = -200;
        if(lattice != NULL) lattice->buried_num++;
    }
    else if (froze)
    {
        theta = 0;
        phi = 90;
        E = -300;
        if(lattice != NULL) lattice->froze_num++;
    }
    else if (off_edge)
    {
        theta = 0;
        phi = 90;
        E = -400;
        if(lattice != NULL) lattice->left_num++;
    }
    else if (discont)
    {
        theta = 0;
        phi = 90;
        E = -500;
        if(lattice != NULL) lattice->err_num++;
    }
    else if (ion.site_site_intersects)
    {
        theta = 0;
        phi = 90;
        E = -900;
        if(lattice != NULL) lattice->intersections++;
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
            if(lattice != NULL) lattice->trapped_num++;
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
}

void Detector::log(std::ofstream &out_file, Site &ion, Lattice *lattice,
                   bool stuck, bool buried, bool froze, bool off_edge, bool discont,
                   bool ignore_bounds)
{
    LogCheck check;
    check.checkConditionsForLog(ion, lattice, stuck, buried, froze, off_edge, discont, ignore_bounds);
    double E = check.E;
    double theta = check.theta;
    double phi = check.phi;

    bool did_hit = E > -10;

    if (did_hit && !hit(E, theta, phi) && !ignore_bounds)
    {
        theta = 0;
        phi = 90;
        E = -5;
        if(lattice != NULL) lattice->undetectable_num++;
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

int SpectrumDetector::getThetaBin(double _theta, double _phi)
{
    double theta_min = theta - dtheta;
    double theta_max = theta + dtheta;

    if(_theta < theta_min || _theta > theta_max) return -1;
    double dT = (2 * dtheta / TRES);
    return (_theta - theta_min) / dT;
}
int SpectrumDetector::getPhiBin(double _theta, double _phi)
{
    double phi_min = phi - dphi;
    double phi_max = phi + dphi;

    if(_phi < phi_min || _phi > phi_max) return -1;
    double dP = (2 * dphi / PRES);
    return (_phi - phi_min) / dP;
}
int SpectrumDetector::getEnergyBin(double E)
{
    float dE = settings.E0 / ERES;
    return (E - e_min) / dE;
}

void SpectrumDetector::log(std::ofstream &out_file, Site &ion, Lattice *lattice,
                bool stuck, bool buried, bool froze, bool off_edge, bool discont,
                bool ignore_bounds)
{
    LogCheck check;
    check.checkConditionsForLog(ion, lattice, stuck, buried, froze, off_edge, discont, ignore_bounds);
    double E = check.E;
    double theta = check.theta;
    double phi = check.phi;

    bool did_hit = E > 0;
    if (did_hit)
    {
        did_hit = false;
        int E_bin = getEnergyBin(E);
        int theta_bin = getThetaBin(theta, phi);
        int phi_bin = getPhiBin(theta, phi);

        if(E_bin < ERES && E_bin >= 0 && 
            theta_bin >= 0 && theta_bin < TRES && 
            phi_bin >= 0 && phi_bin < PRES)
        {
            did_hit = true;
            logNum++;
            total_counts++;
            int num = counts[E_bin][theta_bin][phi_bin] + 1;
            if(num > big_bin[0])
            {
                 big_bin[0] = num;
                 big_bin[1] = E_bin;
                 big_bin[2] = theta_bin;
                 big_bin[3] = phi_bin;
            }
            counts[E_bin][theta_bin][phi_bin] = num;
            //if(num > 3) std::cout << "Another one! "<<num<<"\n" << std::flush;
        }
        else
        {
            if(lattice != NULL) lattice->undetectable_num++;
        }
        if(logNum > saveNum)
        {
            save();
            logNum = 0;
        }
    }
    // We use this as a after saving.
    ion.index = did_hit;
}

void Detector::finish(std::ofstream &out_file)
{
}

void Detector::start()
{
}

void SpectrumDetector::start()
{
    char buffer[1000];
    const char *format = 
    "--------------------------------------------------------\n"
    "\n"
    "SAFARI Spectra File, Ranges in this file are as follows:"
    "\n"
    "Energy range: %.2f to %.2f eV\n"
    "Theta range: %.2f to %.2f Degrees\n"
    "Phi range: %.2f to %.2f Degrees\n"
    "\n"
    "This file is split into blocks for each energy section\n"
    "Each block contains a header row and column, which states\n"
    "the theta and phi angles for those blocks respectively\n"
    "\n";

    sprintf(buffer, format, 
    e_min, settings.E0,
    (theta - dtheta), (theta + dtheta),
    (phi - dphi), (phi + dphi));

    file_header = buffer;
}

void SpectrumDetector::save()
{
    // Then save it
    mutx.lock();

    std::ofstream file;
    std::string filename = settings.output_name + ".spec";
    file.open(filename, std::ofstream::trunc);
    file << file_header;
    char buffer[2048];

    float dE = settings.E0 / ERES;
    double dT = (2 * dtheta / TRES);
    double dP = (2 * dphi / PRES);

    sprintf(buffer, 
    "Total Counts: %d\n"
    "Largest Bin: %d\n
    \t(%.2feV, Theta: %.1f Degrees, Phi: %.1f Degrees)\n"
    "\n"
    "--------------------------------------------------------\n", 
    total_counts, big_bin[0], (dE * big_bin[1] + e_min), (dT * big_bin[2] + theta - dtheta), (dP * big_bin[3] + phi - dphi));

    file << buffer;


    // We will do this the long way for now, can speed it up later I guess?
    for(int i = 0; i < ERES; i++)
    {
        file << (dE * i + e_min);
        for(int j = 0; j < TRES + 1;j++)
        {
            int t_bin = j-1;
            file << "\n\t";
            if(j != 0)
            {
                sprintf(buffer, "%.0f", (dT * t_bin + theta - dtheta));
                file << buffer;
            }
            for(int k = 0; k < PRES; k++)
            {
                file << '\t';
                // Header row
                if(j==0)
                {
                    sprintf(buffer, "%.0f", (dP * k + phi - dphi));
                    file << buffer;
                }
                else
                {
                    file << counts[i][t_bin][k];
                }
            }
        }
        file << '\n';
    }
    file << std::flush;
    mutx.unlock();
}
