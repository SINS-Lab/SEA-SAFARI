#pragma once
#include "ion.h"   // Ion and Lattice
#include "safio.h" // For settings
#include <fstream> // For iofstream
#include <cmath>   // For fabs

class Detector
{
public:

    bool hit(double E, double theta, double phi)
    {
        //These are originally in -180 to 180 range,
        //+360 shift allows easier comparisons without edge cases.
        double phi_test = phi + 360;
        double phi_ref = settings.PHI0 + 360;
        double diff = fabs(phi_test - phi_ref);
        return E > 0 && diff < 10;
    }

    void log(Ion &ion, Lattice &lattice, double E, double theta, double phi)
    {
        bool did_hit = true;
        //Detectors should generally be 1 degree resolution,
        //A difference of 10 is far outside the allowed range.
        if (!hit(E, theta, phi))
        {
            theta = 0;
            phi = 90;
            E = -5;
            lattice.undetectable_num++;
            did_hit = false;
        }
        if(E > 0 || settings.save_errored)
        {       
            /**
             * This uses the default saving behaviour
             */
            char buffer[200];
            //first stuff it in the buffer
            sprintf(buffer, "%f\t%f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\n",
                    ion.r_0[0], ion.r_0[1], ion.r_0[2],
                    E, theta, phi,
                    ion.index, ion.weight,
                    ion.max_n, ion.r_min, ion.steps, ion.time);
            //Then save it
            out_file << buffer;
        }
        //We use this as a after saving.
        ion.index = did_hit;
    }
};

class SpotDetector : public Detector
{
public:
    double e_min;
    double theta;
    double phi;
    double size;

    SpotDetector(double e_min_, double theta_, double phi_, double size_)
    {
        e_min = e_min_;
        theta = theta_;
        phi = phi_;
        size = size_;
    }

    bool hit(double E, double theta_in, double phi_in)
    {
        return E > e_min &&
               theta_in < theta + size && theta_in > theta - size &&
               phi_in < phi + size && phi_in > phi - size;
    }
};

extern Detector default_detector;