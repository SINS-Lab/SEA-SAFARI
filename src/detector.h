#pragma once
#include "ion.h"   // Ion and Lattice
#include "safio.h" // For settings
#include <fstream> // For iofstream
#include <cmath>   // For fabs

class Detector
{
public:
    double e_min;
    double theta;
    double phi;
    double dtheta = 90;
    double dphi = 5;

    bool hit(double E, double theta_in, double phi_in)
    {
        // These are originally in -180 to 180 range,
        // +360 shift allows easier comparisons without edge cases.
        double test = phi_in + 360;
        double ref = phi + 360;
        double diff_phi = fabs(test - ref);
        double diff_theta = fabs(theta_in - theta);
        return E > e_min && diff_phi < dphi && diff_theta < dtheta;
    }

    void log(std::ofstream &out_file, Site &ion, Lattice *lattice, 
             bool stuck, bool buried, bool froze, bool off_edge, bool discont, 
             bool ignore_bounds);
};

extern Detector default_detector;