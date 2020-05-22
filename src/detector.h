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

    int modulo = 360;

    void init(double e_min_, double theta_, double phi_, double dtheta_, double dphi_)
    {
        e_min = e_min_;
        theta = theta_;
        dtheta = dtheta_;
        dphi = dphi_;
        if(theta - dtheta < 0)
        {
            modulo = 180;
        }
        else
        {
            modulo = 360;
        }
        phi = remainder(phi_ + 360, modulo);
    }

    bool hit(double E, double theta_in, double phi_in)
    {
        // The +360 to avoid issues with negative comparisons
        double test = remainder(phi_in + 360, modulo);
        double diff_phi = fabs(test - phi);
        double diff_theta = fabs(theta_in - theta);
        bool inE = E > e_min;
        bool inTheta = diff_theta < dtheta;
        bool inPhi = diff_phi < dphi;
        return inE && inTheta && inPhi;
    }

    void log(std::ofstream &out_file, Site &ion, Lattice *lattice, 
             bool stuck, bool buried, bool froze, bool off_edge, bool discont, 
             bool ignore_bounds);
};

extern Detector default_detector;