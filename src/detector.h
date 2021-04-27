#pragma once
#include "ion.h"   // Ion and Lattice
#include "safio.h" // For settings
#include <fstream> // For iofstream
#include <cmath>   // For fabs

struct LogCheck
{
public:
    double E = 0;
    double theta = 0;
    double phi = 0;

    void checkConditionsForLog(Site &ion, Lattice *lattice,
                   bool stuck, bool buried, bool froze, bool off_edge, bool discont,
                   bool ignore_bounds);
};

struct Detector
{
public:
    double e_min = 0;
    double theta = 45;
    double phi = 0;
    double dtheta = 90;
    double dphi = 5;

    int modulo = 360;

    virtual void init(double e_min_, double theta_, double phi_, double dtheta_, double dphi_)
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

    virtual bool hit(double E, double theta_in, double phi_in)
    {
        // The +360 to avoid issues with negative comparisons
        double test = remainder(phi_in + 360, modulo);
        double diff_phi = fabs(test - phi);
        double diff_theta = fabs(theta_in - theta);
        bool inE = E > e_min;
        bool inTheta = diff_theta < dtheta;
        bool inPhi = diff_phi < dphi;

        // This checks for wrapping around of phi
        inPhi = inPhi || fabs(remainder(360 - test, modulo) - phi) < dphi;

        return inE && inTheta && inPhi;
    }

    virtual void log(std::ofstream &out_file, Site &ion, Lattice *lattice, 
             bool stuck, bool buried, bool froze, bool off_edge, bool discont, 
             bool ignore_bounds);

    virtual void finish(std::ofstream &out_file);
    virtual void start();
};

#define ERES 250
#define TRES 90
#define PRES 90

struct SpectrumDetector: public Detector
{
public:
    int logNum = 0;
    int saveNum = 10000;
    int total_counts = 0;
    int big_bin[4];
    int counts[ERES][TRES][PRES];

    std::string file_header;

    SpectrumDetector()
    {

    }

    void log(std::ofstream &out_file, Site &ion, Lattice *lattice, 
            bool stuck, bool buried, bool froze, bool off_edge, bool discont, 
            bool ignore_bounds);

    int getPhiBin(double theta, double phi);

    int getThetaBin(double theta, double phi);

    int getEnergyBin(double E);

    void save();

    void start();

    void finish(std::ofstream &out_file)
    {
        if(logNum != 0) save();
    }
};

extern Detector* default_detector;
extern SpectrumDetector* spec_detector;