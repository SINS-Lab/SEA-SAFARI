#pragma once
#include "lattice.h"
#include "detector.h"
#include "safio.h"
#include "traj.h"
#include <random>

class ScatRoutine
{
public:
    Detector* spot_detector;
    Detector* area_detector;
    /**
    * This is the main routine to call for the scattering.
    * Sub-classes of this one will use this to then call
    * their specific functions for the actual scattering.
    * Sub classes of this need to call setupDetector and
    * cleanUpDetector themselves at the end of this run
    */
    virtual void run(Lattice *lattice, int *num){}
    /**
    * This creates the relevant detectors, and also initializes them.
    */
    virtual void setupDetector()
    {
        Ion base_ion;
        base_ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, 0, 0);

        double e_min = settings.EMIN;
        double theta = settings.detect_parameters[0];
        double phi = settings.PHI0;
        double dtheta = settings.detect_parameters[1];
        double dphi = settings.detect_parameters[2];

        if(settings.main_detector)
        {
            spot_detector = new Detector();
            spot_detector->init(e_min, theta, phi, dtheta, dphi);
            spot_detector->start(base_ion);
        }
        else spot_detector = NULL;

        if(settings.main_detector)
        {
            area_detector = new SpectrumDetector();
            area_detector->init(e_min, theta, phi, dtheta, dphi);
            area_detector->start(base_ion);
        }
        else area_detector = NULL;
    }
    /**
    * This runs finish for the detectors, and also deletes them.
    */
    void cleanUpDetector()
    {
        if(spot_detector != NULL)
        {
            spot_detector->finish();
            delete spot_detector;
        }
        if(area_detector != NULL)
        {
            area_detector->finish();
            delete area_detector;
        }
    }
    /**
    * Fires the given ion at the given lattice->
    * This also initializes the ions kinetic energy, and sets the
    * impact parameters to the given coordinates. 
    * The ion is assigned an index as well.
    * 
    * @param lattice - The lattice to shoot at
    * @param ion - The ion to shoot at the lattice
    * @param x - the x-impact parameter for the ion
    * @param y - the y-impact parameter for the ion
    * @param index - the index for this ion
    * @param log - whether to log the trajectory
    * @param xyz - whether to log xyz file version of the trajectory.
    * 
    */
    bool fire(Lattice *lattice, Ion &ion, double x, double y, int index, bool log, bool xyz)
    {
        if (!lattice->mask.inside(x, y))
        {
            lattice->out_of_mask++;
            return false;
        }
        ion.set_KE(settings.E0, settings.THETA0, settings.PHI0, x, y);
        ion.index = index;
        ion.thermal_seed = index;
        traj(ion, lattice, log, xyz, spot_detector, area_detector);
        return true;
    }
};

/**
 * Fires ions randomly at the surface, rng is used to
 * determine the target location, so this results in the
 * same random locations for each unique run with the same
 * settings.SEED.
 *
 * These random locations are in the ranges specified by
 * settings.XSTART/STOP and settings.YSTART/STOP
 *
 * The number of trajectories is stuffed in *num
 */
class MonteCarlo : public ScatRoutine
{
    void run(Lattice *lattice, int *num);
    void montecarloscat(Lattice *lattice, int ionStart, int numcha, double seed);
};
/**
 * Fires ions at the surface in a grid. This grid goes from
 * xstart to xstep in steps of xstep, and y for the equivalent 
 * using ystart, etc.
 *
 * The number of trajectories is stuffed in *num
 * 
 * It will then do a sub-scattering for sections which were
 * detected, where it divides the area around each hit up
 * into 4 quadrants, and does scattering off there, it will
 * do these divisions depth times.
 * 
 * For the initial set, current_depth should be 0.
 * 
 * The "area" flag in the data will be scaled based on current_depth.
 */
class AdaptiveGrid : public ScatRoutine
{
    void run(Lattice *lattice, int *num);
    void adaptivegridscat(double xstart, double xstep, double xstop,
                      double ystart, double ystep, double ystop,
                      Lattice *lattice,
                      int max_depth, int current_depth, int *num, int *index, int iter);
};
/**
 * Fires ions at the surface in a grid. This grid goes from
 * x = settings.XSTART to settings.XSTOP in steps of
 * settings.XSTEP, and y for the equivalent using YSTART, etc.
 *
 * The number of trajectories is stuffed in *num
 */
class GridScat : public ScatRoutine
{
    void run(Lattice *lattice, int *num);
};
/**
 * Fires ions in a line, starting at XSTART, YSTART, ending at
 * XSTOP, YSTOP. It fires NUMCHA particles in this line, in even
 * spacing between particles.
 */
class ChainScat : public ScatRoutine
{
    void run(Lattice *lattice, int *num);
};
/**
 * Fires a single ion at the surface, and logs the entire trajectory to .xyz
 * and .traj files
 */
class SingleShot : public ScatRoutine
{
    void run(Lattice *lattice, int *num);
};
