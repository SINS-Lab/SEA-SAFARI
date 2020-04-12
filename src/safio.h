#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <map>
#include "particles.h"
#include "string_utils.h"

/**
 * This is the main configuration and settings header for Sea-Safari.
 * 
 * Units used:
 *   
 *   Energy    - eV
 *   Distance  - Angstroms
 *   Mass      - AMU
 *   Angles    - Degrees
 * 
 *   Momentum  - sqrt(eV * AMU)
 *   Velocity  - sqrt(eV / AMU)
 *   Time      - Angstrom * sqrt(AMU / eV)
 *   Force     - eV / Angstrom
 */

struct Safio
{
    // Name of the safio input file (minus the extension)
    std::string input_name;

    // Name of the files to output to (minus the extension)
    std::string output_name;

    // Parameters for the beam
    // Initial Energy
    double E0;
    // Initial Theta
    double THETA0;
    // Phi Orientation of Crystal
    double PHI0;
    // Ion Mass
    double MASS;
    // Ion Atomic Symbol
    char *SYMION;

    // Ion's Atom object
    Atom ion; // Note this is not read in from the file.

    // Parameters for the detector, currently unused

    double EMIN;  // Min Energy
    double EMAX;  // Max Energy
    double ESIZE; // Energy Resolution
    double ASIZE; // Angular Resolution
    // Detector index
    int detector_type;
    // Parameters for detector
    double *detect_parameters;

    // Whether to save the errored trajectories,
    // Otherwise only the total counts of each error is saved.
    bool save_errored = true;

    // Integration parameters
    // Minimum time step
    double DELLOW;
    // Maximum time step
    double DELT0;
    // Exponent on position error check
    double error_exponent;
    double DEMIN; // unused
    // Reference value for position error check
    double error_scale;

    int NPART;     // Target interaction number
    bool RECOIL;   // Whether ions recoil on impact
    double Z1;     // Initial height of ions
    int MAX_STEPS; // Maximum integration steps before failing

    double R_MAX;      // Maximum R value for interactions
    double rr_max;     // ^ squared, this is calculated from R_MAX
    double DR_MIN_TAB; // Step for cached values

    double ZMIN;  // Start point for values caches
    double ZSTEP; // Step for cached values

    // Radial distance from the ion to search for cells of lattice sites
    int DIST_SEARCH;
    // Average change in energy over last 3 steps to consider a failure
    int FAILED_DE;

    // Flag controlling whether scat or test mode, scat mode is 666
    int SCAT_FLAG;
    // Flag controlling scat type,
    // 666 = montecarlo, 777 = gridscat 888 = chainscat
    int SCAT_TYPE;

    // This is how many AX and AY
    // to build the lattice for,
    // The "radius" of the slab

    // Radius for Y
    double RAY;
    // Radius for X
    double RAX;

    // Site Potential Values
    int binary_potential_type;
    double *binary_potential_parameters;
    // Extra stuff needed for lennard jones things
    // If this is 0, we don't use L_J, we use springs
    int lattice_potential_type = 0;
    // This tells us what index starts the L_J parameters
    int lattice_potential_start;

    // The following indicate as are named, and are set based
    // on the value of lattice_potential_type
    bool useEinsteinSprings = true; // 0
    bool useAtomSpings = false;     // 1
    bool useLennardJones = false;   // 2
    bool rigidBounds = false;       // 4

    // Image Potential Values
    int image_potential_type;
    double *image_parameters;

    // Temperature to use for the simulation
    double TEMP;
    // Seed for RNG for the simulation, same seed will result in same random numbers
    double SEED;

    // initial index of ion, setting this allows same-trajectory thermal runs.
    int ion_index;

    // Whether to use image potentials
    bool use_image;

    // Energy to consider "stuck"
    double SENRGY;
    // Distance to consider "buried"
    double BDIST;

    // Parameters for the target
    // Lattice Constants
    double AX;
    double AY;
    double AZ;

    // Basis Coordinates

    // Number in the basis
    int NBASIS;
    // Sites in the basis (need scaling by AX,AY,AZ before use)
    std::vector<Site> BASIS;

    // Surface face
    double *face;
    // If true, will load the crystal from appropriate file name.
    bool load_crystal;
    // The surface face of the loaded crystal
    double *loaded_face;

    // Basis Atoms

    // Number of types of atoms in basis
    int NTYPES;
    // Atoms in the basis
    std::vector<Atom> ATOMS;
    // Whether the lattice sites are springy
    bool CORR;

    // Values for non-einstein springs
    double ATOMK;        // spring constant between neighbours
    double max_spring_V; // Springs break if the potential is above this.
    int neighbour_count; // Number of nearest neighbours to consider

    // Number of trajectories to run.
    // If this is 1, it will run a single gridscat,
    // located at XSTART, YSTART
    int NUMCHA;
    // X-range for impact parameters
    double XSTART;
    double XSTEP;
    double XSTOP;
    // Y-range for impact parameters
    double YSTART;
    double YSTEP;
    double YSTOP;

    // These are the vertices of the lattice mask.
    // If these are null, no mask is used.
    double *x_mask_points;
    double *y_mask_points;
    // This is how many vertices we have.
    int n_x_mask = 0;
    int n_y_mask = 0;

    // Parameters for electronic frictional forces.

    // Linear Coefficient
    // Units: sqrt(eV * AMU) * Angstrom^2
    double F_a = 0;

    // Quadratic Coefficient
    // Units: AMU * Angstrom^2
    double F_b = 0;

    // These were used for corregated image potentials
    // in the old version of Safari, that hasn't
    // been implemented yet, so these are currently unused.

    int NBZ;
    double TOL;
    int NBG;
    double GTOL;
    double *GMAX;
    double *NG;
    double *NZ;
    double *ZMAX;

    // This is used as a flag for whether scattering has started,
    // for helping with debugging
    bool scat_started = false;

    /**
     * Loads input from file of the given name in the args
     * Main should auto-populate this with the correct arguments.
     * This also opens the relvant output streams,
     * for debug, data, xyz, etc.
     * 
     * Valid args:
     *      
     *      -i [inputfile] - same format that goes in safari.input, required
     *      -o [outputfile] - optional, uses input file value if not found.
     *      -t [temperature] - the temperature for this run
     *      -f - if present, will open crys and data files
     *      -s - Enables single shot mode (NUMCHA=1, SCAT_FLAG=666)
     *         -x [value] - sets x-start (if -s is present)
     *         -y [value] - sets y-start (if -s is present)
     *         -r - if present, will only output nearby sites (--restrict)
     * 
     * @param args - a map containing program arguments.
     */
    void load(std::map<std::string, ArgValue> &args);
};

// Global settings variable, we only need this loaded once anyway.
extern Safio settings;
// This is the .data file to store the outputs.
extern std::ofstream out_file;
// This is the .traj file to ion trajectory to.
extern std::ofstream traj_file;
// This is the .traj file to ion xyz to.
extern std::ofstream xyz_file;
// This is the .dbug file to store extra info.
extern std::ofstream debug_file;
// This is the .dbug file to lattice info in
extern std::ofstream crystal_file;
