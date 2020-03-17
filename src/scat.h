#pragma once
#include "lattice.h"
#include "detector.h"
#include <random>

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
void montecarloscat(Lattice &lattice, int *num);
/**
 * Fires ions at the surface in a grid. This grid goes from
 * x = settings.XSTART to settings.XSTOP in steps of
 * settings.XSTEP, and y for the equivalent using YSTART, etc.
 *
 * The number of trajectories is stuffed in *num
 */
void gridscat(Lattice &lattice, int *num);

/**
 * Fires ions in a line, starting at XSTART, YSTART, ending at
 * XSTOP, YSTOP. It fires NUMCHA particles in this line, in even
 * spacing between particles.
 */
void chainscat(Lattice &lattice, int *num);
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
void adaptivegridscat(double xstart, double xstep, double xstop,
                      double ystart, double ystep, double ystop,
                      Lattice &lattice, Detector &detector,
                      int max_depth, int current_depth, int *num, int *index, int iter);