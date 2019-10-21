#ifndef SCAT_H_INCLUDED
#define SCAT_H_INCLUDED
#include "lattice.h"
#include <random>

extern std::default_random_engine rng;

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

#endif // SCAT_H_INCLUDED