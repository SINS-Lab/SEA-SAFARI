#ifndef SCAT_H_INCLUDED
#define SCAT_H_INCLUDED
#include "lattice.h"
#include <random>

extern std::default_random_engine rng;

void montecarloscat(Lattice &lattice, int *num);
void gridscat(Lattice &lattice, int *num);
void chainscat(Lattice &lattice, int *num);

#endif // SCAT_H_INCLUDED
