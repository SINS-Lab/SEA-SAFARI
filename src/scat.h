#ifndef SCAT_H_INCLUDED
#define SCAT_H_INCLUDED
#include "lattice.h"
#include <random>

extern std::default_random_engine rng;

void montecarloscat(lattice &lattice, int *num);
void gridscat(lattice &lattice, int *num);
void chainscat(lattice &lattice, int *num);

#endif // SCAT_H_INCLUDED
