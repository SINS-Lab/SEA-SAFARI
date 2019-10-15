#ifndef SCAT_H_INCLUDED
#define SCAT_H_INCLUDED
#include "lattice.h"

void montecarloscat(lattice &lattice, int *num);
void gridscat(lattice &lattice, int *num);
void chainscat(lattice &lattice, int *num);

#endif // SCAT_H_INCLUDED
