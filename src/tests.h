#ifndef TESTS_H_INCLUDED
#define TESTS_H_INCLUDED
#include "lattice.h"

/**
 * Tests the speed of cached potential lookups
 * to compare it to directly invoking the function
 */
void test_cache();

/**
 * Tests the speed of copying of the entire lattice.
 * If this is too slow, then it is better to just run
 * multiple instances of safari.
 */
void test_lattice_copy(Lattice &lattice);

#endif // TESTS_H_INCLUDED