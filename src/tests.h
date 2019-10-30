#pragma once
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

/**
 * Tests some RNG related things
 * 
 */ 
void test_rngs();