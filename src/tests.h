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
 * This tests initialization and performance
 * of the inter-lattice forces.
 */ 
void test_lattice_springs(Lattice &lattice);

/**
 * Tests some RNG related things
 * 
 */ 
void test_rngs();

/**
 * Tests that the 1d loop is converting
 * to a 3d loop in correct order.
 */ 
void test_mask(Lattice &lattice);