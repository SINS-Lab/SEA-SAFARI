#pragma once
#include "ion.h"
#include "vec_math.h"
#include <random>

/**
 * This is the RNG which is used for randomizing the lattice sites.
 * The intial seed for this is settings.SEED
 */ 
extern std::default_random_engine temperature_rng;

/**
 * Returns a random, floating point
 * number from 0 to 1;
 * 
 * @param rng - the rng to use to generate this number.
 */ 
double frand(std::default_random_engine &rng);

/**
 * Generates a seed for an RNG, given the double value.
 * This allows converting things like coordinates and
 * indices to seeds for things like the temperature.
 */ 
uint_fast32_t make_seed(double value);

/**
 * This function should initialize the RNG, as well as any other setup
 * needed for setting the temperatures of the lattice, it should be
 * called before any scattering is done, to ensure things are setup
 * 
 * 
 * This function initializes the std deviations in position and
 * momentum for the lattice.
 */ 
void init_temps();

/**
 * This function will set the r and p values for the site, based on
 * a thermal distribution, where the temperature is defined by
 * settings.TEMP, this function uses the site's ion index,
 * so that should probably be set before calling it.
 * 
 * @param site - the site to set to a random thermal state.
 */ 
void thermaize(Site *site);

/**
 * This adjusts the ion's energy, based on the value of settings.ESIZE
 * 
 * @param ion - the ion to adjust energy for.
 */ 
void thermaize_ion(Ion &ion);