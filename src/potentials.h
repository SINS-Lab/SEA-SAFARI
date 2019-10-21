#ifndef POTENTIALS_H_INCLUDED
#define POTENTIALS_H_INCLUDED
#include <vector>

const double eqsr = 14.398;

void init_potentials();

/**
 * Computes the potential for the
 * given radial separation
 *
 * @param r the radial separation
 * @param n the lattice index
 * @return the value of the potential
 */
double Vr_r(double r, int n);
//version of ^ that doesn't use table
double Vr_r_init(double r, int n);

/**
 * Computes the derivative of the potential
 * with respect to the position, for the given
 * radial separation.
 *
 * @param r the radial separation
 * @param n the lattice index
 * @return the derivative of the potential
 */
double dVr_dr(double r, int n);
//version of ^ that doesn't use table
double dVr_dr_init(double r, int n);

/**
 * Computes the potential for the
 * given radial separations
 *
 * @param r the radial separations
 * @param n the lattice indices
 * @return the value of the potential
 */
double Vr_r(double r[], int n[], int num);

/**
 * Computes the derivative of the potential
 * with respect to the position, for the given
 * radial separations.
 *
 * @param r the radial separations
 * @param n the lattice indices
 * @return an array of the derivatives.
 */
double * dVr_dr(double r[], int n[], int num);

/**
 * Computes the potential for the
 * image charge for the given height
 *
 * @param z the displacement from the surface
 * @param q the charge of the ion
 * @return the value of the potential
 */
double Vi_z(double z, int q);

/**
 * Computes derivative of the potential for the
 * image charge for the given height
 *
 * @param z the displacement from the surface
 * @param q the charge of the ion
 * @return the value of the potential
 */
double dVi_dz(double z, int q);

#endif // POTENTIALS_H_INCLUDED
