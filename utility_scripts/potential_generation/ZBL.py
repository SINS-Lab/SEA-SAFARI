#!/usr/bin/env python3

from mendeleev import element # For looking up atomic numbers, etc
import numpy as np

from base_pot import PotProvider

# Scaling factor for ZBL potential in our units
a_0 = 0.5292

# fundamental charge squared in our units
esqr = 14.398

def a_ij_func(Z_i, Z_j):
    'ZBL Length function'
    a_ij = 0.8853 * a_0 / (Z_i ** 0.23 + Z_j ** 0.23)
    return a_ij

def phi(x):
    'ZBL Screening Function'
    phi = 0.1818 * np.exp(-3.2000 * x) + \
        0.5099 * np.exp(-0.9423 * x) + \
        0.2802 * np.exp(-0.4029 * x) + \
        0.02817 * np.exp(-0.2016 * x)
    return phi

def phiprime(x):
    'Derivative of ZBL Screening Function'
    phiprime = -0.58176 * np.exp(-3.2000 * x) - \
        0.48048 * np.exp(-0.9423 * x) - \
        0.112893 * np.exp(-0.4029 * x) - \
        0.005679 * np.exp(-0.2016 * x)
    return phiprime

def zbl(Z_i, Z_j, r):
    a_ij = a_ij_func(Z_i, Z_j)
    # Scale to the length
    x = r / a_ij
    # Compute phi from scaled length
    phi_vals = phi(x)
    zz = Z_i * Z_j
    return zz * esqr * phi_vals / r

def d_zbl_dr(Z_i, Z_j, r):
    a_ij = a_ij_func(Z_i, Z_j)
    # Scale to the length
    x = r / a_ij
    # Compute phi from scaled length
    phi_vals = phi(x)
    phiprime_vals = phiprime(x)
    zz = Z_i * Z_j
    return -zz * esqr * (phiprime_vals * x - phi_vals) / (r*r)

class ZBL(PotProvider):

    def V_r(self, A, B, r):
        Z_i = element(A).atomic_number
        Z_j = element(B).atomic_number
        return zbl(Z_i, Z_j, r)

    def dV_dr(self, A, B, r):
        Z_i = element(A).atomic_number
        Z_j = element(B).atomic_number
        return d_zbl_dr(Z_i, Z_j, r)

    def name(self):
        return "ZBL"

    def abbrv(self):
        return "ZBL"