#!/usr/bin/env python3

from mendeleev import element # For looking up atomic numbers, etc
import numpy as np

import matplotlib                                    # Main plotting
import matplotlib.pyplot as plt                      # More plotting stuff

r_max = 10
r_step = 0.002
r_min = r_step

# These are the atoms to generate for, by atomic number
Z_0 = 22
Z_1 = 79

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

def print_pots(name, V_r, dV_dr):
    fmt = '{}\t{}\t{:.5f}\t{:.5f}\n'
    pots_file = open(name+'.pots', 'w')
    A = element(Z_0).symbol
    B = element(Z_1).symbol
    for i in range(len(V_r)):
        pots_file.write(fmt.format(A, B, V_r[i], dV_dr[i]))
    pots_file.close()


r = np.arange(r_min, r_max, r_step)
V_r = zbl(Z_0, Z_1, r)
dV_dr = d_zbl_dr(Z_0, Z_1, r)

fig, (ax1, ax2) = plt.subplots(1,2)

ax1.plot(r, V_r)
ax1.set_title("V_r")
ax1.set_xlabel('r (Angstroms)')
ax1.set_ylabel('V (eV)')
ax1.tick_params(direction="in", which='both')

ax2.plot(r, dV_dr)
ax2.set_title("dV_dr")
ax2.set_xlabel('r (Angstroms)')
ax2.set_ylabel('dV_dr (eV/A)')
ax2.tick_params(direction="in", which='both')


fig.show()

print_pots('test_pots', V_r, dV_dr)

input("Enter to exit")