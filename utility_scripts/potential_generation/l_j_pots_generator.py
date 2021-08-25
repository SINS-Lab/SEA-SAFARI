#!/usr/bin/env python3

from mendeleev import element # For looking up atomic numbers, etc
import numpy as np

import matplotlib                                    # Main plotting
import matplotlib.pyplot as plt                      # More plotting stuff

r_max = 20
r_step = 0.004
r_min = r_step

# Elements to use
Z_0 = 'Au'
Z_1 = 'Au'

# This is a map of pairs to parameters
# Parameters are in order {epsilon, sigma}
# Epsilon in eV, sigma in A
LJ_Params = {}

# Here are some example parameters
LJ_Params['Au-Au'] = [0.4129, 2.646] # https://doi.org/10.1063/1.4921075
LJ_Params['Ag-Ag'] = [0.3506, 2.648] # https://doi.org/10.1063/1.4921075
LJ_Params['Au-Ag'] = [0.3848, 2.647] # https://doi.org/10.1063/1.4921075

# fundamental charge squared in our units
esqr = 14.398

def V_r(A, B, r):
    pair = A+'-'+B
    params = LJ_Params[pair]
    epsilon = params[0]
    sigma = params[1]

    A = 4 * epsilon * sigma**12
    B = 4 * epsilon * sigma**6

    return (A/(r**12)) - (B/(r**6))

def dV_dr(A, B, r):
    pair = A+'-'+B
    params = LJ_Params[pair]
    epsilon = params[0]
    sigma = params[1]

    A = -48 * epsilon * sigma**12
    B = -24 * epsilon * sigma**6

    return - (A/(r**13)) + (B/(r**7))

def print_pots(name, A, B, r):
    V_r_ = V_r(A, B, r)
    dV_dr_ = dV_dr(A, B, r)
    fmt = '{}\t{}\t{:.5e}\t{:.5e}\n'
    pots_file = open(name+'.pots', 'w')
    header = '# Generated Lennard Jones Potentials for {}-{}.\n#\n'+\
             '# Table Range in file (min, step, max): {}-{}-{}\n#\n'+\
             '# Atom1\tAtom2\tV_r\t-dV_dr\n'
    header = header.format(A, B, r_min, r_step, r_max)
    pots_file.write(header)
    for i in range(len(V_r_)):
        pots_file.write(fmt.format(A, B, V_r_[i], dV_dr_[i]))
    pots_file.close()


r = np.arange(r_min, r_max, r_step)
V_r_ = V_r(Z_0, Z_1, r)
dV_dr_ = dV_dr(Z_0, Z_1, r)

fig, (ax1, ax2) = plt.subplots(1,2)

ax1.plot(r, V_r_)
ax1.set_title("V_r")
ax1.set_xlabel('r (Angstroms)')
ax1.set_ylabel('V (eV)')
ax1.tick_params(direction="in", which='both')
# LJ potentials go too hard wall below 1
y_max = V_r(Z_0, Z_1, 1)
ax1.set_xlim([1, 5])
ax1.set_ylim([min(V_r_)-1, y_max])

ax2.plot(r, dV_dr_)
ax2.set_title("dV_dr")
ax2.set_xlabel('r (Angstroms)')
ax2.set_ylabel('dV_dr (eV/A)')
ax2.tick_params(direction="in", which='both')
# LJ potentials go too hard wall below 1
y_max = dV_dr(Z_0, Z_1, 1)
ax2.set_xlim([1, 5])
ax2.set_ylim([min(dV_dr_)-1, y_max])


fig.show()

print_pots('tests/LJ_pots_{}_{}'.format(Z_0, Z_1), Z_0, Z_1, r)

input("Enter to exit")