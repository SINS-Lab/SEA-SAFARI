#!/usr/bin/env python3

from mendeleev import element # For looking up atomic numbers, etc
import numpy as np

from base_pot import PotProvider

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

class LJ(PotProvider):

    def V_r(self, A, B, r):
        pair = A+'-'+B
        if not pair in DE_Params:
            pair = B+'-'+A
        if not pair in DE_Params:
            print("Unknown Pair {}".format(pair))
            exit()
        params = LJ_Params[pair]
        epsilon = params[0]
        sigma = params[1]

        A = 4 * epsilon * sigma**12
        B = 4 * epsilon * sigma**6

        return (A/(r**12)) - (B/(r**6))

    def dV_dr(self, A, B, r):
        pair = A+'-'+B
        if not pair in DE_Params:
            pair = B+'-'+A
        if not pair in DE_Params:
            print("Unknown Pair {}".format(pair))
            exit()
        params = LJ_Params[pair]
        epsilon = params[0]
        sigma = params[1]

        A = -48 * epsilon * sigma**12
        B = -24 * epsilon * sigma**6

        return - (A/(r**13)) + (B/(r**7))

    def name(self):
        return "Lennard Jones"

    def abbrv(self):
        return "LJ"