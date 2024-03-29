#!/usr/bin/env python3

from mendeleev import element # For looking up atomic numbers, etc
import numpy as np

from base_pot import PotProvider

# This is a map of pairs to parameters
# Parameters for the standard SAFARI Double exponential
DE_Params = {}

# Here are some example parameters
DE_Params['Au-Au'] = [44.691, 1.164, 40987.591, 4.537]
DE_Params['Na-Au'] = [4153.6, 3.625, 27017.57, 7.286]
DE_Params['Na-Cu'] = [2051.73, 3.753, 11163.2, 6.877]

class DE(PotProvider):

    def V_r(self, A, B, r):
        pair = A+'-'+B
        if not pair in DE_Params:
            pair = B+'-'+A
        if not pair in DE_Params:
            print("Unknown Pair {}".format(pair))
            exit()
        params = DE_Params[pair]
        a = params[0]
        b = params[1]
        c = params[2]
        d = params[3]
        return a * np.exp(-b * r) + c * np.exp(-d * r)

    def dV_dr(self, A, B, r):
        pair = A+'-'+B
        if not pair in DE_Params:
            pair = B+'-'+A
        if not pair in DE_Params:
            print("Unknown Pair {}".format(pair))
            exit()
        params = DE_Params[pair]
        a = params[0]
        b = params[1]
        c = params[2]
        d = params[3]
        return b * a * np.exp(-b * r) + d * c * np.exp(-d * r)

    def name(self):
        return "Double Exponential"

    def abbrv(self):
        return "DE"