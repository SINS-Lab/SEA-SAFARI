#!/usr/bin/env python3

from mendeleev import element # For looking up atomic numbers, etc
import numpy as np
import argparse               # Parsing command line arguments

import matplotlib                                    # Main plotting
import matplotlib.pyplot as plt                      # More plotting stuff

# Here we have the list of currently supported potentials
from LJ import LJ
from ZBL import ZBL
from DE import DE

supported = [LJ, ZBL, DE]

pots = {}
for var in supported:
    obj = var(1,10)
    pots[obj.abbrv()] = var
    pots[obj.abbrv().lower()] = var
    pots[obj.name()] = var

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--potential", help="potential type to use")
    parser.add_argument("-A", '--A', help="Atom A")
    parser.add_argument("-B", '--B', help="Atom B")
    parser.add_argument("--r_max", help="Maximum Distance")
    parser.add_argument("--r_step", help="Step Distance")
    args = parser.parse_args()

    # Some default values
    r_max = 20 if args.r_max is None else float(args.r_max)
    r_step = 0.004 if args.r_step is None else float(args.r_step)
    r_min = r_step

    if not args.potential:
        print("Argument (-p) for potential to use is required!")
        exit()
    if not args.potential in pots:
        print("Unknown potential type {}.".format(arg.potential))
        print("Supported types are: {}".format(pots.keys()))
        exit()
    if args.A is None or args.B is None:
        print("Required arguments -A and -B for atoms to use are missing!")
        exit()

    pot = pots[args.potential](r_max, r_step)

    A = args.A
    B = args.B

    r = np.arange(r_min, r_max, r_step)
    V_r = pot.V_r(A, B, r)
    dV_dr = pot.dV_dr(A, B, r)

    fig, (ax1, ax2) = plt.subplots(1,2)

    ax1.plot(r, V_r)
    ax1.set_title("V_r")
    ax1.set_xlabel('r (Angstroms)')
    ax1.set_ylabel('V (eV)')
    ax1.tick_params(direction="in", which='both')
    y_max = pot.V_r(A, B, 0.5)
    ax1.set_xlim([0.5, 5])
    ax1.set_ylim([min(V_r)-1, y_max])

    ax2.plot(r, dV_dr)
    ax2.set_title("dV_dr")
    ax2.set_xlabel('r (Angstroms)')
    ax2.set_ylabel('dV_dr (eV/A)')
    ax2.tick_params(direction="in", which='both')
    y_max = pot.dV_dr(A, B, 0.5)
    ax2.set_xlim([0.5, 5])
    ax2.set_ylim([min(dV_dr)-1, y_max])

    fig.show()

    pot.print_pots('tests/{}_pots_{}_{}'.format(pot.abbrv(), A,B), A, B, r)

    input("Enter to exit")