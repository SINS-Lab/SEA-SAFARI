import argparse
import numpy as np
import os
import math
import matplotlib
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import scipy.constants as consts

class Traj:

    def time_unit(self):
        amu = consts.physical_constants["atomic mass unit-kilogram relationship"][0]
        eV = consts.physical_constants["atomic unit of charge"][0]
        A = consts.angstrom
        return A*math.sqrt(amu/eV)

    def load(self, filename):

        #Coordinates
        self.x = []
        self.y = []
        self.z = []

        #Momenta
        self.px = []
        self.py = []
        self.pz = []

        #Times
        self.t = []  #Total Time
        self.dt = [] #Current Time step

        #Energies
        self.T = []  #Kinetic
        self.V = []  #Potential
        self.E = []  #Total
        self.dV = [] #Potential Change

        #Counters
        self.n = []  #Time step
        self.near=[] #Number nearby

        traj_file = open(filename, 'r')
        num = 0
        for line in traj_file:
            if num == 0:
                num = num + 1
                continue
            num = num + 1
            args = line.split()
            self.x.append(float(args[0]))
            self.y.append(float(args[1]))
            self.z.append(float(args[2]))

            self.px.append(float(args[3]))
            self.py.append(float(args[4]))
            self.pz.append(float(args[5]))

            self.t.append(float(args[6]))
            self.n.append(int(args[7]))

            self.T.append(float(args[8]))
            self.V.append(float(args[9]))
            self.E.append(float(args[10]))

            self.near.append(int(args[11]))
            self.dt.append(float(args[12]))
            #self.dr_max.append(float(args[13])) We dont use this
            self.dV.append(float(args[14]))
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.z = np.array(self.z)

        self.px = np.array(self.px)
        self.py = np.array(self.py)
        self.pz = np.array(self.pz)

        self.t = np.array(self.t)
        #Convert to femtoseconds
        self.t = self.t * self.time_unit() * 1e15

        self.T = np.array(self.T)
        self.V = np.array(self.V)
        self.E = np.array(self.E)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file")
    parser.add_argument("-s", "--save", help="Whether to save the graphs", action='store_true')
    args = parser.parse_args()

    traj = Traj()
    traj.load(args.input)

    fig, ax = plt.subplots()
    ax.plot(traj.t, traj.V, label="Potential Energy")
    ax.plot(traj.t, traj.T, label="Kinetic Energy")
    ax.plot(traj.t, traj.E, label="Total Energy")

    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Energy (eV)')
    ax.legend()
    fig.show()
    if args.save:
        fig.savefig(args.input.replace('.traj', '_energy.png'))
    input("Enter to exit.")
