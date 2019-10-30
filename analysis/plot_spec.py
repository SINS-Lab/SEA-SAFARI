import argparse
import numpy as np
import os
import math
import matplotlib
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import scipy.constants as consts

class Spec:
    def load(self, filename, compare_file, e0):
        self.Energies = []
        self.Counts = []
        self.times = []
        spec_file = open(filename, 'r')
        started = False
        for line in spec_file:
            if not started:
                started = line.startswith('! Energy(eV)')
                continue
            args = line.split()
            self.Energies.append(float(args[0]))
            #These are read as floats for normalizing later
            self.Counts.append(float(args[1]))
            self.times.append(float(args[2]))
        spec_file.close()
        self.sE = None
        self.sI = None

        #Load in the comparison data.
        if compare_file is not None:
            spec_file = open(compare_file, 'r')
            self.sE = []
            self.sI = []
            started = False
            for line in spec_file:
                if not started:
                    #Skip header line
                    started = True
                    continue
                args = line.split()
                self.sE.append(float(args[0]) * float(e0))
                self.sI.append(float(args[1]))
            self.sE = np.array(self.sE)
            self.sI = np.array(self.sI)
            spec_file.close()

        self.Energies = np.array(self.Energies)
        self.Counts = np.array(self.Counts)
        self.times = np.array(self.times)

    def plot(self, save):
        fig, ax = plt.subplots()
        print("Total Counts: {}".format(np.sum(self.Counts)))
        #Detector efficiency scales with E
        intensity = self.Counts/self.Energies
        #Normalize intensity
        intensity = intensity / np.max(intensity)

        if self.sE is not None:
            ax.plot(self.Energies, intensity, label="Data")
            ax.plot(self.sE, self.sI, label="Simulation")
            ax.legend()
        else:
            ax.plot(self.Energies, intensity)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Intensity (Arbitrary)')
        fig.show()
        if save:
            fig.savefig('spec00.png')

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file")
    parser.add_argument("-c", "--compare", help="Simulation graph to compare to")
    parser.add_argument("-e", "--e0", help="Scaling factor for simulation data")
    parser.add_argument("-s", "--save", help="Whether to save the graphs", action='store_true')
    args = parser.parse_args()

    spec = Spec()
    spec.load(args.input, args.compare, args.e0)
    spec.plot(args.save)
    input("Enter to exit")