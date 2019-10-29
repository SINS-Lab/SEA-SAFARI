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
    def load(self, filename):
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
        self.Energies = np.array(self.Energies)
        self.Counts = np.array(self.Counts)
        self.times = np.array(self.times)

    def plot(self):
        fig, ax = plt.subplots()
        print("Total Counts: {}".format(np.sum(self.Counts)))
        #Detector efficiency scales with E
        intensity = self.Counts/self.Energies
        #Normalize intensity
        intensity = intensity / np.max(intensity)
        ax.plot(self.Energies, intensity)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Intensity (Arbitrary)')
        fig.show()

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file")
    parser.add_argument("-s", "--save", help="Whether to save the graphs", action='store_true')
    args = parser.parse_args()

    spec = Spec()
    spec.load(args.input)
    spec.plot()
    input("Enter to exit")