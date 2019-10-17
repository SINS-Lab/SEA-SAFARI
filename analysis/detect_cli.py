import os
import math
import time
import numpy as np
import argparse
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import safari_input
import subprocess
import xyz_postprocess as xyz_p
import detect_processor as detect

def frange(start, end, step):
    return np.arange(start, end, step)

def azimuthal_spectrum(dir, theta, size=3, emin=0, emin_rel=0):
    if dir != '.':
        dir = os.path.join('.',dir)

    output = open(os.path.join(dir, "azimuthal_spectrum_{}_{}_{}.txt".format(theta, size, emin)), 'w');

    for filename in os.listdir(dir):
        if filename.endswith('.data'):
            file = os.path.join(dir, filename)
            safio = safari_input.SafariInput(file.replace('.data', '.input'))
            print('loading: '+filename)
            data = detect.load(file)
            print('data loaded')

            # Setup the spectrum object for this file
            spectrum = detect.Spectrum()
            spectrum.plots = False
            spectrum.name = ""
            spectrum.pics = False
            spectrum.safio = safio
            spectrum.safio.DTECTPAR[0] = theta
            spectrum.detector = None
            spectrum.detector = detect.SpotDetector(theta, safio.PHI0, size)
            
            if emin_rel!=0:
                emin = emin_rel * safio.E0
            spectrum.clean(data, emin=emin)
            output.write("{}\t{}\n".format(safio.PHI0, len(spectrum.detector.detections)));
        
    output.close();
    input('Press Enter to exit')

def e_theta_loop(dir, theta1, theta2, theta_step):
    if dir != '.':
        dir = os.path.join('.',dir)

    for filename in os.listdir(dir):
        if filename.endswith('.data'):
            file = os.path.join(dir, filename)
            safio = safari_input.SafariInput(file.replace('.data', '.input'))
            print('loading data')
            data = detect.load(file)
            print('data loaded')
            fig, ax = plt.subplots()
            num = 0
            for theta in frange(theta1, theta2, theta_step):
                print('Theta: '+str(theta))
                spectrum = detect.Spectrum()
                spectrum.plots = False
                spectrum.name = filename.replace('.data','')
                spectrum.safio = safio
                spectrum.safio.DTECTPAR[0] = theta
                spectrum.detector = None
                spectrum.clean(data)
                energy, intensity = spectrum.detector.spectrum(res=spectrum.safio.ESIZE)
                intensity = intensity + num
                ax.plot(energy, intensity, label=str(theta))
                num = num + 1
                spectrum = None
            ax.legend()
            ax.set_title("Intensity vs Energy")
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('Intensity')
            fig.show()
    input('Press Enter to exit')

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", help="Directory to run from")
parser.add_argument("-s", "--size", help="Angular size of spot detector")
parser.add_argument("-t", "--theta", help="Theta angle for detector")
parser.add_argument("-e", "--emin", help="Minimum energy to consider")
parser.add_argument("-r", "--emin_rel", help="Relative Minimum energy to consider")
args = parser.parse_args()

dir = input('Input Directory: ') if not args.directory else args.directory
theta = float(input('Detector Theta: ')) if not args.theta else float(args.theta)
size = float(input('Detector Size: ')) if not args.size else float(args.size)
emin = float(input('Minimum Energy: ')) if not args.emin else float(args.emin)
emin_rel = 0 if not args.emin_rel else float(args.emin_rel)
azimuthal_spectrum(dir, theta, size, emin, emin_rel);


#theta1 = float(input('Initial Theta: '))
#theta2 = float(input('Final Theta: '))
#theta_step = float(input('Theta Step: '))
#e_theta_loop(dir, theta1, theta2, theta_step)
