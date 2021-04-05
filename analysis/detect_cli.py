#!/usr/bin/env python3

import os                          # Path related things
import numpy as np                 # used to make the frange
import argparse                    # Parsing arguments
from functools import cmp_to_key   # Used to sort files by phi
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib                  # Plotting
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt    # Plotting
import safari_input                # parsing the input files
import detect_processor as detect  # Main detector code

def frange(start, end, step):
    return np.arange(start, end, step)

def compare_file_name(file1, file2):
    args1 = file1.split('_')
    args2 = file2.split('_')
    o1 = 0
    o2 = 0
    try:
        o1 = float(args1[len(args1)-1])
    except:
        o1 = 2e5
    try:
        o2 = float(args2[len(args2)-1])
    except:
        o2 = 2e5
    return o1-o2

def process_from_file(filename, size):
    axis_orig = []
    areas = []
    the_file = open(filename, 'r')
    n = 0
    for line in the_file:
        n = n + 1
        if n == 1:
            continue
        var = line.split()
        axis_orig.append(float(var[0]))
        areas.append(float(var[1]))
    areas = np.array(areas)
    areas = areas / np.max(areas)
    axis_orig = np.array(axis_orig)

    for i in range(len(areas)):
        print('{}, {}'.format(axis_orig[i], areas[i]))

    process(1000, axis_orig, areas, size, filename.replace('.txt', '_proc.txt'))

def process(num, axis_orig, areas, size, filename):
    amax = np.max(axis_orig)
    amin = np.min(axis_orig)
    step = (amax - amin) / num
    axis = np.array([(amin + x*step) for x in range(num)])
    points = detect.integrate(num, 1.0/size, axis_orig, areas, axis)    
    output = open(filename, 'w')
    output.write('{}\t{}\n'.format('Phi', 'Intensity'))
    for i in range(len(points)):
        output.write('{}\t{}\n'.format(str(axis[i]),str(points[i])))

def azimuthal_scan(dir, theta, size=512, res=3, emin=750, emax=1000, norm=True):

    if dir != '.':
        dir = os.path.join('.',dir)

    file_name= "azimuthal_spectrum_{}_{}_{}".format(theta, size, emin)

    datafiles = []

    for filename in os.listdir(dir):
        if filename.endswith('.data'):
            name = filename.replace('.data','')
            args = name.split('_')
            if len(args) < 3:
                continue
            t_arg = args[len(args)-2]
            if(float(t_arg)!=theta):
                continue
            datafiles.append(name)

    datafiles.sort(key=cmp_to_key(compare_file_name))

    plot = []
    axis_orig = []

    img = np.zeros((len(datafiles),size))

    p_min = 400
    p_max = -400
    e_max = -1
    e_min = 1e20

    if os.path.isfile(os.path.join(dir, file_name+'.npy')):
        data = np.load(os.path.join(dir, file_name+'.npy'), allow_pickle=True)
        # data[0] = img
        # data[1] = [e_min, e_max, p_max, p_min]
        img = data[0]
        b = data[1]
        e_min = b[0]
        e_max = b[1]
        p_max = b[2]
        p_min = b[3]
    else:
        i = 0
        for filename in datafiles:
            file = os.path.join(dir, '{}.data'.format(filename))
            safio = safari_input.SafariInput(file.replace('.data', '.input'))
            print('loading: '+filename)
            data = detect.load(file.replace('.data',''))
            data = detect.load(file.replace('.data',''))
            print('data loaded')

            # Setup the spectrum object for this file
            spectrum = detect.Spectrum()
            spectrum.plots = False
            spectrum.name = ""
            spectrum.pics = False
            spectrum.safio = safio
            spectrum.safio.DTECTPAR[0] = theta
            spectrum.detector = None
            phi = safio.PHI0

            spectrum.detector = detect.SpotDetector(theta, phi, res)
            p_max = max(phi, p_max)
            p_min = min(phi, p_min)
            if emin_rel!=0:
                emin = emin_rel * safio.E0
            e_min = min(emin, e_min)
            e_max = max(e_max, safio.E0)
            spectrum.clean(data, emin=emin)
            axis_orig.append(phi)
            plot.append(len(spectrum.detector.detections)*1.0)
            energy, intenisty, scale = spectrum.detector.spectrumE(safio.ESIZE,size,False)
            print("Scale of {}".format(scale))
            if norm:
                scale = 1.0
            img[i] = intenisty * scale
            i = i + 1

        output = open(os.path.join(dir, file_name + "_raw.txt"), 'w')
        output.write('{}\t{}\n'.format('Phi', 'Counts'))
        for i in range(len(plot)):
            output.write('{}\t{}\n'.format(str(axis_orig[i]),str(plot[i])))
        output.close()
        data = []
        data.append(img)
        data.append([e_min, e_max, p_max, p_min])
        data = np.array(data, dtype='object', copy=False)
        np.save(os.path.join(dir, file_name), data)

    img = img / np.max(img)

    fig, ax = plt.subplots()
    im = ax.imshow(img, interpolation="bicubic", extent=(e_min, e_max, p_max, p_min))
    ax.invert_yaxis()
    del_p = (p_max - p_min) # len(datafiles) *
    del_e = (e_max - e_min) # size * (
    ax.set_aspect(aspect=del_e/del_p)
    fig.colorbar(im, ax=ax)
    ax.set_title("Energy vs Phi")
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Phi (Degrees)')
    fig.show()

    input('Press Enter to exit')
    return img

def azimuthal_spectrum(dir, theta, size=3, emin=0, emin_rel=0):
    if dir != '.':
        dir = os.path.join('.',dir)
    
    phimin = 400
    phimax = -400
    datafiles = []

    for filename in os.listdir(dir):
        if filename.endswith('.data'):
            datafiles.append(filename.replace('.data',''))

    datafiles.sort(key=cmp_to_key(compare_file_name))

    plot = []
    axis_orig = []
    areas = []

    for filename in datafiles:
        file = os.path.join(dir, '{}.data'.format(filename))
        safio = safari_input.SafariInput(file.replace('.data', '.input'))
        print('loading: '+filename)
        data = detect.load(file.replace('.data',''))
        print('data loaded')

        # Setup the spectrum object for this file
        spectrum = detect.Spectrum()
        spectrum.plots = False
        spectrum.name = ""
        spectrum.pics = False
        spectrum.safio = safio
        spectrum.safio.DTECTPAR[0] = theta
        spectrum.detector = None
        phi = safio.PHI0
        spectrum.detector = detect.SpotDetector(theta, phi, size)
        phimax = max(phi, phimax)
        phimin = min(phi, phimin)
        if emin_rel!=0:
            emin = emin_rel * safio.E0
        spectrum.clean(data, emin=emin)
        axis_orig.append(phi)
        plot.append(len(spectrum.detector.detections)*1.0)
        areas.append(0)

    output = open(os.path.join(dir, "azimuthal_spectrum_{}_{}_{}_raw.txt".format(theta, size, emin)), 'w')
    output.write('{}\t{}\n'.format('Phi', 'Counts'))
    for i in range(len(plot)):
        output.write('{}\t{}\n'.format(str(axis_orig[i]),str(plot[i])))
    
    output.close()
    input('Press Enter to exit')

def e_theta_loop(dir, theta1, theta2, theta_step):
    if dir != '.':
        dir = os.path.join('.',dir)

    for filename in os.listdir(dir):
        if filename.endswith('.data'):
            file = os.path.join(dir, filename)
            safio = safari_input.SafariInput(file.replace('.data', '.input'))
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
                spectrum.clean()
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
parser.add_argument("-f", "--filename", help="Used for post processing mode")
parser.add_argument("-s", "--size", help="Angular size of spot detector")
parser.add_argument("-t", "--theta", help="Theta angle for detector")
parser.add_argument("-e", "--emin", help="Minimum energy to consider")
parser.add_argument("-m", "--mode", help="run mode (a,p,t)")
parser.add_argument("-r", "--emin_rel", help="Relative Minimum energy to consider")
args = parser.parse_args()

size = float(input('Detector Size: ')) if not args.size else float(args.size)

mode = 'a'
if not args.mode is None:
    mode = args.mode

if mode == 'a':
    theta = float(input('Detector Theta: ')) if not args.theta else float(args.theta)
    emin = float(input('Minimum Energy: ')) if not args.emin else float(args.emin)
    emin_rel = 0 if not args.emin_rel else float(args.emin_rel)
    dir = input('Input Directory: ') if not args.directory else args.directory
    # azimuthal_spectrum(dir, theta, size, emin, emin_rel)
    azimuthal_scan(dir, theta, 512, size, emin, emin_rel)

if mode == 'p':
    filename = args.filename
    process_from_file(filename, size)

