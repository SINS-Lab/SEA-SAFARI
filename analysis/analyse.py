#!/usr/bin/env python3

import math
import argparse                     # Parsing command line arguments
import numpy as np                  # Array manipulation/maths
import matplotlib                   # Plotting
import os                           # Path related stuff
import scipy.signal as signal       # Peak finding
from spec_loader import Log
from spec_loader import TansTbl
from spec_loader import Spec

#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt     # Plotting

from scipy.optimize import curve_fit# Fitting the gaussians
from scipy.stats import linregress  # for R-value on the plot

def gaussian(x, a, sigma, mu):
    dx = x-mu
    dx2 = dx*dx
    s2 = sigma*sigma
    return a*np.exp(-dx2/(2*s2))

# This function is a linear + any number of gaussians
def multiples(x, *params):
    y = x * params[0] + params[1]
    for i in range(2, len(params), 3):
        a = params[i]
        sigma = params[i+1]
        mu = params[i+2]
        y = y + gaussian(x, a, sigma, mu)
    return y

def n_poly(x, *params):
    y = np.zeros_like(x)
    for i in range(len(params)):
        y = y + params[i] * np.power(x, i)
    return y

def fit_esa(values, axis, prom=0.1, minh=0, wid=2, plot=False):
    matched, properties = signal.find_peaks(values, height=prom, width=wid)

    params = []

    if len(matched) == 0:
        return False, params

    width = properties['widths']
    height = properties['peak_heights']

    max_h = np.max(values)
    min_h = np.min(values)
    peaks = []
    widths = []
    heights = []

    dx = axis[len(axis)-1] - axis[0]
    dy = values[len(axis)-1] - values[0]

    m = dy/dx
    b = values[0] - m * axis[0]

    params.append(m)
    params.append(b)

    for i in range(len(matched)):
        index = matched[i]
        h = values[index] - min_h
        if(h >= minh):
            peaks.append(axis[index])
            widths.append((axis[index]-axis[index-1])*width[i])
            heights.append(height[i])
            n = len(heights) - 1
            h = values[index]
            # We want half-width for guess at sigma
            w = widths[n] / 2
            u = axis[index]
            params.append(h)
            params.append(w)
            params.append(u)

    x_min = np.min(axis)
    x_max = np.max(axis)

    try:
        popt, pcov = curve_fit(multiples, axis, values, p0=params)
    except:
        return False, params

    x_0 = axis
    y_0 = multiples(x_0, *params)
    y_1 = multiples(x_0, *popt)

    m,b,r,p,err = linregress(values, y_1)

    fit_label = 'R={:.5f}\nLinear: {:.2e}x+{:.2f}\n'.format(r, popt[0], popt[1])
    
    for i in range(2, len(popt), 3):
        fit_label = fit_label + 'Peak: I={:.2f},E={:.2f}eV,sigma={:.2f}eV\n'.format(popt[i], popt[i+2], abs(popt[i+1]))

    if plot:
        fig,ax = plt.subplots()
        ax.plot(axis, values, label='Data')
        ax.plot(x_0, y_0, label='Initial Guess')
        ax.plot(x_0, y_1, label=fit_label)
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Intensity (Arbitrary)')
        ax.legend()
        fig.show()

    return True, params

def fit_esa_calib(values, axis, peaks, widths, heights, calib=False):
    matched, properties = signal.find_peaks(values, prominence=0.25, width=2)
    width = properties['widths']
    height = properties['prominences']
    if(len(matched)==0):
        print("Error with reading from file {}, no peaks found!".format(actualname))
        return
    if len(matched) > 1:
        print("Error with reading from file {}, more than 1 peak found!".format(actualname))
    
    # Will populate with the initial guess values
    p0 = [0,0,0]
    
    for i in range(len(matched)):
        index = matched[i]
        peaks.append(axis[index])
        widths.append((axis[index]-axis[index-1])*width[i])
        heights.append(height[i])
        h = height[i]
        w = widths[i]
        u = axis[index]
        if h > p0[0]:
            p0[0] = h
            p0[1] = w
            p0[2] = u

    x_min = np.min(axis)
    x_max = np.max(axis)

    popt, pcov = curve_fit(gaussian, axis, values, p0=p0)

    if calib:
        x_0 = axis
        y_0 = gaussian(x_0, *p0)
        y_1 = gaussian(x_0, *popt)

        fig,ax = plt.subplots()
        ax.plot(axis, values, label='Data')
        ax.plot(x_0, y_0, label='Initial Guess')
        ax.plot(x_0, y_1, label='Fit mu:{:.2f} sigma:{:.2f}'.format(popt[2], popt[1]))
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Intensity (Arbitrary)')
        ax.legend()
        fig.show()
    return popt

def interp(y_b, y_a, x_b, x_a, x):
    if y_b == y_a:
        return y_b
    dy_dx = (y_a-y_b)/(x_a-x_b)
    dy = dy_dx * (x - x_b)
    return y_b + dy

def toRange(min_e, max_e, size, vals, val_min_e, val_max_e):
    out = np.zeros(size)
    dout = (max_e - min_e) / size


    val_num = len(vals)
    dval = (val_max_e-val_min_e) / val_num

    val_index = 0
    val_index_max = val_num - 1

    for i in range(val_num):
        e_here = val_min_e + i * dval
        e_next = val_min_e + (i+1) * dval
        if e_next > min_e and val_index == 0:
            val_index = i
            if val_index > 0:
                val_index = val_index - 1
        if e_next > max_e:
            val_index_max = i
            break


    val = vals[val_index]


    prev_val = vals[val_index]
    prev_val_index = val_index

    j = 0

    for i in range(size):
        e_here = min_e + i*dout

        val_prev_e = val_min_e + prev_val_index * dval
        while val_prev_e > e_here:
            prev_val_index  = prev_val_index - 1
            val_prev_e = val_min_e + prev_val_index * dval


        val_next_e = val_min_e + val_index * dval
        while val_next_e < e_here:
            val_index  = val_index + 1
            val_next_e = val_min_e + val_index * dval

        val_index = min(val_index, val_num - 1)
        prev_val_index = min(prev_val_index, val_num - 1)
        val_index = max(val_index, 0)
        prev_val_index = max(prev_val_index, 0)

        val = vals[val_index]
        prev_val = vals[prev_val_index]

        out[i] = interp(prev_val, val, val_prev_e,\
                                       val_next_e,\
                                       e_here)
    return out

def plot(entry, directory, translation, calib=False, stacked=False, compound=False, integrated=False, esa_fit=False):
    # Check the calibration files for what the initial beam energy is.
    peaks = []
    widths = []
    heights = []
    for i in range(len(entry.calibrations)):
        filename = entry.calibrations[i]
        actualname = translation.files[filename]
        specfile = os.path.join(directory, actualname)
        spec = Spec()
        spec.load(specfile)
        if(len(spec.intensity)==0):
            spec.intensity = spec.Counts/spec.Energies
        spec.intensity = spec.intensity / np.max(spec.intensity)
        fit_esa_calib(spec.intensity, spec.Energies, peaks, widths, heights, calib=calib)

    mean = np.sum(peaks) / len(peaks)

    # This is the stacked spectra
    fig_stacked, ax_stacked = plt.subplots()
    # This is the integrated energy vs angle plot
    fig_integrated, ax_integrated = plt.subplots()
    # This is the coloured by intensity energy vs angle
    fig_compound, ax_compound = plt.subplots()

    specs = []

    # These are used for the compound plot
    thetas = []
    intensities = []
    energies = []

    # These are used for the energy vs theta total sum plots
    totals = []
    angles = []
    max_val = 1

    num_plots = len(entry.files)

    max_e = 0
    min_e = 1e20

    # Find max value for a calibration point.
    for i in range(num_plots):
        filename = entry.files[i]
        I_fc = entry.I_fc[i]
        theta_out = entry.theta_out[i]
        angles.append(theta_out)
        totals.append(np.sum(spec.Counts))
        actualname = translation.files[filename]
        specfile = os.path.join(directory, actualname)
        spec = Spec()
        specs.append(spec)
        spec.load(specfile)
        if spec.hasCounts:
            spec.Counts = spec.Counts / I_fc
            spec.intensity = spec.Counts / spec.Energies
        max_val = max(max_val, np.max(spec.intensity))

    img = np.zeros(())

    indecies = [x for x in range(num_plots)]
    sorted_indecies = [x for _,x in sorted(zip(entry.theta_out,indecies))]

    max_sig_e = 0
    min_sig_e = 1e20
    max_num = 0

    # Plot all of the spectra
    for j in range(num_plots):
        i = sorted_indecies[j]
        filename = entry.files[i]
        theta_out = entry.theta_out[i]
        spec = specs[i]

        # Add the esa fit if we want to use it
        if esa_fit:
            fit_esa(spec.intensity / np.max(spec.intensity), spec.Energies)

        # Normalize to the maximum intensity spectrum
        spec.intensity = spec.intensity / max_val

        min_k = 0
        max_k = 0

        min_test = 1e20
        max_test = 0

        for k in range(len(spec.intensity)):
            val = spec.intensity[k]
            if val > 0.01:
                test_e = spec.Energies[k]/mean
                max_sig_e = max(test_e, max_sig_e)
                min_sig_e = min(test_e, min_sig_e)

                if test_e > max_test:
                    max_test = test_e
                    max_k = k
                if test_e < min_test:
                    min_test = test_e
                    min_k = k

        max_num = max(max_num, max_k - min_k)

        # Rescale the plot so that it will fit in the range
        scale = 1
        max_i = np.max(spec.intensity)
        if max_i > 0:
            while max_i < 0.1:
                scale *= 10
                max_i *= 10
        intensity = spec.intensity * scale + i
        energy = spec.Energies/mean

        # Calculate the min/max for energy values here
        spec.max_e = np.max(energy)
        spec.min_e = np.min(energy)
        max_e = max(max_e, spec.max_e)
        min_e = min(min_e, spec.min_e)

        # Add the plot if we will be using it
        if stacked:
            label = "{}".format(theta_out)
            if scale != 1:
                label = "{} (x{})".format(theta_out, scale)
            ax_stacked.plot(energy, intensity, label=label)

    size = max_num
    img = np.zeros((size,num_plots ))
    # Fix up the e-theta plot to be constant dimensions
    for j in range(num_plots):
        i = sorted_indecies[j]
        filename = entry.files[i]
        theta_out = entry.theta_out[i]
        spec = specs[i]
        intensity = toRange(min_sig_e, max_sig_e, size, spec.intensity, spec.min_e, spec.max_e)
        for k in range(size):
            img[k][j] = intensity[k]

    title = 'Theta In: {}, E0: {}'.format(entry.theta_in, mean)

    if stacked:
        ax_stacked.set_xlim(min_e,max_e)
        ax2 = ax_stacked.twiny()
        ax2.set_xlabel('Energy (eV)')
        ax2.set_xlim(0,max_e * mean)
        ax2.tick_params(axis="x", direction="in")
        ax_stacked.legend()
        ax_stacked.set_title(title)
        ax_stacked.set_xlabel('Energy (E/E0)')
        ax_stacked.set_ylabel('Intensity (Arbitrary)')
        ax_stacked.tick_params(which='both', direction="in")
        fig_stacked.show()

    if compound:
        ax_compound.imshow(img, interpolation="bicubic", extent=(min_sig_e * mean, max_sig_e * mean, np.min(angles),np.max(angles)))
        ax_compound.invert_yaxis()
        ax_compound.set_xlabel('Outgoing angle (Degrees)')
        ax_compound.set_ylabel('Outgoing Energy (eV)')
        ax_compound.set_title(title)
        fig_compound.show()

    if integrated:
        ax_integrated.plot(angles, totals)
        ax_integrated.set_title(title)
        ax_integrated.set_xlabel('Outgoing angle (Degrees)')
        ax_integrated.set_ylabel('Intensity (Arbitrary)')
        fig_integrated.show()

def fit_img(file, emin, emax, tmin, tmax):
    img = matplotlib.image.imread(file)

    shape = img.shape
    axis = np.arange(0,shape[0])

    de = emax - emin
    dt = tmax - tmin

    de_dp = de / shape[1]
    dt_dp = dt / shape[0]

    x=[]
    y=[]
    s=[]

    imgsig = np.zeros((shape[0], shape[1]))

    max_points = 0
    img_threshold = 0.7

    print("Initializing Fits")
    for i in range(shape[0]):
        slyce = img[:,i].transpose()
        vals = slyce[0]+slyce[1]+slyce[2]
        noerr, params = fit_esa(vals, axis, prom=0.02, minh=0.1, wid=1)
        if noerr:
            max_points = max(max_points,(len(params) - 2)/3)
            for k in range(2, len(params), 3):
                x.append(i * dt_dp)
                mu = params[k+2]
                sigma = params[k+1]
                s.append(sigma * 2)
                y.append(mu * de_dp)
                if vals[mu] >= img_threshold:
                    imgsig[mu][i] = imgsig[mu][i] + 1

    lines = []
    for i in range(500):
        lines.append([[],[]])
    print("Finding Lines on Fits")

    theshold = 10
    numlines = len(lines)

    next_line = 1

    max_ds = 0

    for i in range(shape[0]):
        slyce = imgsig[:,i].transpose()
        for j in range(len(slyce)):
            if slyce[j] == 0:
                continue
            line = lines[0]
            if len(line[0]) == 0:
                lines[0][0].append(i)
                lines[0][1].append(j)
                continue
        
            line_i = 0
            min_ds = 1e5

            for k in range(numlines):
                line = lines[k]
                if len(line[1]) == 0:
                    continue
                ind = len(line[1])

                for k1 in range(ind):
                    i_s = ind-1
                    xt = line[0][i_s - k1]
                    yt = line[1][i_s - k1]
                    dy = abs(yt - j)
                    dx = abs(xt - i)
                    ds = math.sqrt(dx*dx + dy*dy)
                    if ds < min_ds:
                        line_i = k
                        min_ds = ds
            
            if min_ds >= theshold:
                line_i = next_line
                next_line = next_line + 1
                next_line = min(next_line, numlines-1)
            else:
                max_ds = max(max_ds, min_ds)

            line = lines[line_i]
            line[0].append(i)
            line[1].append(j)

    fig, ax = plt.subplots()

    s = np.array(s)

    imgplot = ax.imshow(img, interpolation="bicubic", extent=(tmin, tmax, emax, emin))
    fig.colorbar(imgplot, ax=ax)
    ax.invert_yaxis()
    ax.set_aspect(aspect=dt/de)

    index = 0
    for n in range(len(lines)):
        if len(lines[n][0]) < 10:
            continue
        # x1 = np.array(lines[n][0]) * dt_dp
        # y0 = np.array(lines[n][1]) * de_dp

        # ax.plot(x1,y0, label="{} Points".format(index))

        # popt, pcov = curve_fit(n_poly, x1, y0, p0=[1,1])
        # y1 = n_poly(x1, *popt)
        # m,b,r,p,err = linregress(y0, y1)
        # ax.plot(x1,y1, label="{} 1st order R={:.4f}, {:.3e} {:.3e}".format(index, r, *popt))

        # popt, pcov = curve_fit(n_poly, x1, y0, p0=[1,1,1])
        # y1 = n_poly(x1, *popt)
        # m,b,r,p,err = linregress(y0, y1)
        # ax.plot(x1,y1, label="{} 2nd order R={:.4f}, {:.3e} {:.3e} {:.3e}".format(index, r, *popt))

        # popt, pcov = curve_fit(n_poly, x1, y0, p0=[1,1,1,1])
        # y1 = n_poly(x1, *popt)
        # m,b,r,p,err = linregress(y0, y1)
        # ax.plot(x1,y1, label="{} 3rd order R={:.4f}".format(index, r))

        # popt, pcov = curve_fit(n_poly, x1, y0, p0=[1,1,1,1,1])
        # m,b,r,p,err = linregress(y0, y1)
        # y1 = n_poly(x1, *popt)
        # ax.plot(x1,y1, label="{} 4th order R={:.4f}".format(index, r))

        index = index + 1
    ax.scatter(x,y, c="y", s=s,label="Simulation Peaks")

    other = open('./tests/vals.tab', 'r')

    dat_x = []
    dat_y = []

    for line in other:
        vars = line.split()
        dat_x.append(float(vars[0]))
        dat_y.append(float(vars[1]))
    ax.scatter(dat_x,dat_y, c="r", s=1,label="Data Peaks")

    ax.set_xlabel('Outgoing angle (Degrees)')
    ax.set_ylabel('Outgoing Energy (E/E0)')
    ax.tick_params(direction="in", which='both')

    ax.legend()
    fig.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="Directory to run from")
    parser.add_argument("--emin", help="Starting energy")
    parser.add_argument("--emax", help="Final Energy")
    parser.add_argument("--tmin", help="Starting Angle")
    parser.add_argument("--tmax", help="Final Angle")
    args = parser.parse_args()
    directory = '.' 
    if args.directory:
        directory = args.directory


    emin = 0
    emax = 1

    tmin = 0
    tmax = 90

    if args.emin:
        emin = float(args.emin)
    if args.emax:
        emax = float(args.emax)

    if args.tmin:
        tmin = float(args.tmin)
    if args.tmax:
        tmax = float(args.tmax)

    fit_img('./tests/img.png', emin, emax, tmin, tmax)

    # filetable = TansTbl()
    # filetable.load(directory)
    # log = Log()
    # log.load(directory)

    # for entry in log.entries:
    #     plot(entry, directory, filetable, calib=False, stacked=False, compound=False, integrated=False, esa_fit=True)

    input("Enter to exit")
    