from PyQt5.QtWidgets import QWidget, QApplication
from PyQt5.QtWidgets import QGridLayout, QHBoxLayout, QVBoxLayout, QComboBox
from PyQt5.QtWidgets import QLineEdit, QLabel, QPushButton
from scipy.io import FortranFile
import os
import math
import time
import numpy as np
import platform
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import safari_input
import subprocess
import xyz_postprocess as xyz_p

# Used for shift-click functionality
shift_is_held = False

def loadCrystal(name):
    f = open(name+'.crys', 'r')
    data = []
    for line in f:
        arr = line.split()
        # x, y, z
        # Z, mass
        data.append([float(arr[0]),float(arr[1]),float(arr[2]),\
                     float(arr[3]),float(arr[4])])
    return data

def loadFromText(file):
    f = open(file, 'r')
    n = 0
    data = []
    for line in f:
        arr = line.split()
        if n == 0:
            # Skip, this is header line
            pass
        else:
            # x, y, z_min
            # E, Theta, Phi
            # level, weight
            data.append([float(arr[0]), float(arr[1]),float(arr[2]),\
                         float(arr[3]),float(arr[4]),float(arr[5]),\
                         float(arr[6]),float(arr[7])])
        n = n + 1
    return data

def load(file):
    return loadFromText(file+'.data')

def kinematicFactor(theta_final, theta_inc, massProject, massTarget):
    mu = massProject/massTarget
    theta_tsa = 180 - theta_inc - theta_final
    cos_tsa = math.cos(math.radians(theta_tsa))
    sin_tsa = math.sin(math.radians(theta_tsa))
    k = ((mu/(1+mu))**2) * (cos_tsa + (1 / (mu**2) - sin_tsa**2)**0.5)**2
    return k

def unit(theta, phi):
    th = theta * math.pi / 180
    ph = phi * math.pi / 180
    sinth = math.sin(th)
    x = sinth * math.cos(ph)
    y = sinth * math.sin(ph)
    z = math.cos(th)
    s = math.sqrt(x*x + y*y + z*z)
    return np.array([x/s, y/s, z/s])

# x is an array containing the values to do the gaussian for.
def gauss(x, winv):
    return np.exp(-x*x*2.*winv*winv)*winv*0.7978845608

def integrate(numpoints, winv, points, areas, axis):
    # Initializing the array to 0 breaks for some reason.
    intensity = np.array([1e-60 for x in range(numpoints)])
    
    zero = np.sum(areas)==0
    
    # We vectorize the maths here, so it only needs 1 loop.
    for i in range(numpoints):
        # eArr - energy[i] is the coordinate for the gaussian
        # Intensity of gaussian at this point
        if zero:
            intensity[i] = np.sum(gauss(points - axis[i], winv))
        else:
            intensity[i] = np.sum(gauss(points - axis[i], winv) * areas)
        
        # Cull out values that dont play nicely in excel
        if intensity[i] < 1e-60:
            intensity[i] = 0
            
    m = np.max(intensity)
    if m != 0:
        intensity /= m
    return intensity

class Detector:

    def __init__(self, *args, **kwargs):
        self.detections = np.zeros((0,8))
        self.outputprefix = 'spectrum'
        self.tmax = 180
        self.tmin = -180
        self.emin = 1e20
        self.emax = -1e20
        self.safio = None
        self.plots = True
        self.pics = False

    def clear(self):
        self.detections = np.zeros((0,8))
        
    def addDetection(self, line):
        if line[3] < 0:
        #These shouldn't be in detector
            return
        self.detections = np.vstack((self.detections, line))
        e = line[3]
        if e < self.emin:
            self.emin  = e
        if e > self.emax:
            self.emax  = e
            
    def spectrumT(self, res, numpoints=1000):
        step = (self.tmax - self.tmin) / numpoints
        winv = 1/res
        angles = np.array([(self.tmin + x*step) for x in range(numpoints)])
        
        tArr = self.detections[...,4]
        aArr = self.detections[...,7]

        intensity = integrate(numpoints, winv, tArr, aArr, angles)

        out = open(self.outputprefix\
                  + 'Theta-'\
                  + str(self.tmin) + '-'\
                  + str(self.tmax)+'_'\
                  + str(res)+'.txt', 'w')
        out.write(str(len(aArr))+'\n')
        #writes the angle 
        for i in range(numpoints):
            out.write(str(angles[i])+'\t'+str(intensity[i])+'\n')
        out.close()
        
        if self.plots or self.pics:
            fig, ax = plt.subplots()
            ax.plot(angles, intensity)
            ax.set_title("Intensity vs Theta, Detections: "+str(len(aArr)))
            ax.set_xlabel('Angle (Degrees)')
            ax.set_ylabel('Intensity')
            if self.plots:
                fig.show()
        #The following saves the plot as a png file
            if self.pics:
                fig.savefig('thetaplot.png')
        return angles, intensity
        
    def spectrumE(self, res, numpoints=1000):
    
        res = res/self.safio.E0
        step = 1/numpoints
        winv = 1/res
        energy = np.array([(x*step) for x in range(numpoints)])
        
        eArr = self.detections[...,3]/self.safio.E0
        aArr = self.detections[...,7]
        
        #Convert the points into gaussians
        intensity = integrate(numpoints, winv, eArr, aArr, energy)
        
        #Calculate the kinematic factor
        k = kinematicFactor(self.tmax, self.safio.THETA0,\
                            self.safio.MASS,self.safio.ATOMS[0][0])
                            
        file = self.outputprefix\
                  + 'Energy-'\
                  + str(self.emin) + '-'\
                  + str(self.emax)+'_'\
                  + str(res)+'.txt'

        out = open(file, 'w')
        out.write('energy\tintensity\tcounts\tk-factor\n')
        #This writes out the energy into a text file
        for i in range(numpoints):
            if i == 0:
                out.write(str(energy[i])+'\t'+str(intensity[i])+'\t'+\
                          str(len(aArr))+'\t'+str(k)+'\n')
            else:
                out.write(str(energy[i])+'\t'+str(intensity[i])+'\n')
        out.close()
        
        if self.plots or self.pics:
            fig, ax = plt.subplots()
            
            ax.plot(energy, intensity)
            #ax.set_xlim(0,self.safio.E0)
            ax.set_ylim(0,1)
            ax.set_xlim(0,1)
            
            ax2 = ax.twiny()
            ax2.set_xlim(0,self.safio.E0)
            
            kplot, = ax.plot([k,k],[-1,2], label='k-Factor')
            
            ax.set_title("I_E, Detections: "+str(len(aArr)))
            
            ax.set_xlabel('Energy (E/E0)')
            ax2.set_xlabel('Energy (eV)')

            ax.set_ylabel('Intensity')
            plt.legend(handles=[kplot])
            if self.plots:
                fig.show()
            #The following saves the plot as a png file
            if self.pics:
                fig.savefig('energyplot.png')
        return energy, intensity
        
    def impactParam(self, basis=None, dx=0, dy=0):
        x = self.detections[..., 0]
        y = self.detections[..., 1]
        c = self.detections[..., 3]
        
        fig, ax = plt.subplots(figsize=(8.0, 6.0))
        patches = []
        colours = []
        
        maxX = dx
        maxY = dy
        
        if dx == 0 and dy == 2:
            maxX = np.max(x)
            maxY = np.max(y)
            
        z_threshold = -self.safio.BDIST;
        
        if basis is not None:
            minz = 1e6
            for site in basis:
                if site[2] < z_threshold:
                    continue;
                colours.append(site[2])
                circle = Circle((site[0], site[1]), 1)
                patches.append(circle)
        
        ax.set_xlim(self.safio.XSTART, self.safio.XSTOP)
        ax.set_ylim(self.safio.YSTART, self.safio.YSTOP)

        p = PatchCollection(patches, alpha=0.4)
        p.set_array(np.array(colours))
        
        #Draw the basis
        ax.add_collection(p)
        
        #Draw the points
        scat = ax.scatter(x, y, c=c)
        fig.colorbar(scat, ax=ax)
        
        #Add a heightmap
        fig.colorbar(p, ax=ax)

        #Add selected point label
        text = fig.text(0.1, 0.95, 'None Selected',fontsize=9)
        
        ax.set_title("Detections: "+str(len(x)))
        ax.set_xlabel('X Target (Angstroms)')
        ax.set_ylabel('Y Target (Angstroms)')
        fig.text(0.6, 0.9, "Left Click: View Point\nDouble Left Click: Open Normal-Colored VMD\nDouble Right Click: Open Nearest-Colored VMD\nShift + Left Click: Open Velocity-Colored VMD", fontsize=9)

        self.p, = ax.plot(0,0,'r+')

        def onclick(event):
            if event.xdata is None:
                return

            close = [1e20, 1e20]
            distsq = close[0]**2 + close[1]**2
            index = -1

            for i in range(len(x)):
                dxsq = (x[i]-event.xdata)**2
                dysq = (y[i]-event.ydata)**2
                if distsq > dxsq + dysq:
                    distsq = dxsq + dysq
                    close[0] = x[i]
                    close[1] = y[i]
                    index = i

            if event.dblclick and event.button == 1 and not shift_is_held:
                # Setup a single run safari for this.
                self.safio.fileIn = self.safio.fileIn.replace('_mod.input', '_ss.input')
                self.safio.setGridScat(True)
                self.safio.NUMCHA = 1
                self.safio.XSTART = close[0]
                self.safio.YSTART = close[1]
                self.safio.genInputFile(fileIn=self.safio.fileIn)
                
                imp_x = round(close[0],2)
                imp_y = round(close[1],2)
                input_file = self.safio.fileIn
                output_file = self.safio.fileIn.replace('.input', '') + \
                              '{},{}.xyz'.format(imp_x, imp_y)
                cmd = ['python3 detect_impact.py -i {} -o {}'.format(input_file, output_file)]
                sub = subprocess.Popen(cmd, shell=True)

            if event.dblclick and event.button == 3:
                # Setup a single run safari using nearness colored data
                print("Setting up a safari run for a nearness colored dataset")
                self.safio.fileIn = self.safio.fileIn.replace('_mod.input', '_ss.input')
                self.safio.setGridScat(True)
                self.safio.NUMCHA = 1
                self.safio.XSTART = close[0]
                self.safio.YSTART = close[1]
                self.safio.genInputFile(fileIn=self.safio.fileIn)                
                imp_x = round(close[0],2)
                imp_y = round(close[1],2)
                input_file = self.safio.fileIn
                output_file = self.safio.fileIn.replace('.input', '') + \
                              '{},{}.xyz'.format(imp_x, imp_y)
                cmd = ['python3 detect_impact.py -i {} -o {} -c nearest'.format(input_file, output_file)]
                sub = subprocess.Popen(cmd, shell=True)

            if event.button == 1 and shift_is_held:
                # Setup a single run safari using velocity colored data
                print("Setting up a safari run for a velocity colored dataset")
                self.safio.fileIn = self.safio.fileIn.replace('_mod.input', '_ss.input')
                self.safio.setGridScat(True)
                self.safio.NUMCHA = 1
                self.safio.XSTART = close[0]
                self.safio.YSTART = close[1]
                self.safio.genInputFile(fileIn=self.safio.fileIn)            
                imp_x = round(close[0],2)
                imp_y = round(close[1],2)
                input_file = self.safio.fileIn
                output_file = self.safio.fileIn.replace('.input', '') + \
                              '{},{}.xyz'.format(imp_x, imp_y)
                cmd = ['python3 detect_impact.py -i {} -o {} -c velocity'.format(input_file, output_file)]
                sub = subprocess.Popen(cmd, shell=True)

            close[0] = round(close[0], 5)
            close[1] = round(close[1], 5)
            energy = round(self.detections[index][3], 2)
            value = '{}, {}eV ({})'.format(close, energy, round(energy/self.safio.E0,3))
            text.set_text(value)
            self.p.set_xdata([close[0]])
            self.p.set_ydata([close[1]])
            fig.canvas.draw()

        def on_key_press(event):
            if event.key == 'shift':
                global shift_is_held
                shift_is_held = True

        def on_key_release(event):
            if event.key == 'shift':
                global shift_is_held
                shift_is_held = False

        fig.canvas.mpl_connect('key_press_event', on_key_press)
        fig.canvas.mpl_connect('key_release_event', on_key_release)
        fig.canvas.mpl_connect('button_press_event', onclick)
        
        fig.show()
        
class StripeDetector(Detector):
    
    def __init__(self, theta1, theta2, phi, width):
        super().__init__()
        width = abs(width)
        self.tmin = min(theta1, theta2)
        self.tmax = max(theta1, theta2)
        self.phiMax = phi + width/2
        self.phiMin = phi - width/2

    def isInDetector(self, theta, phi, e):
        #These are failed trajectories, shouldn't be here!
        if e < 0:
            return False;
        return theta > self.tmin and theta < self.tmax\
             and phi > self.phiMin and phi < self.phiMax

class SpotDetector(Detector):

    def __init__(self, theta, phi, size):
        super().__init__()
        self.theta = theta
        self.tmin = theta
        self.tmax = theta
        self.phi = phi
        self.size = size
        self.dir = unit(theta, phi)
        self.quadDots = []
        self.quadDots.append(self.dir.dot(unit(theta - size/2, phi)))
        self.quadDots.append(self.dir.dot(unit(theta + size/2, phi)))
        self.quadDots.append(self.dir.dot(unit(theta, phi - size/2)))
        self.quadDots.append(self.dir.dot(unit(theta, phi + size/2)))

    def isInDetector(self, theta, phi, e):
        #These are failed trajectories, shouldn't be here!
        if e < 0:
            return False;
        dir = unit(theta, phi)
        dotdir = dir.dot(self.dir)
        for dot in self.quadDots:
        # This would mean it is more aligned
        # to the centre than the corner is.
            if dotdir >= dot:
                return True
        return False

    def spectrum(self, res, numpoints=1000):
        return self.spectrumE(res=res, numpoints=numpoints)

class Spectrum:

    def __init__(self):
        self.detector = None
        self.box_emin = None
        self.safio = None
        self.name = None
        self.plots = True
        self.pics = False
        self.rawData = []
        self.stuck = []
        self.buried = []
        self.other_failed = []
        self.data = []
        self.crystal = []

    def clear(self):
        self.detector = None
        self.box_emin = None
        self.safio = None
        self.rawData = []
        self.stuck = []
        self.buried = []
        self.data = []
        self.crystal = []
        self.other_failed = []

    def clean(self, data, detectorType=-1, emin=-1e6, emax=1e6,\
                                           lmin=1, lmax=20, \
                                           phimin=-1e6, phimax=1e6, \
                                           thmin=-1e6, thmax=1e6):
        self.rawData = data
        self.data = []
        for i in range(0, 4):
            self.data.append(data[i])

        # If this is not the case, detector is defined elsewhere.
        if self.detector is None:
            self.detectorType = self.safio.NDTECT
            self.detectorParams = self.safio.DTECTPAR
            if self.detectorType == 1:
                self.detector = SpotDetector(self.detectorParams[0],\
                                             self.safio.PHI0,\
                                             self.detectorParams[2])
        self.detector.safio = self.safio
        self.detector.plots = self.plots
        self.detector.pics = self.pics
        self.detector.outputprefix = self.name+'_spectrum_'
        
        self.detector.emin = emin
        self.detector.emax = emax
        
        self.detector.clear()
        for i in range(len(data)):
            traj = data[i]
            e = traj[3]
            t = traj[4]
            p = traj[5]
            l = traj[6]
            # Stuck
            if e == -100:
                self.stuck.append(traj)
                continue
            if e == -200:
                self.buried.append(traj)
                continue
            if e < 0:
                self.other_failed.append(traj)
                continue
            if e < emin or e > emax or l > lmax or l < lmin\
            or t > thmax or t < thmin or p > phimax or p < phimin:
                continue
            if self.detector.isInDetector(t, p, e):
                self.detector.addDetection(traj)
            self.data.append(traj)

    def plotThetaE(self):
        
        # X Coord on graph
        x = []
        # Y Coord on graph
        y = []
        # Dot Colour, scaled by phi.
        c = []
        
        for i in range(len(self.data)):
            line = self.data[i]
            x.append(line[4])
            y.append(line[3])
            c.append(line[5])
        
       # c = np.log(c)
        
        if np.min(c) != np.max(c):
            c = c - np.min(c)
        c = c / np.max(c)
        
        colour = [(var,0,0) for var in c]
        
        if self.plots or self.pics:
            fig, ax = plt.subplots()
            ax.scatter(x, y, c=colour)
            ax.set_title("Energy vs Theta, Detections: "+str(len(x)))
            ax.set_xlabel('Angle (Degrees)')
            ax.set_ylabel('Energy (eV)')
            if self.plots:
                fig.show()
            #The following saves the plot as a png file
            if self.pics:
                fig.savefig('energythetaplot.png')

    def plotPhiTheta(self):
        
        # X Coord on graph
        x = []
        # Y Coord on graph
        y = []
        # Dot Colour, scaled by energy.
        c = []
        
        for i in range(len(self.data)):
            line = self.data[i]
            x.append(line[5])
            y.append(line[4])
            c.append(line[3])
        
        if np.min(c) != np.max(c):
            c = c - np.min(c)
        c = c / np.max(c)
        
        colour = [(var,0,0) for var in c]
        
        fig, ax = plt.subplots()
        ax.scatter(x, y, c=colour)
        ax.set_title("Theta vs Phi, Detections: "+str(len(x)))
        ax.set_xlabel('Phi Angle (Degrees)')
        ax.set_ylabel('Theta Angle (Degrees)')
        fig.show()
        #The following saves the plot as a png file
        #fig.savefig('thetaphiplot.png')
