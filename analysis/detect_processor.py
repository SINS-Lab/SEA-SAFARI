import math                                          # Used for sin/cos/etc
import numpy as np                                   # General array stuff.
import platform                                      # Linux vs Windows Checks
import os                                            # Path related stuff
import shutil                                        # Used to copy files.
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib                                    # Main plotting
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt                      # More plotting stuff
from matplotlib.patches import Circle                # Cirlces on impact plot
from matplotlib.collections import PatchCollection   # Also for the circles
import subprocess                                    # For calling XYZ processor
import time

# Used for shift-click functionality
shift_is_held = False

#Used to toggle tooltips on and off
tooltips = True

def round_n(x, n):
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    return float(as_string)

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
    f = open(file, 'r', errors='ignore')
    n = 0
    data = []
    errors = 0

    try:
        for line in f:
            arr = line.split()
            if n == 0:
                # Skip, this is header line
                pass
            else:
                if len(arr) < 10:
                #   print("Error on line {}".format(n))
                    errors = errors + 1
                    continue
                try:
                    # x, y, z_min
                    # E, Theta, Phi
                    # index, weight
                    data.append([float(arr[0]), float(arr[1]),float(arr[2]),\
                                float(arr[3]),float(arr[4]),float(arr[5]),\
                                float(arr[6]),1.0])# 1.0 was float(arr[7])
                except:
                    # print("Error on line {}".format(n))
                    # errors = errors + 1
                    pass
            n = n + 1
    except:
        print("File errored {}, {}".format(file, n))
  #  if errors!=0 :
  #      print("Total Errored Lines: {} ({}%)".format(errors, round_n(errors * 100.0/n,2)))
    return data

def getDataFile(file):
    filename = file
    if not (filename.endswith('.data') or filename.endswith('.sptr')):
        filename = file+'.data'
        if not os.path.exists(filename):
            filename = file+'.sptr'
    return filename

def load(file):
    return loadFromText(getDataFile(file))

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

def sort_z(val):
    return val[2]

def sort_basis(basis):
    basis.sort(key=sort_z)
    return basis
    
# x is an array containing the values to do the gaussian for.
def gauss(x, winv):
    return np.exp(-x*x*2.*winv*winv)*winv*0.7978845608

def gaussian_at(x, sigma, mu):
    dx = x-mu
    dx2 = dx*dx
    s2 = 2*sigma*sigma
    a = 0.7978845608 / sigma
    return np.exp(-dx2/s2) * a

def interp(y_b, y_a, x_b, x_a, x):
    if y_b == y_a:
        return y_b
    dy_dx = (y_a-y_b)/(x_a-x_b)
    dy = dy_dx * (x - x_b)
    return y_b + dy

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
        if intensity[i] <= 1e-60:
            intensity[i] = 0
            
    m = np.max(intensity)
    if m != 0:
        intensity /= m
    return intensity, m

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
        self.pics = True

        # These two are used for adjusting output file names.
        self.centre = 0
        self.width = 0

        self.tmp = []

    def clear(self):
        self.detections = np.zeros((0,8))
        
    def addDetection(self, line):
        if line[3] < 0:
        #These shouldn't be in detector
            return
        self.tmp.append(line)
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

        # # Weight the trajectories based on outoging velocity
        # eArr = self.detections[...,3]
        # # eArr = np.sqrt(eArr)
        # tfArr = np.multiply(tArr, math.pi/180)
        # tfArr = np.sin(tfArr)
        # aArr = np.multiply(eArr, tfArr)

        intensity, scale = integrate(numpoints, winv, tArr, aArr, angles)

        file_name = self.outputprefix\
                  + 'Theta-'\
                  + str(self.tmin) + '-'\
                  + str(self.tmax)+'_'\
                  + str(res)+"_"+str(scale)
        out = open(file_name+'.txt', 'w')
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
            ax.tick_params(direction="in", which='both')
            if self.plots:
                fig.show()
        #The following saves the plot as a png file
            if self.pics:
                fig.savefig(file_name+'.png')
        return angles, intensity
        
    def spectrumE(self, res, numpoints=1000, write_file=True):
    
        res = res / self.safio.E0
        step = (self.safio.E0 - self.emin)/(numpoints * self.safio.E0)
        winv = 1/res
        e_min = self.emin / self.safio.E0
        energy = np.array([(e_min + x*step) for x in range(numpoints)])
        
        eArr = self.detections[...,3]/self.safio.E0
        aArr = self.detections[...,7]
        
        #Convert the points into gaussians
        intensity, scale = integrate(numpoints, winv, eArr, aArr, energy)
        
        #Calculate the kinematic factor
        k = kinematicFactor(self.tmax, self.safio.THETA0,\
                            self.safio.MASS,self.safio.ATOMS[0][0])

        # Prefix, emin-emax_eres-theta_thetasize
        template = "{}Energy-{}-{}_{}-{}_{}"
        
        if write_file:
            file_name = template.format(self.outputprefix, self.emin, self.emax, res, self.centre, self.width)
            out = open(file_name+'.txt', 'w')
            out.write('energy\tintensity\tcounts\tk-factor\tscale\n')
            #This writes out the energy into a text file
            for i in range(numpoints):
                if i == 0:
                    out.write("{}\t{}\t{}\t{}\t{}\n".format(energy[i],\
                        intensity[i], len(aArr), k, scale))
                else:
                    out.write("{}\t{}\n".format(energy[i], intensity[i]))
            out.close()
        
        if self.plots or self.pics:
            fig, ax = plt.subplots(figsize=(8.0, 6.0))
            
            ax.plot(energy, intensity)
            #ax.set_xlim(0,self.safio.E0)
            ax.set_ylim(0,1)
            ax.set_xlim(0,1)
            
            ax2 = ax.twiny()
            ax2.set_xlim(0,self.safio.E0)
            
            identification = self.outputprefix
            fig.text(0.0, 0.975, identification, fontsize=9)

            print(identification)
            print("I_E, Detections: "+str(len(aArr)))
            
            ax.set_xlabel('Energy (E/E0)')
            ax.tick_params(direction="in", which='both')
            ax2.set_xlabel('Energy (eV)')
            ax2.tick_params(axis="x", direction="in")

            ax.set_ylabel('Intensity')

            self.added_kplot = True

            if self.plots:

                self.kplot, = ax.plot([k,k],[-1,2], label='k-Factor', c='orange')
                ax.legend(handles=[self.kplot], loc='upper left')

                def on_click(event):
                    if self.added_kplot:
                        self.kplot.remove()
                        ax.get_legend().remove()
                    else:
                        self.kplot, = ax.plot([k,k],[-1,2], label='k-Factor', c='orange')
                        ax.legend(handles=[self.kplot], loc='upper left')
                    self.added_kplot = not self.added_kplot
                    fig.canvas.draw()
                    return

                fig.canvas.mpl_connect('button_press_event', on_click)
                fig.show()
                
            #The following saves the plot as a png file
            if self.pics:
                fig.savefig(file_name+'.png')
        return energy, intensity, scale

    def run_single_shot(self, close, index, args):
        #things default nicely to py on windows, the linux machine like python3
        if platform.system() != 'Linux':
            args = args.replace('python3', 'py')
        new_input = self.safio.fileIn.replace('_mod.input', '_ss.input')

        if self.safio.load_crystal:
            old_crys = new_input.replace('_ss.input', '.crys_in')
            new_crys = old_crys.replace('.crys_in', '_ss.crys_in')
            try:
                shutil.copy(old_crys, new_crys)
            except:
                print("Error copying the crystal file over!")
                print(crys_file_in)
                print(crys_file_out)
        
        self.safio.fileIn = new_input
        self.safio.setGridScat(True)
        self.safio.NUMCHA = 1
        self.safio.XSTART = close[0]
        self.safio.YSTART = close[1]
        self.safio.Ion_Index = index
        self.safio.genInputFile(fileIn=self.safio.fileIn)
        imp_x = round(close[0],3)
        imp_y = round(close[1],3)

        input_file = self.safio.fileIn.replace('.input', '')
        output_file = self.safio.fileIn.replace('.input', '') + \
                        '{}_{}'.format(imp_x, imp_y)

        cmd = [args.format(input_file, output_file, close[0], close[1], index)]
        if platform.system() != 'Linux':
            cmd = cmd[0] #Not sure why this was needed on windows...
        subprocess.Popen(cmd, shell=True)
        
    def impactParam(self, basis=None, dx=0, dy=0):
        x = self.detections[..., 0]
        y = self.detections[..., 1]
        c = self.detections[..., 3]
        
        fig, ax = plt.subplots(figsize=(8.0, 6.0))
        patches = []
        colours = []
        
        z_threshold = -self.safio.BDIST*1e3
        
        if basis is not None:
            
            basis = sort_basis(basis)
            
            for site in basis:
                if site[2] < z_threshold:
                    continue
                colours.append(site[2])
                circle = Circle((site[0], site[1]), 1)
                patches.append(circle)

            p = PatchCollection(patches, cmap=plt.get_cmap('BuGn'))
            p.set_array(np.array(colours))
            
            #Draw the basis
            ax.add_collection(p)
        
            #Add a heightmap
            fig.colorbar(p, ax=ax)
        
        ax.set_xlim(self.safio.XSTART, self.safio.XSTOP)
        ax.set_ylim(self.safio.YSTART, self.safio.YSTOP)
        
        #Draw the points
        scat = ax.scatter(x, y, c=c, cmap=plt.get_cmap('plasma'))
        fig.colorbar(scat, ax=ax)


        tool_text = "Left Click: View Point\nDouble Left Click: Open Normal-Colored VMD\nDouble Right Click: Open Nearest-Colored VMD\nShift + Left Click: Open Velocity-Colored VMD"
        select_text = 'None Selected'
        #Add selected point label
        text_selected = fig.text(0.1, 0.95, select_text,fontsize=9)
        
        ax.set_title("Detections: "+str(len(x)))
        ax.set_xlabel('X Target (Angstroms)')
        ax.set_ylabel('Y Target (Angstroms)')

        text_tooltip = fig.text(0.6, 0.9, tool_text, fontsize=9)

        #Make the selected item indicator
        px = 0
        py = 0
        self.p, = ax.plot(px,py,'r+')

        #Add an arrow indicating the beam direction
        dir = [1,0]
        lx = (self.safio.XSTOP - self.safio.XSTART) / 10
        ly = (self.safio.YSTOP - self.safio.YSTART) / 10
        #Adds some padding to start location
        start = [self.safio.XSTART + lx/2,self.safio.YSTART + ly/2]
        dir[0] = math.cos(math.radians(self.safio.PHI0))
        dir[1] = math.sin(math.radians(self.safio.PHI0))

        if dir[0] < 0:
            #Also padd ends
            start[0] = self.safio.XSTOP - lx/2
        if dir[1] < 0:
            #Also padd ends
            start[1] = self.safio.YSTOP - ly/2
        d_arrow = [0,0]
        d_arrow[0] = dir[0] * lx
        d_arrow[1] = dir[1] * ly
        w =  max(lx, ly)/10
        hw = max(lx, ly)/5
        hl = max(lx, ly)/5

        self.arrow = ax.arrow(start[0], start[1], d_arrow[0], d_arrow[1], width=w, head_width=hw, head_length=hl, color='c', visible=len(x)>0)
        
        def onclick(event):
            if event.xdata is None or not tooltips:
                return

            close = [1e20, 1e20]
            distsq = close[0]**2 + close[1]**2
            index = -1
            ion_index = -1

            for i in range(len(x)):
                dxsq = (x[i]-event.xdata)**2
                dysq = (y[i]-event.ydata)**2
                if distsq > dxsq + dysq:
                    distsq = dxsq + dysq
                    close[0] = x[i]
                    close[1] = y[i]
                    index = i
                    ion_index = self.detections[..., 6][i]

            if event.dblclick and event.button == 1 and not shift_is_held:
                print("Setting up a safari run for a single shot")
                # Setup a single run safari for this.
                self.run_single_shot(close, ion_index,\
                                'python3 detect_impact.py -i {} -o {} -x {} -y {} -s {} -r')
            if event.dblclick and event.button == 3:
                # Setup a single run safari using nearness colored data
                print("Setting up a safari run for a nearness colored dataset")
                # Setup a single run safari for this.
                self.run_single_shot(close, ion_index,\
                                'python3 detect_impact.py -i {} -o {} -x {} -y {} -s {} -r -c nearby')
            if event.button == 1 and shift_is_held:
                # Setup a single run safari using velocity colored data
                print("Setting up a safari run for a velocity colored dataset")
                # Setup a single run safari for this.
                self.run_single_shot(close, ion_index,\
                                'python3 detect_impact.py -i {} -o {} -x {} -y {} -s {} -r -c velocity')
            
            close[0] = round(close[0], 5)
            close[1] = round(close[1], 5)
            energy = round(self.detections[index][3], 2)
            select_text = '{}, {}eV ({})'.format(close, energy, round(energy/self.safio.E0,3))
            text_selected.set_text(select_text)
            px = close[0]
            py = close[1]
            self.p.set_xdata([px])
            self.p.set_ydata([py])
            fig.canvas.draw()

        def on_key_press(event):
            if event.key == 'shift':
                global shift_is_held
                shift_is_held = True
            if event.key == 'alt':
                global tooltips
                tooltips = not tooltips
                if not tooltips:
                    text_selected.set_text('')
                    text_tooltip.set_text('')
                    self.p.set_xdata([-1e20])
                    self.p.set_ydata([0])
                    fig.canvas.draw()
                else:
                    text_selected.set_text(select_text)
                    text_tooltip.set_text(tool_text)
                    self.p.set_xdata([px])
                    self.p.set_ydata([py])
                    fig.canvas.draw()


        def on_key_release(event):
            if event.key == 'shift':
                global shift_is_held
                shift_is_held = False

        fig.canvas.mpl_connect('key_press_event', on_key_press)
        fig.canvas.mpl_connect('key_release_event', on_key_release)
        fig.canvas.mpl_connect('button_press_event', onclick)
        
        fig.show()

        if self.pics:
            file_name = self.outputprefix\
                  + 'Imapct-'\
                  + str(self.emin) + '-'\
                  + str(self.emax)+'_'\
                  + str(res)
            fig.savefig(file_name+'.png')
        
class StripeDetector(Detector):
    
    def __init__(self, theta1, theta2, phi, width):
        super().__init__()
        self.width = abs(width)
        self.tmin = min(theta1, theta2)
        self.tmax = max(theta1, theta2)
        self.phi = (phi + 360) % 360

    def isInDetector(self, theta, phi, e):
        #These are failed trajectories, shouldn't be here!
        if e < 0:
            return False
        inTheta = theta > self.tmin and theta < self.tmax
        if not inTheta:
            return False
        phi = (phi + 360) % 360
        dphi = abs(phi - self.phi)
        if dphi < self.width:
            return True
        dphi = abs(((360-phi)%360) - self.phi)
        if dphi < self.width:
            return True
        return False

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

        self.centre = theta
        self.width = size

    def isInDetector(self, theta, phi, e):
        #These are failed trajectories, shouldn't be here!
        if e < 0:
            return False
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
        self.pics = True
        self.stuck = []
        self.buried = []
        self.other_failed = []
        self.crystal = []

    def clear(self):
        self.detector = None
        self.box_emin = None
        self.safio = None
        self.stuck = []
        self.buried = []
        self.crystal = []
        self.other_failed = []

    def clean(self, detectorType=-1, emin=-1e6, emax=1e6,\
                                     phimin=-1e6, phimax=1e6, \
                                     thmin=-1e6, thmax=1e6):
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

        self.t_min = thmin
        self.t_max = thmax

        self.e_min = emin
        self.e_max = emax

        self.p_min = phimin
        self.p_max = phimax
        
        self.detector.clear()
        start = time.time()
        print("Collecting points")
        filename = getDataFile(self.safio.filename)
        f = open(filename, 'r', errors='ignore')
        print("Loading from: "+filename)
        n = 0
        tested = 0
        hit = 0
        for line in f:
            if n == 0:
                n = 1
                continue
            arr = line.split()
            # Errored line?
            if len(arr)<10:
                n = n + 1
                print("Error on line {}".format(n))
                continue
            try:
                # x, y, z_min
                # E, Theta, Phi
                # index, weight
                traj = [float(arr[0]), float(arr[1]),float(arr[2]),\
                        float(arr[3]),float(arr[4]),float(arr[5]),\
                        float(arr[6]),1.0]
                tested = tested + 1
            except:
                print("Error on line {}".format(n))
                # errors = errors + 1
                n = n + 1
                continue

            n = n + 1

            e = traj[3]
            t = traj[4]
            p = traj[5]
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
            if e < emin or e > emax\
            or t > thmax or t < thmin\
            or p > phimax or p < phimin:
                continue
            if self.detector.isInDetector(t, p, e):
                self.detector.addDetection(traj)
                hit = hit + 1
        print("Collected points, sorting now. {} out of {} were in detector".format(hit, tested))
        self.detector.detections = np.array(self.detector.tmp)
        self.tmp = []
        end = time.time()
        print("Time to process data: {:.3f}s".format(end - start))

    def plotThetaE(self):
        
        size = 1024
        img = np.zeros((size,size))
        e_max = self.e_max
        e_min = self.e_min
        t_min = self.t_min
        t_max = self.t_max

        del_e = e_max-e_min
        del_t = t_max-t_min

        de = del_e/size
        dt = del_t/size

        print("bounds: {} {} {} {}".format(e_min, e_max, t_min, t_max))
        x = 0
        for i in range(len(self.detector.detections)):
            line = self.detector.detections[i]
            e = line[3]
            t = line[4]
            val = 0
            i_e = 0
            i_t = 0
            t = t - t_min
            e = e - e_min
            if e >= 0 and e < del_e and t >= 0 and t < del_t:
                i_e = math.floor(e/de)
                i_t = math.floor(t/dt)
                img[i_e][i_t] = img[i_e][i_t] + 1
                x = x + 1

        max_intensity = np.max(img)
        while max_intensity < 100 and size > 2:
            size = int(size / 2)
            img = np.zeros((size,size))
            de = del_e/size
            dt = del_t/size
            print("bounds: {} {} {} {}".format(e_min, e_max, t_min, t_max))
            for i in range(len(self.detector.detections)):
                line = self.detector.detections[i]
                e = line[3]
                t = line[4]
                val = 0
                i_e = 0
                i_t = 0
                t = t - t_min
                e = e - e_min
                if e >= 0 and e < del_e and t >= 0 and t < del_t:
                    i_e = math.floor(e/de)
                    i_t = math.floor(t/dt)
                    img[i_e][i_t] = img[i_e][i_t] + 1
            max_intensity = np.max(img)
        
        if self.plots:
            fig, ax = plt.subplots()
            im = ax.imshow(img, interpolation="bicubic", extent=(t_min, t_max, e_max, e_min))
            ax.invert_yaxis()
            ax.set_aspect(aspect=del_t/del_e)
            fig.colorbar(im, ax=ax)
            ax.set_title("Energy vs Theta, Counts: {}, Size: {}".format(x, size))
            ax.set_xlabel('Outgoing angle (Degrees)')
            ax.set_ylabel('Outgoing Energy (eV)')
            fig.show()
            p_max = self.p_max
            p_min = self.p_min
            formatting = "{}ETheta-{}eV-{}eV_{}-{}_{}-{}_{}.png"
            file_name = formatting.format(self.detector.outputprefix, e_min, e_max, t_min, t_max, p_min, p_max, size)
            matplotlib.image.imsave(file_name, img)

    def plotPhiTheta(self):

        size = 1024

        img = np.zeros((size,size))

        p_max = self.p_max
        p_min = self.p_min
        t_min = self.t_min
        t_max = self.t_max

        del_p = p_max-p_min
        del_t = t_max-t_min
        
        dp = del_p/size
        dt = del_t/size

        print("bounds: {} {} {} {}".format(t_min, t_max, p_min, p_max))
        x = 0
        for i in range(len(self.detector.detections)):
            line = self.detector.detections[i]
            p = line[5]
            t = line[4]
            val = 0
            i_p = 0
            i_t = 0
            t = t - t_min
            p = p - p_min
            if p >= 0 and p < del_p and t >= 0 and t < del_t:
                i_p = math.floor(p/dp)
                i_t = math.floor(t/dt)
                img[i_t][i_p] = img[i_t][i_p] + 1
                x = x + 1

        max_intensity = np.max(img)
        while max_intensity < 100 and size > 2:
            size = int(size / 2)
            img = np.zeros((size,size))
            dp = del_p/size
            dt = del_t/size
            print("bounds: {} {} {} {}".format(t_min, t_max, p_min, p_max))
            for i in range(len(self.detector.detections)):
                line = self.detector.detections[i]
                p = line[5]
                t = line[4]
                val = 0
                i_p = 0
                i_t = 0
                t = t - t_min
                p = p - p_min
                if p >= 0 and p < del_p and t >= 0 and t < del_t:
                    i_p = math.floor(p/dp)
                    i_t = math.floor(t/dt)
                    img[i_t][i_p] = img[i_t][i_p] + 1
            max_intensity = np.max(img)
        
        if self.plots:
            fig, ax = plt.subplots()
            im = ax.imshow(img, interpolation="bicubic", extent=(p_max, p_min, t_min, t_max))
            ax.invert_yaxis()
            ax.set_aspect(aspect=del_p/del_t)
            fig.colorbar(im, ax=ax)
            ax.set_title("Theta vs Phi, Counts: {}, Size: {}".format(x, size))
            ax.set_xlabel('Phi Angle (Degrees)')
            ax.set_ylabel('Theta Angle (Degrees)')
            fig.show()
            formatting = "{}PTheta-{}eV-{}eV_{}-{}_{}.png"
            file_name = formatting.format(self.detector.outputprefix, p_max, p_min, t_min, t_max, size)
            matplotlib.image.imsave(file_name, img)
        
