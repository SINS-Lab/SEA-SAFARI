import numpy as np                  # Array manipulation/maths
import matplotlib                   # Plotting
import os                           # Path related stuff
import scipy.signal as signal       # Peak finding
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt     # Plotting

def isFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

class LogEntry:
    def __init__(self):
        self.theta_in = 0
        self.files = []
        self.theta_out = []
        self.I_fc = []
        self.T_i = []
        self.T_f = []
        self.calibrations = []

    def calib(self, calibration):
        for i in range(len(calibration)-1):
            if calibration[i+1] != '' and calibration[i+1] != '\n':
                self.calibrations.append(calibration[i+1])

    def entry(self, line):
        if len(line) != 5:
            print("Error with line: "+str(line))
            return
        self.files.append('spec{}'.format(line[0]))
        self.theta_out.append(float(line[1]))
        I_fc = 0
        
        if isFloat(line[2]):
            I_fc = float(line[2])
        else:
            I_fc = self.I_fc[len(self.I_fc)-1]
        self.I_fc.append(I_fc)
        
        T_i = 0
        if isFloat(line[3]):
            T_i = float(line[3])
        T_f = 0
        if isFloat(line[4]):
            T_f = float(line[4])
        if T_i == 0:
            T_i = T_f
        if T_f == 0:
            T_f = T_i
        self.T_f.append(T_f)
        self.T_i.append(T_i)

class Log:
    def __init__(self):
        self.entries = []
    
    def load(self, directory):
        tblfile = open(os.path.join(directory, 'log.tab'), 'r')
        entry = None
        prev = None
        for line in tblfile:
            args = line.split('\t')
            # Start of an entry
            if args[0] == 'Theta_In':
                entry = LogEntry()
                entry.theta_in = float(args[1])
                self.entries.append(entry)
            # Not got an entry yet
            elif entry is None:
                continue
            elif args[0] == 'Stability':
                pass
            elif args[0] == 'Calibration':
                entry.calib(args)
            else:
                entry.entry(args)

        for entry in self.entries:
            if prev is not None and len(entry.calibrations) == 0:
                entry.calibrations = prev.calibrations
            prev = entry


# This is the translation table for the directory
# It is mostly a dictionary of names from the ones in the log
# to the ones in the directory itself.
class TansTbl:
    def load(self, directory):
        tblfile = open(os.path.join(directory, 'TRANS.TBL'), 'r')
        self.files = {}
        for line in tblfile:
            args = line.split()
            # The file contains lines of the following format:
            # F <newname>.;1 <oldname>
            newname = args[1].replace('.;1', '')
            oldname = args[2]
            self.files[oldname] = newname
        tblfile.close()

# This is an ESA Spectrum
class Spec:
    def load(self, filename):
        self.Energies = []
        self.Counts = []
        self.times = []
        self.intensity = []
        spec_file = open(filename, 'r')
        started = False
        self.hasCounts = False
        for line in spec_file:
            if not started:
                started = line.startswith('! Energy(eV)')
                continue
            args = line.split()
            if len(args) > 2:
                #These are read as floats for normalizing later
                self.Energies.append(float(args[0]))
                self.Counts.append(float(args[1]))
                self.times.append(float(args[2]))
                self.hasCounts = True
            else:# Simlulated ones don't output this info
                i = float(args[1])
                if(i>1e-5):
                    self.Energies.append(float(args[0])*254.6)
                    self.times.append(0)
                    self.intensity.append(i)
                    self.Counts.append(0)
                    self.hasCounts = False
            if len(args) > 3:
                self.intensity.append(float(args[3]))
        spec_file.close()
        self.Energies = np.array(self.Energies)
        self.Counts = np.array(self.Counts)
        self.hasCounts = np.sum(self.Counts) > 0.5
        self.times = np.array(self.times)
        if len(self.intensity) > 0:
            self.intensity = np.array(self.intensity)