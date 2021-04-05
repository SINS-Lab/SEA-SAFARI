#!/usr/bin/env python3

from PyQt5.QtWidgets import QWidget, QApplication
from PyQt5.QtWidgets import QGridLayout, QHBoxLayout, QVBoxLayout, QComboBox
from PyQt5.QtWidgets import QLineEdit, QLabel, QPushButton
# ^^ Gui components
from functools import cmp_to_key        # Used to sort files in the dropdown box
import os                               # Path related stuff
import time                             # sleep for delays
import argparse                         # Parsing command line arguments
#if you utilize the following two lines you will be able to run 
#the figures in here. This requires changing the backend of the fig.show()
#for more backend choices please see https://matplotlib.org/tutorials/introductory/usage.html#what-is-a-backend
import matplotlib                       # Plotting
#Qt5Agg is the backend
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt         # Plotting
import safari_input                     # Loading input files
import detect_processor as detect       # Main detect code

class Spectrum(detect.Spectrum):

    def __init__(self):
        super().__init__()
        self.popup = None
        
    def detectorSelection(self):
        '''This makes the dropdown menu for detectors.'''
        dropdown = QComboBox()
        dropdown.addItem('Spot')
        dropdown.addItem('Stripe')
        self.detectorDropdown = dropdown
        return dropdown
    
    def getDetectorType(self):
        return self.detectorDropdown.currentText()
        
    def detectorSettings(self):
        '''When the button is pressed, this should
           Provide settings for the type of detector selected'''
        layout = QVBoxLayout()
        dtype = self.getDetectorType()
        
        self.popup2 = QWidget()
        window = self.popup2
        
        # Initialize the box settings
        if self.box_emin is None:
            self.box_emin = QLineEdit(str(self.safio.EMIN))
            self.box_emax = QLineEdit(str(self.safio.E0))
            self.box_phimin = QLineEdit(str(self.safio.PHI0))
            self.box_phimax = QLineEdit(str(self.safio.PHI0))
            self.box_thetamin = QLineEdit(str(self.safio.DTECTPAR[0]))
            self.box_thetamax = QLineEdit(str(self.safio.DTECTPAR[0]))
            self.box_eres = QLineEdit(str(self.safio.ESIZE))
            self.box_ares = QLineEdit(str(self.safio.DTECTPAR[2]))
        
        # Handle the detector types
        if dtype == 'Spot':
            layout.addWidget(QLabel('phi'))
            layout.addWidget(self.box_phimin)
            layout.addWidget(QLabel('theta'))
            layout.addWidget(self.box_thetamin)
            layout.addWidget(QLabel('resolution'))
            layout.addWidget(self.box_ares)
        elif dtype == 'Stripe':
            print(dtype)
            
            
        # Button to close the window
        close = QPushButton('Done')
        def done():
            # Update the detector
            if dtype == 'Spot':
                phi = float(self.box_phimin.displayText())
                theta = float(self.box_thetamin.displayText())
                ares = float(self.box_ares.displayText())
                self.detector = detect.SpotDetector(theta, phi, ares)
            elif dtype == 'Stripe':
                print(dtype)
            
            window.close()
        close.clicked.connect(done)
        layout.addWidget(close)
        window.setLayout(layout)
        window.show()
        return
    
    def run(self):
        self.popup = QWidget()
        window = self.popup
        layout = QVBoxLayout()
        sublayout = QHBoxLayout()
        dtectlayout = QHBoxLayout()
        ellayout = QHBoxLayout()
        anglelayout = QHBoxLayout()
        
        # Dropdown selector for detector types
        dtectlayout.addWidget(self.detectorSelection())
        
        dtectbutton = QPushButton('Detector Settings')
        def run():
            try:
                self.detectorSettings()
            except Exception as e:
                print(e)
                pass
        dtectbutton.clicked.connect(run)
        dtectlayout.addWidget(dtectbutton)
        
        
        layout.addLayout(dtectlayout)
        

        # Fields to enter values for stuff
        label = QLabel('emin')
        emin = QLineEdit(str(self.safio.EMIN))
        sublayout.addWidget(label)
        sublayout.addWidget(emin)
        ellayout.addLayout(sublayout)

        sublayout = QHBoxLayout()
        label = QLabel('emax')
        emax = QLineEdit(str(self.safio.EMAX))
        sublayout.addWidget(label)
        sublayout.addWidget(emax)
        ellayout.addLayout(sublayout)

        sublayout = QHBoxLayout()
        label = QLabel('eres')
        eres = QLineEdit(str(self.safio.ESIZE))
        sublayout.addWidget(label)
        sublayout.addWidget(eres)
        ellayout.addLayout(sublayout)

        layout.addLayout(ellayout)
        
        sublayout = QHBoxLayout()
        label = QLabel('thmin')
        thmin = QLineEdit('0')
        sublayout.addWidget(label)
        sublayout.addWidget(thmin)
        anglelayout.addLayout(sublayout)

        sublayout = QHBoxLayout()
        label = QLabel('thmax')
        thmax = QLineEdit('90')
        sublayout.addWidget(label)
        sublayout.addWidget(thmax)
        anglelayout.addLayout(sublayout)
        
        sublayout = QHBoxLayout()
        label = QLabel('phimin')
        phimin = QLineEdit(str(self.safio.PHI0 - 5))
        sublayout.addWidget(label)
        sublayout.addWidget(phimin)
        anglelayout.addLayout(sublayout)

        sublayout = QHBoxLayout()
        label = QLabel('phimax')
        phimax = QLineEdit(str(self.safio.PHI0 + 5))
        sublayout.addWidget(label)
        sublayout.addWidget(phimax)
        anglelayout.addLayout(sublayout)

        sublayout = QHBoxLayout()
        label = QLabel('ares')
        ares = QLineEdit(str(self.safio.ASIZE))
        sublayout.addWidget(label)
        sublayout.addWidget(ares)
        anglelayout.addLayout(sublayout)

        layout.addLayout(anglelayout)

        #Button to run the spectrum stuff.
        runbutton = QPushButton('I vs Energy')
        def runIE():
            _emin=float(emin.displayText())
            _emax=float(emax.displayText())
            _phimin=float(phimin.displayText())
            _phimax=float(phimax.displayText())
            _thmin=float(thmin.displayText())
            _thmax=float(thmax.displayText())

            # phi = float(self.box_phimin.displayText())
            # ares = float(self.box_ares.displayText())
            # for i in range(50):
            #     self.detector = detect.SpotDetector(i, phi, ares)

            try:
                self.clean(emin=_emin,emax=_emax,\
                        phimin=_phimin,phimax=_phimax,\
                        thmin=_thmin,thmax=_thmax)
                self.detector.spectrum(res=float(eres.displayText()))
            except Exception as e:
                print(e)
                pass

        runbutton.clicked.connect(runIE)
        layout.addWidget(runbutton)
        
        #Button to run the spectrum stuff.
        tibutton = QPushButton('I vs Theta')
        def runIT():
            try:
                print('init detector')
                detector = self.detector
                dphi =  float(phimax.displayText()) - float(phimin.displayText())
                self.detector = detect.StripeDetector(float(thmin.displayText()),\
                                               float(thmax.displayText()),\
                                               float(phimin.displayText()),\
                                               dphi)
                print('clean data')
                self.detector.safio = self.safio
                self.clean(emin=float(emin.displayText()),\
                           emax=float(emax.displayText()),\
                           phimin=float(phimin.displayText()),\
                           phimax=float(phimax.displayText()),\
                           thmin=float(thmin.displayText()),\
                           thmax=float(thmax.displayText())\
                           )
                print('make spectrum')
                self.detector.spectrumT(res=float(ares.displayText()))
                self.detector = detector
                print('done')
            except Exception as e:
                print(e)
                pass
        tibutton.clicked.connect(runIT)
        layout.addWidget(tibutton)
        
        #Button to run the spectrum stuff.
        etbutton = QPushButton('E vs Theta')
        def runET():
            try:
                self.clean(emin=float(emin.displayText()),\
                           emax=float(emax.displayText()),\
                           phimin=float(phimin.displayText()),\
                           phimax=float(phimax.displayText()),\
                           thmin=float(thmin.displayText()),\
                           thmax=float(thmax.displayText())\
                           )
                self.plotThetaE()
            except Exception as e:
                print(e)
                pass
        etbutton.clicked.connect(runET)
        layout.addWidget(etbutton)
        
        #Button to run the spectrum stuff.
        tpbutton = QPushButton('Theta vs Phi')
        def runTP():
            try:
                print('init detector')
                detector = self.detector
                dphi =  float(phimax.displayText()) - float(phimin.displayText())
                self.detector = detect.StripeDetector(float(thmin.displayText()),\
                                               float(thmax.displayText()),\
                                               float(phimin.displayText()),\
                                               dphi)
                self.detector.safio = self.safio
                self.clean(emin=float(emin.displayText()),\
                           emax=float(emax.displayText()),\
                           phimin=float(phimin.displayText()),\
                           phimax=float(phimax.displayText()),\
                           thmin=float(thmin.displayText()),\
                           thmax=float(thmax.displayText())\
                           )
                self.plotPhiTheta()
                self.detector = detector
            except Exception as e:
                print(e)
                pass
        tpbutton.clicked.connect(runTP)
        layout.addWidget(tpbutton)
        
        #Button to run the spectrum stuff.
        impbutton = QPushButton('Impact Plot')
        def runIMP():
            try:
                self.clean(emin=float(emin.displayText()),\
                           emax=float(emax.displayText()),\
                           phimin=float(phimin.displayText()),\
                           phimax=float(phimax.displayText()),\
                           thmin=float(thmin.displayText()),\
                           thmax=float(thmax.displayText())\
                           )
                self.detector.impactParam(self.crystal,
                                         self.safio.AX, 
                                         self.safio.AY)
            except Exception as e:
                print(e)
                pass
        impbutton.clicked.connect(runIMP)
        layout.addWidget(impbutton)
        
        #Button to run the spectrum stuff.
        impbutton2 = QPushButton('Impact Plot No Basis')
        def runIMP2():
            try:
                self.clean(emin=float(emin.displayText()),\
                           emax=float(emax.displayText()),\
                           phimin=float(phimin.displayText()),\
                           phimax=float(phimax.displayText()),\
                           thmin=float(thmin.displayText()),\
                           thmax=float(thmax.displayText())\
                           )
                self.detector.impactParam(None,
                                         self.safio.AX, 
                                         self.safio.AY)
            except Exception as e:
                print(e)
                pass
        impbutton2.clicked.connect(runIMP2)
        layout.addWidget(impbutton2)

        # Button to close the window
        close = QPushButton('Done')
        def done():
            window.close()
        close.clicked.connect(done)
        layout.addWidget(close)
        window.setLayout(layout)
        window.show()
        return

def addDropdownItems(directory, subdirectories, input_list):
    for filename in os.listdir(directory):
        # Load in the inputs, but not _mod or _ss as those are temporary ones.
        if filename.endswith('.input') and not filename==('safari.input')\
           and not (filename.endswith('_mod.input') or filename.endswith('_ss.input')):
            input_list.append(os.path.join(directory, filename).replace('.input', ''))
        if filename.endswith('.dbug'):
            hasInput = False
            for filename2 in os.listdir(directory):
                hasInput = hasInput or filename2.replace('.input', '') == filename.replace('.dbug', '')
            if not hasInput:
                input_list.append(os.path.join(directory, filename).replace('.dbug', ''))
        # Also add any input files in the next level down from here.
        newDir = os.path.join(directory, filename)
        if subdirectories and os.path.isdir(newDir):
            addDropdownItems(newDir, True, input_list)

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


def fileSelection(directory='.'):
    dropdown = QComboBox()
    input_list = []
    addDropdownItems(directory, True, input_list)
    input_list.sort(key=cmp_to_key(compare_file_name))
    for line in input_list:
        dropdown.addItem(line)
    return dropdown
        
def run(directory='.'):
    app = QApplication([])
    window = QWidget()
    layout = QGridLayout()
    app.spectrums = []
    sublayout = QHBoxLayout()
    label = QLabel('input file name')
    filebox = fileSelection(directory)
    
    sublayout.addWidget(label)
    sublayout.addWidget(filebox)
    
    x = 0
    y = 0

    layout.addLayout(sublayout, x, y)

    #Make a button for running 
    run = QPushButton('Spectrum')
    def push():
        spectrum = Spectrum()
        spectrum.crystal = detect.loadCrystal(filebox.currentText())
        try:
            app.spectrums.append(spectrum)
            file = filebox.currentText()
            spectrum.name = file
            if file.endswith('.data'):
                file = file.replace('.data', '.input')
            elif file.endswith('.txt'):
                file = file.replace('.txt', '.input')
            elif file.endswith('.undata'):
                file = file.replace('.undata', '.input')
            elif file.endswith('.npy'):
                file = file.replace('.npy', '.input')
            else:
                fileA = file + '.input'
                if not os.path.exists(fileA):
                    fileA = file + '.dbug'
                file = fileA
            spectrum.safio = safari_input.SafariInput(file)
            spectrum.run()
            spectrum.popup.setWindowTitle(file)
        except Exception as e:
            print(e)
            pass
    run.clicked.connect(push)
    
    box = QVBoxLayout()
    box.addWidget(run)

    #Make a button for running all of them
    run_all = QPushButton('Run All')
    def push_2():
        AllItems = [filebox.itemText(i) for i in range(filebox.count())]
        for file in AllItems:
            try:
                spectrum = Spectrum()
                app.spectrums.append(spectrum)
                spectrum.crystal = detect.loadCrystal(filebox.currentText())
                spectrum.name = file
                if file.endswith('.data'):
                    file = file.replace('.data', '.input')
                elif file.endswith('.txt'):
                    file = file.replace('.txt', '.input')
                elif file.endswith('.undata'):
                    file = file.replace('.undata', '.input')
                elif file.endswith('.npy'):
                    file = file.replace('.npy', '.input')
                else:
                    fileA = file + '.input'
                    if not os.path.exists(fileA):
                        fileA = file + '.dbug'
                    file = fileA
                spectrum.safio = safari_input.SafariInput(file)
                spectrum.run()
                spectrum.popup.setWindowTitle(file)
            except Exception as e:
                print(e)
                pass
            #Wait a second for it to start running
            time.sleep(1)

    run_all.clicked.connect(push_2)
    
    box.addWidget(run_all)

    layout.addLayout(box, x + 10, y)
    
    # Button to clear spectrums
    clear = QPushButton('Clear')
    def clear_func():
        for spectrum in app.spectrums:
            if spectrum.popup is not None:
                spectrum.popup.close()
        app.spectrums = []
    clear.clicked.connect(clear_func)
    layout.addWidget(clear)
    
    # Button to close the window
    close = QPushButton('Close')
    def done():
        plt.close("all")
        window.close()
        for spectrum in app.spectrums:
            if spectrum.popup is not None:
                spectrum.popup.close()
    close.clicked.connect(done)
    layout.addWidget(close)
    window.setLayout(layout)
    window.setWindowTitle('Detect')
    window.show()
    app.exec_()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="Directory to run from")
    args = parser.parse_args()
    directory = '.' 
    if args.directory:
        directory = args.directory
    run(directory)
