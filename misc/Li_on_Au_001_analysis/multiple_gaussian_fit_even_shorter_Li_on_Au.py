"""
This looks at 25 degrees Li on Au. I'm looking at making this code shorter.

"""


from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
#import argparse                     # Parsing command line arguments
import numpy as np                  # Array manipulation/maths
import matplotlib.pyplot as plt     # Plotting
#import os                           # Path related stuff
import math
import scipy.signal as signal       # Peak finding


#---------------------------------------------------------------------
#--------------------open and read files-25 deg-----------------------------
#---------------------------------------------------------------------


#declare an empty variable. This creats an empty list named txt_data:
txt_data25nofric=[]
txt_read25nofric=[]
trial=[]
#I might need an indext for the text data
cnt=0

#Using the with statement to open a file
with open('Li_on_Au_001_25Deg_nofric.txt') as file:
    #need to stay indented when we are working with the file
    
    #read the data and put it into txt_data: 0 is all lines, 
    #1 and 2 are the first line
    #txt_data=file.readlines(100)
    #In order to split the data with spaces which is how the data was done do the following:
    for line in file:
        #add contents of the file into txt_data format 
        #txt_data.append(line)
        #this only reads each line of data but does not store all of the lines
        txt_read25nofric=line.split()
        #this stores all of the lines of data
        txt_data25nofric.append(txt_read25nofric)
    #end for loop (in python this done by unindenting the for loop as seen in the next line below)    
    #print(txt_data[0][2]) #this was useful for looking at data

    
#since we opened the text file with the "with" statement above we want to close it
file.close()

#small friction open and close file
#declare an empty variable. This creats an empty list named txt_data:
txt_data25smallfric=[]
txt_read25smallfric=[]
trialsm=[]
#I might need an indext for the text data
cntsm=0

#Using the with statement to open a file
with open('Li_on_Au_001_25Degsmaller_fric.txt') as filesm:

    for line in filesm:
        txt_read25smallfric=line.split()
        txt_data25smallfric.append(txt_read25smallfric)
    #end for loop (in python this done by unindenting the for loop as seen in the next line below)    
   
#since we opened the text file with the "with" statement above we want to close it
filesm.close()

#large friction open and close file
#declare an empty variable. This creats an empty list named txt_data:
txt_data25largefric=[]
txt_read25largefric=[]
triallarge=[]
#I might need an indext for the text data
cntlarge=0

#Using the with statement to open a file
with open('Li_on_Au_001_25Degw_fric.txt') as filelarge:
 
    for line in filelarge:

        txt_read25largefric=line.split()
        txt_data25largefric.append(txt_read25largefric)
    #end of for loop

#since we opened the text file with the "with" statement above we want to close it
filelarge.close()


#0.02 friction open and close file
#declare an empty variable. This creats an empty list named txt_data:
txt_data25justrightfric=[]
txt_read25justrightfric=[]
#I might need an indext for the text data
cntlarge=0

#Using the with statement to open a file
with open('Li_on_Au_001_25Deg_justrightfric.txt') as fileright:
 
    for line in fileright:

        txt_read25justrightfric=line.split()
        txt_data25justrightfric.append(txt_read25justrightfric)


    #end of for loop

#since we opened the text file with the "with" statement above we want to close it
fileright.close()


#---------------------------------------------------------------------
#--------------------open and read files------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#--------------------adjust the data------------------------------
#---------------------------------------------------------------------

a=[]
cnt=int(1)
energy_25deg_nofric=[]
intensity_25deg_nofric=[]
energy_25deg_smallfric=[]
intensity_25deg_smallfric=[]
energy_25deg_largefric=[]
intensity_25deg_largefric=[]
shiftedlarge=[]
shiftedsmall=[]
shiftedright=[]
energy_25deg_jr=[]
intensity_25deg_jr=[]

#one normally cannot simply use integers to go in increments by
#thus we will need to change the integer into a range object as show below
for ii  in range(len(txt_data25nofric)-1):
    #ii is the index of the for loop
    #a.append(ii)
    energy_25deg_nofric.append(100*float(txt_data25nofric[cnt][0]))
    intensity_25deg_nofric.append(float(txt_data25nofric[cnt][1]))
    #a saves each of the indicies in the for loop
    #if you don't do anything in the for loop it shows up as an error
    #in spyder
    #small friction
    energy_25deg_smallfric.append(100*float(txt_data25smallfric[cnt][0]))
    intensity_25deg_smallfric.append(float(txt_data25smallfric[cnt][1]))
    #large friction
    energy_25deg_largefric.append(100*float(txt_data25largefric[cnt][0]))
    intensity_25deg_largefric.append(float(txt_data25largefric[cnt][1]))
    #just right friction
    energy_25deg_jr.append(100*float(txt_data25justrightfric[cnt][0]))
    intensity_25deg_jr.append(float(txt_data25justrightfric[cnt][1]))
    #shift the intensity up by 1.01
    shiftedsmall.append(2.00+float(txt_data25smallfric[cnt][1]))
    shiftedright.append(1.00+float(txt_data25justrightfric[cnt][1]))
    shiftedlarge.append(3.00+float(txt_data25largefric[cnt][1]))
    #Now lets do a percent difference of large and no friction
    #to do this i will need to check if the values are zero..if so save as
    #zero...need if statement
    #percentlarge(100*abs(float(txt_data25nofric[cnt][1])-float(txt_data25largefric[cnt][1]))/(float(txt_data25nofric[cnt][1])+float(txt_data25largefric[cnt][1])))
    cnt=cnt+1
    #trial=(txt_data[1000][0])

#get the number of counts that hit the detector
c_sm=[] #txt_data5nofric[1][2] -this is location of counts that were detected at 5 degrees
c_sm.append(int(txt_data25nofric[1][2]))

c_none=[] 
c_none.append(int(txt_data25smallfric[1][2]))

c_large=[] 
c_large.append(int(txt_data25largefric[1][2]))

c_justright=[] 
c_justright.append(int(txt_data25justrightfric[1][2]))


#theshifted large on same plot as no friction with large friction plot
plt.plot(energy_25deg_largefric, shiftedlarge, label="large friction counts = "+str(c_large), color='black')
plt.plot(energy_25deg_largefric, intensity_25deg_nofric, label="no friction counts = "+str(c_none), color='red')
plt.plot(energy_25deg_largefric, shiftedsmall, label="small friction counts = "+str(c_sm), color='green')
plt.plot(energy_25deg_largefric, shiftedright, label="just right friction counts = "+str(c_justright), color='blue')

plt.legend()
plt.xlabel("Energy (eV)")
plt.ylabel("Intensity shift")
plt.title('This is the data from Li ions on Au (001) surface at theta0 = 25 degrees')
plt.show()
   
#Now that we have energy and intensity values 
#let's look at them and grab them
#find peaks and widths of no friction
peaksnofric, propertiesnofric = signal.find_peaks(intensity_25deg_nofric, height=0.2, width=.01)

#The width I'm looking for is under properties in widths- take that value divide by 10. we have 1000 samples and the samples 
#go up to 100 eV hence divide the width given by 10.
actualpeakswidths=0.1*propertiesnofric['widths']
h1 = propertiesnofric['peak_heights']
peak1width=actualpeakswidths[0]
peak2width=actualpeakswidths[1]
amp1=h1[0]
amp2=h1[1]
sigma1=peak1width/2.35;
sigma2=peak2width/2.35;
x1=energy_25deg_nofric[peaksnofric[0]]
x2=energy_25deg_nofric[peaksnofric[1]]


#find peaks and widths of a= .2 friction
peakssmallfric, propertiessmallfric = signal.find_peaks(intensity_25deg_smallfric, height=0.2, width=.01)

#The width I'm looking for is under properties in widths- take that value divide by 10. we have 1000 samples and the samples 
#go up to 100 eV hence divide the width given by 10.
actualpeakswidthsm=0.1*propertiessmallfric['widths']
h1s = propertiessmallfric['peak_heights']
peak1widthsm=actualpeakswidthsm[0]
peak2widthsm=actualpeakswidthsm[1]
amp1s=h1s[0]
amp2s=h1s[1]
sigma1sm=peak1widthsm/2.35;
sigma2sm=peak2widthsm/2.35;
x1s=energy_25deg_smallfric[peakssmallfric[0]]
x2s=energy_25deg_smallfric[peakssmallfric[1]]

#find peaks and widths of a= .02 friction
peaksjr, propertiesjr = signal.find_peaks(intensity_25deg_jr, height=0.2, width=.01)


actualpeakswidthsjr=0.1*propertiesjr['widths']
h1jr = propertiesjr['peak_heights']
peak1widthjr=actualpeakswidthsjr[0]
peak2widthjr=actualpeakswidthsjr[1]
amp1jr=h1jr[0]
amp2jr=h1jr[1]
sigma1jr=peak1widthjr/2.35;
sigma2jr=peak2widthjr/2.35;
x1jr=energy_25deg_jr[peaksjr[0]]
x2jr=energy_25deg_jr[peaksjr[1]]


#---------------------------------------------------------------------------
#----------------fitting gaussian-------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#----------------25 degrees-------------------------------------------
#---------------------------------------------------------------------------

# fitting with the small friction a *SEASAFARI* = 0.0
x_none25=np.array(energy_25deg_nofric)
y_none25=np.array(intensity_25deg_nofric)
n25 = len(x_none25)  

#weird other method utilized- this works
#weird formula i don't understand via stack overflow
#https://stackoverflow.com/questions/19206332/gaussian-fit-for-python
mean25 = sum(x_none25*y_none25)/n25                  
sigma25 = np.sqrt(sum(y_none25*(x_none25)**2)/n25  )

def gaus(x,a, x0, sigma):
    #this is a gaussian fit
    #here our parameters are a, x0, and sigma
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

#https://towardsdatascience.com/basic-curve-fitting-of-scientific-data-with-python-9592244a2509
#above link helps explain curve fitting

#this is for two gaussian fits
def gauss2(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return amp1*(np.exp(-(x-cen1)**2/(2*sigma1**2))) + \
            amp2*(np.exp(-(x-cen2)**2/(2*sigma2**2))) 

#now lets make some guesses
#i'm having problems because of the initial parameters.
amp1=amp1#1
cen1=x1#90.70
sigma1=sigma1#1.361702 #sigma25
amp2=amp2#.307388
cen2=x2#93.70
sigma2=sigma2#0.565745   #sigma25

popt_2gauss, pcov_2gauss = curve_fit(gauss2, x_none25, y_none25, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2])
perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
pars_1 = popt_2gauss[0:3]
pars_2 = popt_2gauss[3:6]
gauss_peak_1 = gaus(x_none25, *pars_1)
gauss_peak_2 = gaus(x_none25, *pars_2)



# This function is a linear + any number of gaussians
def multiples(x, *params):
    y = []
    for i in range(1, len(params), 3):
        a = params[i]
        sigma = params[i+1]
        x0 = params[i+2]
        y = y + gaus(x, a, x0, sigma)
    return 

def fit_esa(values, axis, prom=0.1, minh=0, wid=2, plot=False):
    #values is the intensity; axis is the energy
    matched, properties = signal.find_peaks(values, height=prom, width=wid)

    params = []
    #If there are no peaks found exit the function call
    if len(matched) == 0:
        return False, params
    #find the widths and the peak heigts
    width = properties['widths']
    height = properties['peak_heights']

    #max_h is the maximum intensity value
    #max_h = np.max(values)
    min_h = np.min(values) #this could be non-zero due to noise floor
    #preallocate what we will find for the peaks, widths, and height actual values of what 
    #we wish to fit
    peaks = []
    widths = []
    heights = []

    #if there are found peak locations go through all of the values
    for i in range(len(matched)):
        #index contains the height and width data that found a peak.
        #i represents the number of peaks found for a piece of data
        index = matched[i]
        
        #h is the peak location value
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


    try:
        popt, pcov = curve_fit(multiples, axis, values, p0=params)
    except:
        return False, params

    x_0=axis
    y_0 = multiples(x_0, *params)
    y_1 = multiples(x_0, *popt)
    
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

