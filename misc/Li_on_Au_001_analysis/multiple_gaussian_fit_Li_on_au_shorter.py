 # -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 13:56:53 2020

@author: dcmccal
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 13:53:11 2020

@author: dave
"""


# -*- coding: utf-8 -*-
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

#import scipy.signal as signal       # Peak finding
#from scipy.optimize import curve_fit# Fitting the gaussians
#from scipy.stats import linregress  # for R-value on the plot
#I believe the three comments below only give you the data/
#from spec_loader import Log
#from spec_loader import TansTbl
#from spec_loader import Spec




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
    """
        After some finagling I figured out the data breakdown:
            txt_data is size 1001. All of this content comes out as 
            strings so we will need to properly float data later.
            txt_data[0][0] = "energy"
            txt_data[0][1] = "intensity"
            txt_data[0][2] = "counts"
            txt_data[0][3] = "k-factor"
            [0][4] = out of bounds
            txt_data[1][0] = 0 #this is the first value of "energy" (E/E0)
            txt_data[1][1] = 0 #this is the first value of intensity
            txt_data[1][2] = 3 #this is the total number of counts
            txt_data[1][3] = 0.7318521205746247  #This is the k-factor (kinematic factor)

            txt_data[2][0] = 0.001 #this is the second value of "energy" (E/E0)
            txt_data[2][1] = 0.0 #this is the second value of intensity
            
                This is how to convert string of number to actual number
                trial=float(txt_data[1][3])
            python indices start at 0 so if you see a size of 1001
            the final index value will be 1000 (rule of thumb final index = size-1)
            
            What we need here is only the data of the energy and intensity values
            not the titles of the data, nor the counts (specifcallly) nor the k-factor...for now
            
    """
    
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

"""
I was having a problem with reading and adjusting the data separately 
so I'm working on doing it all at once. Since all of the data is taken
at the same resolution this is fine
"""
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
   

"""
# Fit the dummy exponential data
pars, cov = curve_fit(f=exponential, xdata=x_dummy, ydata=y_dummy, p0=[0, 0], bounds=(-np.inf, np.inf))

f — function used for fitting (in this case exponential)
xdata — array of x-data for fitting
ydata — array of y-data for fitting
p0 — array of initial guesses for the fitting parameters (both a and b as 0)
bounds — bounds for the parameters (-∞ to ∞)
Outputs
pars — array of parameters from fit (in this case [a, b])
cov — the estimated covariance of pars which can be used to determine the standard deviations of the fitting parameters (square roots of the diagonals)
"""

#Now that we have energy and intensity values 
#let's look at them and grab them

'''
    peaks, properties = signal.find_peaks(intensity values, height=number, width=number)


#Peaks tells you the index in intensity_25-deg-nofric of where the peaks are.

https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html

properties contain

‘peak_heights’
If height is given, the height of each peak in x.

‘left_thresholds’, ‘right_thresholds’
If threshold is given, these keys contain a peaks vertical distance to its neighbouring samples.

‘prominences’, ‘right_bases’, ‘left_bases’
If prominence is given, these keys are accessible. See peak_prominences for a description of their content.

‘width_heights’, ‘left_ips’, ‘right_ips’
If width is given, these keys are accessible. See peak_widths for a description of their content.

‘plateau_sizes’, left_edges’, ‘right_edges’
If plateau_size is given, these keys are accessible and contain the indices of a peak’s edges (edges are still part of the plateau) and the calculated plateau sizes.

widths
The widths for each peak in samples - samples = length of data points.

width_heightsndarray
The height of the contour lines at which the widths where evaluated.

left_ips, right_ipsndarray
Interpolated positions of left and right intersection points of a horizontal line at the respective evaluation height.
'''

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

 #just right friction
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


 #large friction
#find peaks and widths of a= 2 friction - only has one peak
peakslg, propertieslg = signal.find_peaks(intensity_25deg_largefric, height=0.2, width=.01)


actualpeakswidthslg=0.1*propertieslg['widths']
h1lg = propertieslg['peak_heights']
peak1widthlg=actualpeakswidthslg[0]
amp1lg=h1lg[0]
sigma1lg=peak1widthlg/2.35;
x1lg=energy_25deg_largefric[peakslg[0]]


#---------------------------------------------------------------------
#--------------------adjust the data------------------------------
#---------------------------------------------------------------------



#---------------------------------------------------------------------------
#----------------fitting gaussian-------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#----------------25 degrees-------------------------------------------
#---------------------------------------------------------------------------

"""


n25 = len(x_none25)                          #length/number of data elements
#since the different energies have different intensities multiplying each energy value
#by its intensity then summing up all the energy values and dividing by the number
#of data elements gives us the average/mean:
mean25 = sum(y_none25)/n25                   #note this correction
#here the standard deviation does not occur with respect to each data point.
#I need to update my definition of sigma to reflect the correct sigma
#https://en.wikipedia.org/wiki/Standard_deviation
sigstart25=[]
simga25=[]
dx=[]
dx2=[]

for ii  in range(len(x_none25)):
    dx=(y_none25[ii]-mean25)    
    dx2=dx*dx
    sigstart25.append(dx2)
    
#    https://www.tutorialspoint.com/python/number_pow.htm
#    math.pow(x,y)=x^y


sigma25 = math.sqrt((sum(sigstart25)/n25))        #note this correction
"""

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


gauss_peak_total=gauss_peak_1+gauss_peak_2

#plot the fits...only problem is it is not finding the peaks correctly
plt.plot(x_none25, gauss_peak_1, "g")
#plt.fill_between(x_none25, gauss_peak_1.min(), gauss_peak_1, facecolor="green", alpha=0.5)
plt.plot(x_none25, gauss_peak_2, "y")
#plt.fill_between(x_none25, gauss_peak_2.min(), gauss_peak_2, facecolor="yellow", alpha=0.5)  
plt.plot(x_none25, gauss_peak_total, "m")
plt.fill_between(x_none25, gauss_peak_total.min(), gauss_peak_total, facecolor="magenta", alpha=0.5)  
plt.plot(x_none25, y_none25, "k+")
plt.xlim(85,105)
plt.title('Fig. 3 - no friction a=0 ')

plt.show()


#Now let's look at small friction
popt_sm, pcov_sm = curve_fit(gauss2, energy_25deg_smallfric, intensity_25deg_smallfric, p0=[amp1s, x1s, sigma1sm, amp2s, x2s, sigma2sm])
perr_sm = np.sqrt(np.diag(pcov_sm))
pars_1sm = popt_sm[0:3]
pars_2sm = popt_sm[3:6]
gauss_peak_1sm = gaus(energy_25deg_smallfric, *pars_1sm)
gauss_peak_2sm = gaus(energy_25deg_smallfric, *pars_2sm)


gauss_peak_totalsm=gauss_peak_1sm+gauss_peak_2sm

#plot the fits...only problem is it is not finding the peaks correctly
plt.plot(energy_25deg_smallfric, gauss_peak_1sm, "g")
#plt.fill_between(x_none25, gauss_peak_1.min(), gauss_peak_1, facecolor="green", alpha=0.5)
plt.plot(energy_25deg_smallfric, gauss_peak_2sm, "y")
#plt.fill_between(x_none25, gauss_peak_2.min(), gauss_peak_2, facecolor="yellow", alpha=0.5)  
plt.plot(energy_25deg_smallfric, gauss_peak_totalsm, "m")
plt.fill_between(energy_25deg_smallfric, gauss_peak_totalsm.min(), gauss_peak_totalsm, facecolor="magenta", alpha=0.5)  
plt.plot(energy_25deg_smallfric, intensity_25deg_smallfric, "k+")
plt.xlim(85,105)
plt.title('Fig. 3 - friction a=.2 ')
plt.show()


#now lets do just right friction
popt_jr, pcov_jr = curve_fit(gauss2, energy_25deg_jr, intensity_25deg_jr, p0=[amp1jr, x1jr, sigma1jr, amp2jr, x2jr, sigma2jr])
perr_jr = np.sqrt(np.diag(pcov_jr))
pars_1jr = popt_jr[0:3]
pars_2jr = popt_jr[3:6]
gauss_peak_1jr = gaus(energy_25deg_jr, *pars_1jr)
gauss_peak_2jr = gaus(energy_25deg_jr, *pars_2jr)


gauss_peak_totaljr=gauss_peak_1jr+gauss_peak_2jr

#plot the fits...only problem is it is not finding the peaks correctly
plt.plot(energy_25deg_jr, gauss_peak_1jr, "g")
#plt.fill_between(x_none25, gauss_peak_1.min(), gauss_peak_1, facecolor="green", alpha=0.5)
plt.plot(energy_25deg_jr, gauss_peak_2jr, "y")
#plt.fill_between(x_none25, gauss_peak_2.min(), gauss_peak_2, facecolor="yellow", alpha=0.5)  
plt.plot(energy_25deg_jr, gauss_peak_totaljr, "m")
plt.fill_between(energy_25deg_jr, gauss_peak_totaljr.min(), gauss_peak_totalsm, facecolor="magenta", alpha=0.5)  
plt.plot(energy_25deg_jr, intensity_25deg_jr, "k+")
plt.xlim(85,105)
plt.title('Fig. 3 -just right friction a=.02 ')
plt.show()

#next for the large friction
popt_lg, pcov_lg = curve_fit(gaus, energy_25deg_largefric, intensity_25deg_largefric, p0=[amp1lg, x1lg, sigma1lg])
perr_lg = np.sqrt(np.diag(pcov_lg))
pars_1lg = popt_lg[0:3]
gauss_peak_1lg = gaus(energy_25deg_largefric, *pars_1lg)

plt.plot(energy_25deg_jr, gauss_peak_1lg, "m")
plt.fill_between(energy_25deg_jr, gauss_peak_1lg.min(), gauss_peak_1lg, facecolor="magenta", alpha=0.5)  
plt.plot(energy_25deg_largefric, intensity_25deg_largefric, "k+")
plt.xlim(85,105)
plt.title('Fig. 3 -just right friction a=2 ')
plt.show()

'''
actualpeakswidthsg=0.1*propertieslargefric['widths']
h1g = propertieslargefric['peak_heights']
peak1widthg=actualpeakswidthsg[0]
peak2widthg=actualpeakswidthsg[1]
amp1g=h1g[0]
amp2g=h1g[1]
sigma1g=peak1widthg/2.35;
sigma2g=peak2widthg/2.35;
x1g=energy_25deg_largefric[peakslargefric[0]]
x2g=energy_25deg_largefric[peakslargefric[1]]
'''


"""
#from Pat's code:
# This function is a linear + any number of gaussians
def multiples(x, *params):
    y = x * params[0] + params[1]
    for i in range(2, len(params), 3):
        a = params[i]
        sigma = params[i+1]
        mu = params[i+2]
        y = y + gaussian(x, a, sigma, mu)
    return y
"""



"""
from running the code pars contains (x0, a, and sigma)
    pars[0] = 0.955118
    pars[1]=90.88 - looks like peak location -x0
    pars[2]=1.501 - sigma =>FWHM /approx 2.35*sigma
    https://en.wikipedia.org/wiki/Gaussian_function
"""

"""
For now i will plot the peak position as a function of friction
Then I will plot the full width half max as afunction of friction

"""
#friction values for the first peak
a_allpeak1=(0, 0.02, 0.2, 2);
#friction values for the second peak
a_allpeak2=(0, 0.02, 0.2);
FWHM_peak1=(actualpeakswidths[0], actualpeakswidthsjr[0], actualpeakswidthsm[0], actualpeakswidthslg[0])
FWHM_peak2=(actualpeakswidths[1], actualpeakswidthsjr[1], actualpeakswidthsm[1])

x1_peak1=(x1, x1jr, x1s, x1lg);
x1_peak2=(x2, x2jr, x2s)


plt.plot(a_allpeak1,FWHM_peak1,'b*', label='main peak')
plt.plot(a_allpeak2,FWHM_peak2,'m+', label='small peak')
plt.title('Fig. 3 - FWHM vs friction value ')
plt.legend()
plt.xlabel('friction value (a) [Sea-SAFARI]')
plt.ylabel('FWHM locations')
plt.show()

plt.plot(a_allpeak1,x1_peak1,'b*', label='main peak')
plt.plot(a_allpeak2,x1_peak2,'m+', label='small peak')
plt.title('Fig. 4 - peak locations vs friction value ')
plt.legend()
plt.xlabel('friction value (a) [Sea-SAFARI]')
plt.ylabel('peak locations (eV)')
plt.show()



#---------------------------------------------------------------------------
#----------------plots-------------------------------------------
#---------------------------------------------------------------------------
"""
x_fricvals25=[0, 0.02, 0.2, 2]
y_peaklocal25=[pars_nofric25[1], parssmall25[1], parlarge25[1], parjr25[1]]
#for details see
#https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html#:~:text=%2819%29%2C%20the%20standard%20deviation%20corresponds%20to%20the%20half,instead.%20This%20is%20somewhat%20larger%20than%20and%20can
FWHM25=[2.53*pars_nofric25[2], 2.53*parssmall25[2], 2.53*parlarge25[2], 2.53*parjr25[2]]


plt.plot(x_fricvals25,y_peaklocal25,'b*', label='25 deg')
plt.title('Fig. 1 - peaks vs friction value ')
plt.legend()
plt.xlabel('friction value (a) [Sea-SAFARI]')
plt.ylabel('Peak locations')
plt.show()

plt.title('Fig. 2 - FWHM vs friction value ')
plt.legend()
plt.xlabel('friction value (a) [Sea-SAFARI]')
plt.ylabel('FWHM')
plt.show()
"""