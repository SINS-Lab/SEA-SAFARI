# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:23:20 2020

@author: dave

Take octave code from deg25.m and put this into python syntax

"""


from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
#import argparse                     # Parsing command line arguments
import numpy as np                  # Array manipulation/maths
import matplotlib.pyplot as plt     # Plotting
#import os                           # Path related stuff
import math

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

"""
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
localmax=[]
localindex=[]
localenergy=[]
localval=[]
localindexup=[]
localvalsdown=[]
localindexdown=[]

for jj  in range(len(energy_25deg_nofric)-2):
    #now lets get the local max values *first if statement helps not look at noise levels
    if intensity_25deg_nofric[jj] >= 0.1:
        if intensity_25deg_nofric[jj]>intensity_25deg_nofric[jj-1] and intensity_25deg_nofric[jj]>intensity_25deg_nofric[jj+1]:
            localmax.append(intensity_25deg_nofric[jj])
            localenergy.append(energy_25deg_nofric[jj])
            localindex.append([jj])
        #end if statement
    #end outer if statement
#end for loop

FWHM=[]
#some funkiness with this compared to octave
for kk in range(len(energy_25deg_nofric)-2):
    for mm in range(len(localmax)):
        #if intensity_25deg_nofric[kk] >= 0.1:
            if (intensity_25deg_nofric[kk] <= (localmax[mm]/2) and  (intensity_25deg_nofric[kk+1] >= (localmax[mm]/2))) or ((intensity_25deg_nofric[kk-1] >= (localmax[mm]/2)) and  (intensity_25deg_nofric[kk] <= (localmax[mm]/2))):
                #we are approximately at the fullwidth half max location of peaks (the actual fwhm is every two data values in here)
                FWHM.append(energy_25deg_nofric[kk])
            #end if statement
       #end if statement
    #end for loop (mm)

    #i establish anything below 0.1 as the noise
    if intensity_25deg_nofric[kk] >= 0.1:
        if  (intensity_25deg_nofric[kk-1]  < intensity_25deg_nofric[kk] ) and  (intensity_25deg_nofric[kk] < intensity_25deg_nofric[kk+1]):
            localval.append(intensity_25deg_nofric[kk])
            localindexup.append(energy_25deg_nofric[kk])


        elif (intensity_25deg_nofric[kk] < intensity_25deg_nofric[kk-1]) and (intensity_25deg_nofric[kk] > intensity_25deg_nofric[kk+1]): 

           localvalsdown.append(intensity_25deg_nofric[kk])
           localindexdown.append(energy_25deg_nofric[kk])
   
       #endif
  #endif
#end for loop kk
            



plt.plot(localindexup, localval , '.k', label='vals up')
plt.plot(localindexdown, localvalsdown , '.m', label='vals down')
plt.plot(localenergy, localmax , 'xb', label='max')
plt.title('Fig. 1 local values ')
plt.legend()
plt.xlabel('energy ')
plt.ylabel('intensity')
plt.xlim(86, 95)
plt.show()



#o get an approximate of the full width
energysave=[]

for nn in range(len(localvalsdown)-1): 
  #attempt to get the first peaks right most data point width
   if localvalsdown[nn] < localvalsdown[nn+1]:
     #then we should be at the next set of values going update
     energysave.append(localindexdown[nn])
   #endif
#endfor

energysave1=[]
for oo in range(len(localval)-1):
  #attempt to get the second peak's leftmost datapoint
   if localval[oo] > localval[oo+1]:
     #then we should be at the next set of values going update
     energysave1.append(localindexup[oo])
  #endif
#endfor

peak1width= (FWHM[2]-FWHM[1])
peak2width= (FWHM[3]-FWHM[2])#(localindexdown[len(localindexdown)-1]-energysave[0])

sigma1=peak1width/2.35;
sigma2=peak2width/2.35;
amp1=localmax[0];
amp2=localmax[1];
x1=localenergy[0]
x2=localenergy[1];

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
plt.show()


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

    
#This is the guess of the peak location. 
x00=1

pars_nofric25,cov_none = curve_fit(gaus,x_none25,y_none25,p0=[x00,mean25,sigma25])

#pars_nofric250,cov_none0 = curve_fit(gaus,x_none25,y_none25,p0=[x00,mean250,sigma250])


plt.plot(energy_25deg_nofric,intensity_25deg_nofric,'b+:',label='data')
#plt.plot(x_none25,gaus(x_none25,*pars_nofric25),'ro:',label='fit')
plt.plot(x_none25,gaus(x_none25,*pars_nofric25),'k*:',label='weird fit')
plt.legend()
plt.title('Fig. 3 -rando code')
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity')
plt.xlim(85,105)
plt.show()


#now we do the same fitting with the small friction a *SEASAFARI* = 0.2
x_small25=np.array(energy_25deg_smallfric)
y_small25=np.array( intensity_25deg_smallfric)
nsmall25 = len(x_small25)                          #length/number of data elements
means25 = sum(x_small25*y_small25)/nsmall25                  
sigmas25 = sum(y_small25*(x_small25)**2)/nsmall25       
parssmall25,cov_sm25 = curve_fit(gaus,x_small25,y_small25,p0=[1,means25,sigmas25])

#now we do the same fitting with the small friction a *SEASAFARI* = 2
x_large25=np.array( energy_25deg_largefric)
y_large25=np.array(   intensity_25deg_largefric)
nlarge25 = len(x_large25)                          
meanlg25 = sum(x_large25*y_large25)/nlarge25                   
sigmalg25 = sum(y_large25*(x_large25)**2)/nlarge25        
parlarge25,cov_lg25 = curve_fit(gaus,x_large25,y_large25,p0=[1,meanlg25,sigmalg25])

#now we do the same fitting with the small friction a *SEASAFARI* = .02
x_jr25=np.array(energy_25deg_jr)
y_jr25=np.array(intensity_25deg_jr)
njr25 = len(x_jr25)                          
meanjr25 = sum(x_jr25*y_jr25)/njr25                   
sigmajr25 = sum(y_jr25*(x_jr25)**2)/njr25        
parjr25,cov_jr25 = curve_fit(gaus,x_jr25,y_jr25,p0=[1,meanjr25,sigmajr25])


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