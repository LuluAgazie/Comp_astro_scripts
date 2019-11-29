#!/usr/bin/env python

#This code was written by Gabriella Agazie with help from Michael Mingyar and Haley Wahl on November 7, 2019
#This code is meant to analyze 21-cm HI data from a galaxy and calculate the HI mass and dynamical mass of the galaxy.
#To run the code you will need a file named UGC328.dat with the first column containing the radio frequency in GHz, the second containing values of intensity on the source, and the third containing intensity values of an off source region of the sky.
#This program requires python 2.7 to run. 

from numpy import loadtxt, savetxt, column_stack
import numpy as np
from matplotlib import pyplot as plt

#The system temperature in units of Kelvin of the telescope used to take the observation. In this case, the Robert C. Byrd Green Bank Telescope
Tsys = 20
#The gain of the telescope receiver in units of K/Jy
gain = 2
#The radio frequencies used to observe the source
freq = loadtxt('UGC328.dat', float, usecols=[0])
#The intensity at each radio frequency of the source in arbitrary units
on_source = loadtxt('UGC328.dat', float, usecols=[1])
#The intensity at each radio frequency of an off source region of the sky in arbitrary units
off_source = loadtxt('UGC328.dat', float, usecols=[2])
#The speed of light in units of m/s
c = 3e8

#Plotting the on source intensity against radio frequency
plt.figure(1)
plt.scatter(freq,on_source, s=2)
plt.title("On source HI intensity")
plt.xlabel("Radio Frequency (GHz)")
plt.xlim(min(freq),max(freq))
plt.ylabel("Intensity")
#Plotting the off source intensity against radio frequency
plt.figure(2)
plt.scatter(freq,off_source, s=2)
plt.title("Off source HI intensity")
plt.xlabel("Radio Frequency (GHz)")
plt.xlim(min(freq),max(freq))
plt.ylabel("Intensity")

plt.show()


#The on source intensity with the off source background subtracted off and then normalized by the off source intensity. Still in arbitrary units
diff = (on_source - off_source)/ off_source
#The normalized flux density of the source converted to units of Jy.
diff = diff*Tsys/gain
#Plotting the flux density of the source against radio frequency
plt.figure(3)
plt.scatter(freq, diff, s=2)
plt.title("Normalized HI Intensity")
plt.xlabel("RadioFrequency (GHz)")
plt.ylabel("Flux Density (Jy)")
plt.xlim(min(freq),max(freq))
plt.show()


#The rest frequency of HI
f_rest = 1.4204057177
#Using the doppler formula to transform frequency scale to velocity scale
v = (c*((f_rest/freq)-1))/1000

#Plotting the flux density of the source against velocity with the minimum and maximum velocities of the emission line marked with red lines and the mean velocity marked with a green line. 
plt.plot(v,diff)
plt.title("HI Emission Line for UGC 328")
plt.xlabel("Velocity (km/s)")
plt.ylabel("Flux Density (Jy)")
plt.xlim(1850,2150)
plt.ylim(-0.05,0.2)
plt.axvline(x=1905, c = 'red')
plt.axvline(x=2075, c = 'red')
plt.axvline(x=1990, c = 'green')
plt.show()


#The mean velocity of the emission line in km/s
v_mid = 1990
#The maximum velocity of the emission line in km/s
v_max = 2075
#The minimum velocity of the emission line in km/s
v_min = 1905
#The mass of the sun in kg
M_sun = 2e30
#The width of the emission line in km/s
W = v_max-v_min
#Hubble's constant in Km/s/Mpc
H = 65 
#The estimated distance of UGC 328 from Earth
D = v_mid/H
#Newton's gravitational constant in unit of (km/s)^2 pc/solar masses
G = 4.3e-3
print "UGC 328 is %d Mpc away" % D
#The radius of UGC 328 as given by Pisano et. al. 2002
R = (20*1000)
#The inclination angle of UGC 328 as given by Hyperleda
inc_r = 24.9*(np.pi/180)
#The rotational velocity of UGC 328 corrected for inclination
V_rot = (W/2)/np.sin(inc_r)
#The dynamical mass of UCG 328
M_dyn = ((V_rot**2)*R)/G
print "The dynamical mass of UGC 328 is %.2E solar masses" % M_dyn


#Velocity values within the emission line shown in the velocity vs flux plot 
v_line=[]
#Corresponding flux values to the velocity values within the emission line.
f_line=[]

#This while loop selects velocity values and the corresponding fluxes that are within the emission line in velocity vs flux plot above and saves them to new arrays to be used to calculate HI mass.
i = 0
while i <= (len(v)-1):
    if v_max>=v[i]>=v_min:
        v_line.append(v[i])
        f_line.append(diff[i])
    i+=1

#Equation to calculate the HI mass of UGC 328. 
M_HI = (2.36e5)*(D**2)*abs(np.trapz(f_line,v_line))
print "The HI mass of UGC 328 is %.2E solar masses" % M_HI




