#!/usr/bin/env python
#This code was written by Gabriella Agazie on November 23, 2019

#The purpose of this code is to calculate the number density of known FRBs at several distance intervals, using data from the FRB catalog (http://frbcat.org/), and compare this to the number density distrbution predicted by the star formation rate equation.  
#We used the 2 sample KS test to compare the two distributions
#To run this code you will need the file frbcat_data.txt, where the first column is the FRB name, the second is the telescope that first detected the FRB, the third is right ascension, the fourth is redshift, the fifth is luminosity distance and the sixth is comoving distance.
#If the file has a header it must be commented out
#You will need python 2.7 to run this code
from numpy import loadtxt, savetxt, column_stack
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

#The name of the FRB
FRB_name = loadtxt('frbcat_data.txt', str, usecols=[0])
#The estimated redshift of the FRB
FRB_z = loadtxt('frbcat_data.txt', float, usecols=[4])
#The maximum luminosity distance in Gpc. 
D_lum = loadtxt('frbcat_data.txt', float, usecols=[5])
#The comoving distance estimated from the redshift in Gpc.
D_comov = loadtxt('frbcat_data.txt', float, usecols=[6])

#Putting the FRB comoving distances into 7 bins, so that the number of FRBs per bin is greater than or equal to one and plotting the distribution 
bining = plt.hist(D_comov, bins = 7)
plt.title("Distribution of Comoving Distance from Earth of FRBs")
plt.xlabel("Distance (Gpc)")
plt.ylabel("Number of FRBs")
plt.show()
#Putting the FRB red shifts into 7 bins, so that the number of FRBs per bin is greater than or equal to one and plotting the distribution
bins = plt.hist(FRB_z,bins = 7)
plt.title("Distribution of Estimated Redshifts of FRBs")
plt.xlabel("Redshift")
plt.ylabel("Number of FRBs")
plt.show()
#Number of FRBs per redshift bins
z_per_bin = bins[0]
#Bounds of each redshift bin
z_bounds = bins[1]
#Number of FRBs per comoving distance bin
N = bining[0]
#Bounds of each distance bin in Gpc
D = bining[1]
#The average comoving distance for each bin
D_mean = []


#Average redshift per bin 
z = []

#This loop calculates the average comoving distance per bin and stores it in an array for further calculations
for i in range(7):
    D_mean.append(np.mean([D[i],D[i+1]]))
D_mean = np.array(D_mean)

#This loop calculates the average redshift per bin and stores it in an array for further calculations
for h in range(7):
    z.append(np.mean([z_bounds[h],z_bounds[h+1]]))
z = np.array(z)
print z
#The volume per bin in Gpc^3
v = []
#This loop calculates the volume per bin in units of Gpc^3 and stores it in an array for further calculations.
for k in range(7):
    v.append(((4./3.)*np.pi*D[k+1]**3)-((4./3.)*np.pi*D[k]**3))
v = np.array(v)
#The calculated number density of known FRBS
number_density = N/v

#The projected number density of FRBs based off of the star formation rate equation from Gardinier et al., 2019
p_frb_z = ((1+z)**(2.7))/(1+((1+z)/2.9)**(5.6))
#The comoving distance associated with each redshift calculated using the cosmological calculator (http://www.astro.ucla.edu/~wright/CosmoCalc.html)
D_z = np.array([0.8371,1.9442,2.873,3.6506,4.3051,4.861,5.3411])


#Plotting both the projected number density of FRBs and the measured number density vs the comoving distance
plt.title("FRB Number Density vs Distance")
plt.ylabel("Number Density")
plt.xlabel("Distance (Gpc)")
plt.scatter(D_z,p_frb_z,label="Projected")
plt.scatter(D_mean, number_density,label="Measured")
plt.legend()
plt.show()


#Running a KS test to determine whether the difference between the measured and projected number densities is statisically significant
KS_test = stats.ks_2samp(number_density,(p_frb_z))
print KS_test






