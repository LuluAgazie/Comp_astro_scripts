#!/usr/bin/env python
#This code was written by Gabriella Agazie with help from Michael Mingyar on Tuesday October 1st, 2019.

#This code is a Monte Carlo simulation in which we create a random population of galaxies with randomly assigned masses,
#To run this simulation you will need python 2.7. There are no additional files needed to run the simulation
from random import random
import numpy as np
import math
from matplotlib import pyplot as plt
import scipy.integrate as integrate

#Our designated local number density of galaxies to be used to calacuate the HIMF of our detectable galaxies
lgd = 0.0048*5  #Over 5 dex
phi_star = 0.0048
#Our maximum distance in Mpc a simulated galaxy can be from Earth. We will not consider galaxies farther away than 200Mpc, although we are aware that it is possible for a galaxy to be further than 200 Mpc.
r_max = 200
#The area of the sky our generated galaxy population covers
surface_area = np.pi
#The total volume over which we are generating galaxies
v = 4/3*np.pi*(r_max**3)

#The maximum number of galaxies we can generate with our designated global number density of galaxies in our total volume
n_galaxies = int(lgd*v)
print n_galaxies
#Right Ascension of generated galaxies will be stored in ra
ra = []
#Declination of generated galaxies will be stored in dec
dec = []
#Distance from Earth to generated galaxies will be stored in dist
dist = []
#Inclination angle of generated galaxies will be stored in incl
incl = []



#i is the running total of generated galaxies within our desired volume. The while loop will run until this number is equivalent to our desired number of generated galaxies.
i = 0
#This loop first generates the x,y,and z coordinates of a galaxy, and then calculates the radial distance from Earth in Mpc
#The if loop determines if the radial distance is less than 200 Mpc and if so, it will calculate and record the right ascension, declination, and distance of the galaxy.
#It will also increase i by one.
#If the radial distance of the generated galaxy is greater than 200 Mpc, it will not be stored and the while loop will restart without increasing i.
while i <= n_galaxies:
    x=random()*2*r_max-r_max
    y=random()*2*r_max-r_max
    z=random()*2*r_max-r_max
    r = np.sqrt(x**2+y**2+z**2)
    if r <= r_max:
        i+=1
        theta = np.arctan2(x,y)
        ra.append(theta*180/np.pi)
        phi = np.arctan2(z,np.sqrt(x**2+y**2))
        dec.append(phi*180/np.pi)
        incl.append(phi*180/np.pi)
        dist.append(r)

#print np.mean(ra)
#print np.mean(dec)
#print np.mean(dist)

#Check if distribution of ra is uniform and dec is gaussian
plt.figure(1)
plt.title("Distribution of Right Ascensions for Simulated Galaxies")
plt.xlabel("Right Ascension (deg)")
plt.ylabel("Number of Galaxies")
plt.hist(ra,bins=500)
plt.figure(2)
plt.title("Distribution of Declinations for Simulated Galaxies")
plt.xlabel("Declination (deg)")
plt.ylabel("Number of Galaxies")
plt.hist(dec,bins=500)


plt.show()



#HI Mass of the generated galaxies that have a flux density above the detection threshold of 0.72 Jy km/s.
M_HI = []
#Flux density of the generated galaxies that are above the detection threshold
S_HI = []
#The distance from Earth of the detectable galaxies from the poplution generated above.
dist_detect = []

#a is the faint end slope of the HIMF
a = -1.33
#M_char is a characteristic mass in terms of the mass of the sun and the lgd
M_char = 10**(9.96)

#j is the index of the particular generated galaxy we are assigning an HI mass to.
j = 0
#This while loop will assign an HI mass (m) to each galaxy we generated in the first while loop and then calculate what the flux density of the galaxy would be from its distance and HI mass.
#If the flux density is above our detection threshold, we will store it as a detecable galaxy.
#There is no need to store the HI masses of nondetectable galaxies
while j <= n_galaxies:
    n = random()*5+6
    m = 10**n
    s = m/((2.36e5)*(dist[j])**2)
    if s >= 0.72:
        M_HI.append(m)
        S_HI.append(s)
        dist_detect.append(dist[j])
    j+=1
#print M_HI
print len(M_HI)
print len(dist)
M_HI = np.array(M_HI)
#HIMF is the HI mass function of the detectable galaxies
HIMF = np.log(10)*phi_star*(M_HI/M_char)**(a+1)*np.exp(-M_HI/M_char)
#print HIMF



#Since M_HI has extremely high values and HIMF has extremely low vales, it is helpful to plot them in logarithmic instead of linear space to better illustrate the relationship between the two
M_HI_log = np.log10(M_HI)
HIMF_log = np.log10(HIMF)
#Plotting the distance of the detectable galaxies versus their HI mass, with the HI mass in log base 10 space
plt.figure(3)
plt.title('Distance vs Mass of Detected Galaxies')
plt.xlabel('Distance(Mpc)')
plt.ylabel('Log10(Galactic HI Mass)')
plt.scatter(dist_detect,M_HI_log,s=1)
#Plotting the HI mass of dectectable galaxies vs the corespeonding HIMF
plt.figure(4)
plt.title('Galactic HI Mass vs HI Mass Function')
plt.xlabel('Log10(HI Mass)')
plt.ylabel('Log10(HIMF)')
plt.scatter(M_HI_log,HIMF_log, s=2)
#Plotting
plt.figure(5)
plt.hist(M_HI_log,bins = 1000)
plt.title('Distribution of Detected Galaxy HI Mass')
plt.xlabel('Log10(HI Mass)')
plt.ylabel('Number of Galaxies')
plt.show()



###This section of the code was written on October 24, 2019 by Gabriella Agazie with help from Michael Mingyar###
#This section is used to conduct a Bayesian analysis of galaxies within a 7 Mpc distance from Earth in order to determine a value of the faint-end slope of the HIMF


#The distance from Earth of all detected galaxies within 7 Mpc of Earth
dist_detect_7 = []
#The HI masses of detctable galaxies within 7 Mpc of Earth
M_HI_7 = []
#The number of detected galaxies
N_detect = len(dist_detect)

#Iterator to be used in the following while loop
k = 0
#This loop looks at each distance of a detected galaxy, dist_detect[k], and determines if it is within 7 Mpc of Earth. If so, the distance and the corresponding HI mass, M_HI[k], are saved in dist_detect_20 and M_HI_20 respectively.
while k <= N_detect-1:
    if dist_detect[k] <= 7:
        dist_detect_7.append(dist_detect[k])
        M_HI_7.append(M_HI[k])
    k+=1

#The number of detected galaxies within 7 Mpc. It will serve as the prior for our Bayesian analysis
m = len(M_HI_7)

#Converting the list of masses to an array for later calculations 
M_HI_7 = np.array(M_HI_7)



#The lowest possible HI mass produced by our simulation
M_min = 10**6

#The minimum and maximum possible values for the HIMF power law function normalization constant, K. These will serve as the bounds of our integral when solving for the posterior
K_min = 1e-6
K_max = 27

#The array where we will store our posteriors calculated for each possible value of the faint-end slope, alpha. 
posteriors = []

#This function takes a value of alpha and integrates the posterior probabality density function for alpha over K. It returns the posterior value for each value of alpha it is given  
def posterior(alpha):
    p_prob = lambda K: ((K*(M_min**(alpha+2)))**m*np.exp(-(K*(M_min**(alpha+2)))))/math.factorial(m)
    post = integrate.quad(p_prob,K_min,K_max)
    post = post[0]
    return post

#This is an array of possible values of alpha that we will calculate posteriors for.
alpha_inputs = np.linspace(-1.0,-1.9,100)

#This loop takes each possible value of alpha and runs the posterior function on it, and then saves it to the array posteriors
for i in alpha_inputs:
    posteriors.append(posterior(i))

#Plotting a bar graph of the posterior values for each value of alpha
plt.bar(alpha_inputs,posteriors, width = 0.05)
plt.title("Probability Density Function of Alpha")
plt.xlabel("Alpha")
plt.ylabel("Posterior")
plt.xlim(-1.9,-1.0)
plt.show()


