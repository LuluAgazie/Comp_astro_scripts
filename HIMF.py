#!/usr/bin/env python
#This code was written by Gabriella Agazie with help from Michael Mingyar on Tuesday October 1st, 2019. 

#This code is a Monte Carlo simulation in which we create a random population of galaxies with randomly assigned masses, and use the rejection method to determine detectable galaxies and calculate the HIMF for each galaxy
#To run this simulation you will need python 2.7. There are no additional files needed to run the simulation
from random import random
import numpy as np
from matplotlib import pyplot as plt

#Our designated local number density of galaxies to be used to calacuate the HIMF of our detectable galaxies
lgd = 0.0048 
#Our maximum distance in Mpc a simulated galaxy can be from Earth. We will not consider galaxies farther away than 200Mpc, although we are aware that it is possible for a galaxy to be further than 200 Mpc. 
r_max = 200
#The area of the sky our generated galaxy population covers
surface_area = np.pi
#The total volume over which we are generating galaxies
v = 4/3*np.pi*r_max**3

#The maximum number of galaxies we can generate with our designated global number density of galaxies in our total volume
n_galaxies = int(lgd*v)

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
plt.hist(ra,bins=500)
plt.figure(2)
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
#print len(M_HI)
#print len(dist)
M_HI = np.array(M_HI)
#HIMF is the HI mass function of the detectable galaxies
HIMF = np.log(10)*lgd*(M_HI/M_char)**(a+1)*np.exp(-M_HI/M_char)
#print HIMF



#Since M_HI has extremely high values and HIMF has extremely low vales, it is helpful to plot them in logarithmic instead of linear space to better illustrate the relationship between the two 
M_HI_log = np.log10(M_HI)
HIMF_log = np.log10(HIMF)
#Plotting the distance of the detectable galaxies versus their HI mass, with the HI mass in log base 10 space
plt.figure(3)
plt.title('Distance vs Mass of Detected Galaxies')
plt.xlabel('Distance(Mpc)')
plt.ylabel('Log10(Galactic HI Mass)')
plt.scatter(dist_detect,M_HI_log)
#Plotting the HI mass of dectectable galaxies vs the corespeonding HIMF
plt.figure(4)
plt.title('Galactic HI Mass vs HI Mass Function')
plt.xlabel('Log10(HI Mass)')
plt.ylabel('Log10(HIMF)')
plt.scatter(M_HI_log,HIMF_log)
#Plotting 
plt.figure(5)
plt.hist(M_HI_log,bins = 1000)
plt.title('Distribution of Detected Galaxy HI Mass')
plt.xlabel('Log10(HI Mass)')
plt.ylabel('Number of Galaxies')
plt.show()






