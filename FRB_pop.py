#!/usr/bin/env python
#Written by Gabriella Agazie on Thursday September 19, 2019

#This code is meant to calculate linear regression of two fast radio burst (FRB) properties (signal to noise ratio and peak flux density) plotted against each other
#To run this code a file named frbcat.txt is needed with the first column containing all the SNR values and the second column containing the corresponding flux values.
#If the file has a header it must be commented out. 
#You will need python 2 to run this code

from numpy import loadtxt,savetxt,column_stack
import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erfc
from scipy.stats import rankdata

#s is the array containing all the SNR ratios of verified FRB events
s = loadtxt('frbcat.txt',float,skiprows=1,usecols=[0])
#f contains all the corresponding flux values in units of Janskys
f = loadtxt('frbcat.txt',float,skiprows=1,usecols=[1])
#N is the number of data points we have for SNR and Flux
N = len(s)

#Calculations of mean, meadian, variance, standard deviation, and variance for the peak flux density
print "Average flux: %s Jy" % np.mean(f)
print "Median flux: %s Jy" % np.median(f)
print "Standard Deviation of flux: %s" % np.std(f)
print "Variance of flux: %s" % np.var(f)

#Plotting SNR vs flux for the FRBs in the catalog
plt.scatter(s,f)
plt.xlabel("Signal to Noise Ratio")
plt.ylabel("Peak Flux Density (Jy)")
plt.title("Signal to Noise Ration vs Peak Flux Density of Verified FRBs")


#The following variables are used to help simplify our linear regression calculations 
#s_2 is the array of squared SNR values
s_2 = s**2
#s_2_mean is the average value of s_2
s_2_mean = np.mean(s_2)
#s_mean is the average value of s
s_mean = np.mean(s)
#sf_mean is the average value of the product of the s and f
sf_mean = np.mean(s*f) 
#f_mean is the average value of f
f_mean = np.mean(f)



###Calculations for our linear regression fit line###
#a is the slope of our best fit line
a = (sf_mean-(s_mean*f_mean))/(s_2_mean-(s_mean**2))
#b is the y intercept of our beest fit line
b = ((s_2_mean*f_mean)-(s_mean*sf_mean))/(s_2_mean-(s_mean**2))

#f_fit is our equation for the best fit line 
f_fit = a*s+b

#S_yx is the scatter of f_fit calculated from s
S_yx = np.sqrt(np.sum((f-f_fit)**2)/(N-2))
#S_y is the scatter of the flux values from the catalog
S_y = np.sqrt(np.sum((f-f_mean)**2)/(N-1))
#S_x is the scatter of the SNR values
S_x = np.sqrt(np.sum((s-s_mean)**2)/(N-1))

#S_a is the scatter of the slope of the best fit line, a
S_a = S_yx/(S_x*np.sqrt(N-1))
#S_b is the scatter of the y-intercept of the best fit line, b
S_b = S_a*f_mean
print "S_a is:%s" % S_a
print "S_b is:%s" % S_b

##Calculation for Pearson correlation coefficient(r) and population correlation coefficient (rho)
r = a * (S_x/S_y)
print "r is %s" % r
rho = erfc((r*np.sqrt(N))/np.sqrt(2))
print "rho is:%s" % rho

#t_r is the signicance of r which should be less than 0.05 to reject the null hypothesis
t_r = r*np.sqrt((N-2)/(1-r**2))
print "Significance of r is: %s" % t_r


#R_i is an array in which every SNR value has been assigned an interger value rank. Any identical SNR values are assigned half interger ranks
R_i = rankdata(s)
#print R_i
#R is the average value of R_i
R = np.mean(R_i)
#S_i is an array containing the rank of each flux value 
S_i = rankdata(f)
#print S_i
#S is the average value of S_i
S = np.mean(S_i)

##Calculation of Spearman rank coefficient r_s if R_i and S_i contain no half integers
#D = np.sum((R_i-S_i)**2)
#print D
#r_s = 1-6*D/(N**3-N)

##Calculation for Spearman rank coefficent r_s since both R_i and S_i contain half integers
r_s = np.sum(((R_i-R)*(S_i-S))/(np.sqrt(np.sum((R_i-R)**2))*np.sqrt(np.sum((S_i-S)**2))))
print "r_s is:%s" % r_s
#t is the significance of r_s which should be less that 0.05 to reject the null hypothesis.
t = r_s*np.sqrt((N-2)/(1-r_s**2))
print "Significance of r_s is: %s" % t
##Calculation of the Fischer Transform using r_s
F_r = np.arctanh(r_s)
print "Fisher transform is:%s" % F_r

##Calculation of z value for given Fischer transform
z = np.sqrt((N-3)/1.03)*F_r
print z
plt.plot(s,f_fit)
plt.show()



