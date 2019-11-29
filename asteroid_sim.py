#!/usr/bin/env python

#Written by Gabriella Agazie with help from Michael Mingyar and Haley Wahl on Wednesday September 18, 2019

#This code is meant to simulate the orbits of ten asteroids that are influenced by the gravitational force of the Sun and Jupiter. 
#The goal is to demonstrate the lack of stable asteroid orbit at the Kirkwood gap regions in the asteroid belt.
#The Sun is treated as a fixed point about which the asteroids and Jupiter orbit.  
#We are treating the gravitational force of the asteroids on each other and the asteroids on Jupiter to be negligible.
#The code is run as a python script. python asteroid_sim.py
#To run this code you will need python 2

#To run this code you will need the numpy and pyplot python packages
import numpy as np
from matplotlib import pyplot as plt


#Constants needed for calculating orbits
#GM_sun is Newton's gravitational constant multiplied by the mass of the sun in units of AU^3/yr^2.
GM_sun = 4*(np.pi)**2  
#GM_jup is Newton's gravitational constant multiplied by the mass of Jupiter, which we can calculate from GM_sun by multiplying it by the mass ratio of the Sun and Jupiter.
#We are making Jupiter's mass 10x its actual mass to better demonstrate the effects of Jupiter's gravitational force on the asteroid orbits.
GM_jup = GM_sun*0.009547
#All positon (x and y) and distance (r and semimajor axis) factors are in units of AU
#All velocity factors are in units of AU/yr
#All time factors are in units of years

###This section is calculating the orbit of jupiter around a fixed sun at (0,0).###
#The initial conditions for Jupiter.
#a is Jupiter's semimajor axis in units of AU
a = 5.2
#x_i is the inital x coordinate of Jupiter's position in an xy-plane
x_i = a
#y_i is the inital y coordinate of Jupiter's position in an xy-plane
y_i = 0 
#v_xi is the initial x-component of Jupiter's velocity vector
v_xi = 0
#v_yi is the initial y-compnent of Jupiter's velocity vector which is calculated by setting the gravitational force of the Sun acting on Jupiter to the centripedal force of Jupiter and then solving for velocity.
v_yi = np.sqrt(GM_sun/a)

#x is the array in which all calculated x-positions of Jupiter will be stored
x = [a]
#y is the array in which all calculated y-positions of Jupiter will be stored
y = [0.]
#v_x is the array in which all calculated x-velocity components of Jupiter will be stored
v_x = [0.]
#v_y is the array in which all calculated y-velocity components of Jupiter will be stored
v_y = [v_yi]

#t is the array containing each time after the start of the simulation that the orbits of Jupiter and the astroids were calculated.
t = [0.]
#dt is the time step over which the orbit is calculated. 
dt = 0.01


#This loop calculates the x and y positions of Jupiter as well as the x and y velocity components, and using the previous values to calculate the next values.
#r is the radius of Jupiter's orbit, which is recalculated for each iteration of the loop.
#i is the particular iteration number the loop is on
#The loop is allowed to run over 100000 iterations
i = 0
while i <= 100000:
    t.append(t[i]+dt)
    r = np.sqrt((x[i]**2)+(y[i]**2))
    v_x.append(v_x[i]-(GM_sun*x[i]*dt/(r**3)))
    v_y.append(v_y[i]-(GM_sun*y[i]*dt/(r**3)))
    x.append(x[i]+v_x[i+1]*dt)
    y.append(y[i]+v_y[i+1]*dt)
    i+=1

#time_span is the length of time the simulation runs over, which is the last value of t array
time_span = t[-1] 
print "Simulated time span: %s" % time_span
    
###This section is calculating the orbits of the astroids around a fixed sun at (0,0) with the orbit of Jupiter affecting the asteroid orbits###

#The function asteroids calculates and plots the orbit of an astroid when given a float for the semimajor axis of the asteroid and a string containing the color the orbit should be plotted with.


#The variables and arrays within the function asteroids:
#aa is the semi-major axis of the asteroid.
#x_a_i is the initial x-position of the asteroid orbit.
#y_a_i is the initial y-position of the asteroid orbit.
#v_y_a_i is the initial y component of the asteroid velocity vector which is calculated by setting the gravitational force of the Sun acting on the asteroid to the centripedal force of the asteroid and then solving for velocity.
#v_x_a_i is the initial x component of the asteroid velocity vector.
#x_ja_i is the initial x component of the separation vector between Jupiter and the asteroid. We calculate this by subtracting x_i and x_a_i
#y_ja_i is the initial y component of the separation vector between Jupiter and the asteroid. We set it to zero to simplify our calculation of x_ja_i
#x_a is the array in which all calculated x positions of the asteroid orbit will be stored.
#y_a is the array in which all calculated y positions of the asteroid orbit will be stored.
#v_x_a is the array in which all calculated x compenents of the astroid velocity will be stored. 
#v_y_a is the array in which all calculated x compenents of the astroid velocity will be stored.
#x_ja is the array in which all calculated x components of separation vector between Jupiter and the asteroid will be stored
#y_ja is the array in which all calculated y components of separation vector between Jupiter and the asteroid will be stored
#r_a is the radius of the asteroid orbit, which is recalculated for each iteration of the loop within the function asteroids.
#r_ja is the magnitude of the separation vector between Jupiter and the asteroid, which is recalculated for each iteration of the loop within the function asteroids.
#j is the number of iterations over which the asteroid orbit is calculated. 

#asteroids calculates the asteroid orbit in the same manner and over the same time scale as Jupiter's orbit. 
#The only difference is that the calculations for the asteroid orbit take into account the effect of Jupiter's gravitational force acting on the asteroid as well as the Sun's.
def asteroids(semi_major_axis=float,asteroid_color=str):
    aa = semi_major_axis
    x_a_i = aa
    y_a_i = 0.0
    v_y_a_i = np.sqrt(GM_sun/aa)
    v_x_a_i = 0.
    x_ja_i = a - aa
    y_ja_i = 0.

    x_a = [aa]
    y_a = [0.]
    v_x_a = [0.]
    v_y_a = [v_y_a_i]
    x_ja = []
    y_ja = []

    j = 0 
    while j <= 100000:
        x_ja.append((x[j]-x_a[j]))
        y_ja.append((y[j]-y_a[j]))
        r_a = np.sqrt((x_a[j]**2)+(y_a[j]**2))
        r_ja = np.sqrt((x_ja[j]**2)+(y_ja[j]**2))
        v_x_a.append(v_x_a[j]-(((GM_sun*x_a[j]/r_a**3)+(GM_jup*x_ja[j]/r_ja**3))*dt))
        v_y_a.append(v_y_a[j]-(((GM_sun*y_a[j]/r_a**3)+(GM_jup*y_ja[j]/r_ja**3))*dt))
        x_a.append(x_a[j]+v_x_a[j+1]*dt)
        y_a.append(y_a[j]+v_y_a[j+1]*dt)
        j+=1
#Here we are plotting all the calculated x and y postions calculated above in order to visually represent the change in the asteroids orbit over a 1000 year time scale.
    plt.plot(x_a,y_a,asteroid_color,linewidth=0.1,label="%s orbit" % str(semi_major_axis))

#Here we run the asteroids function for 11 different semi-major axes and thus are calculating and plotting the orbits of 11 asteroids
asteroids(3.3,'thistle')
asteroids(3.27,'darkcyan') 
asteroids(3.2,'crimson')
asteroids(3.1,'orchid')
asteroids(2.95,'paleturquoise')
asteroids(2.82,'coral')
asteroids(2.7,'greenyellow')
asteroids(2.5,'gold')
asteroids(2.4,'magenta')
asteroids(2.2,'indigo')
asteroids(2.06,'red')


plt.title("Orbit of Jupiter and 11 Astroids over 1000 years")
plt.xlabel("X position of orbit (AU)")
plt.ylabel("Y position of orbit (AU)")
plt.xlim(-6,6)
plt.ylim(-6,6)
plt.axis('equal')
#Below we are plotting the position of the Sun as a yellow dot.
plt.plot(0.,0.,'yo',label='Sun')
#We are plotting the orbit of Jupiter as a blue line.
plt.plot(x,y,'b-',label='Jupiter orbit')
plt.legend()
plt.show()        







