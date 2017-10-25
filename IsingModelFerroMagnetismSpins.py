#!/usr/bin/env python

###################################################

# Ryan Flamm
# This code uses the Ising model which implements the Monte Carlo method to study the changes of spins within a material. As the temperature is increased the probability of ferro-magnetism decreases. In contrast, with low temperature ferro-magnetism will occur as the system settles into having all positive or all negative spins. Two plots are produced: (1) a 2D model showing the spins over time, (2) a plot of total energy over time

###################################################
from random import choice,random,randint,uniform
from pylab import imshow,show,draw,pause,clf,plot,gray
from numpy import exp
import matplotlib.pyplot as plt
###################################################

# Finds the required energy for a given spin to flip
def Eflip(spins,x,y):
	J = -1.
        E = -2 * J * (spins[(x-1)%nx][y] + spins[(x+1)%nx][y] + spins[x][(y-1)%ny] + spins[x][(y+1)%ny]) * spins[x][y]
	return E


# Finds the total energy of the system
def totalE(spins,x,y):
        total = 0
        for i in range(nx):
                for j in range(ny):
                        total += spins[i][j]
        return total


# Variables
nx = 10	 	        # Grid size in x direction
ny = 10		        # Grid size in y direction
t = 0		        # Initial time	
tmax = 5000000		# Max time the program will reach
dt = 2		        # Timestep
kB = 1                  # Boltzmann's Constant
temp = 1                # Temperature


# Set up random initial spins in grid
spins = [[choice([-1,1]) for n in range(nx)] for m in range (ny)]


# Loop that decides whether a given spin will flip or not, according to the energy levels in the surrounding area
while t < tmax:
        for i in range(tmax):           	# For each time step
	        x = randint(0,nx-1)     	# Chooses a random x position	
	        y = randint(0,nx-1)             # Chooses a random y position
	        Energy = Eflip(spins,x,y)
                totE = totalE(spins,x,y)
                t += dt
	        if Energy < 0:
		        spins[x][y] *= -1
		        totE += Energy
	        else:	
		        P = exp(-Energy/(kB*temp))  # Boltzmann factor
		        r = uniform(0,1)	    # Random number between 0 and 1
		        if r <= P:		    # If r is less than or equal to P then the spin will flip
			        spins[x][y] *= -1
			        totE += Energy
                                
	        # Show the spin plot and energy vs. time plot
                if i % 1000 == 0:	
			
			# Spin plot
                        plt.subplot(1,2,2)
                        imshow(spins)
                        gray()
                        draw()
		        pause(10e-6)

                        # Energy vs. time plot that updates in time
                        plt.subplot(1,2,1)
                        plt.xlabel('Time')
                        plt.ylabel('Total Energy')
                        plt.scatter(t,totE)
                        plt.axis([0,50000,-150,150])
                        plt.plot(t,totE,marker='o',linestyle='--')
