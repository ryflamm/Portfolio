#!/usr/bin/env python

##################################################################
#Ryan Flamm
#10.8
#this code explores the variational monte carlo method
##################################################################
 
import matplotlib.pyplot as plt
from random import randint,random,sample,uniform,choice
from numpy import array,linalg,where,sum,exp,sin,abs,pi,arctan, arange,sqrt,cos,linspace,sign,copy
from numpy.linalg import norm

##################################################################


#Lennard-Jones Potential
def V(x):
    sigma = 1
    epsilon = 1
    r6 = (x/sigma)**(-6)
    potential = 4 * epsilon * (r6**2 - r6)
    return potential

#calculates energy 
def calculateEnergy(psi,pot,dx,x):
    E0 = 0
    summ = 0
    for i in range(0,N-1):
         if 0 < psi[x] < 47:
             energy = E0 + dx * V(x)*psi[x]**2 - dx * psi[x]*(psi[x+1] + psi[x-1] - 2*psi[x])/(dx**2)
             summ = summ + psi[x] * psi[x] * dx
             energy = energy / summ
         else:
             energy =  0
    return energy

#updates energy with random change in psi 
def newCalculateEnergy(psi2,pot,dx,x):
    E0 = 0
    summ = 0
    for i in range(0,N-1):
         if 0 < psi2[x] < 47:
             energy = E0 + dx * V(x)*psi2[x]**2 - dx * psi2[x]*(psi2[x+1] + psi2[x-1] - 2*psi2[x])/(dx**2)
             summ = summ + psi2[x] * psi2[x] * dx
             energy = energy / summ
         else:
             energy = 0
    return energy
     
#bounds of problem  
lBound = 0.7
rBound = 5
N = 50
dx = (rBound - lBound)/N
xDomain = arange(lBound,rBound+dx,dx)

#initial conditions
psi = array([0.5 for x in xDomain])  # Initialiize psi, changing value here will change test function
psi[:3] = 0. #sets bounds on the function
psi[46:] = 0.#such that for x<3 and x>46 the function is equal to zero
psi = psi/sqrt(sum(psi * psi * dx))  # normalize psi
pot = array([V(n) for  n in xDomain]) # get discrete version of potential
energy = calculateEnergy(psi,pot,dx,x)  # calculate initial energy
 
 
index = 0
while index < 10000000:  # Main loop
    for i in range (N):
        x = randint(1,N-1) #choose random location in trial function
        deltaPsi = randint(-10,20) #choose random change for psi at x
        psi2 = psi + deltaPsi #Adjust psi
        energy = calculateEnergy(psi,pot,dx,x) #calculate original energy
        newEnergy = newCalculateEnergy(psi2,pot,dx,x) #calculate adjusted energy
        if newEnergy < energy: #replace old energy
            energy = newEnergy
        if newEnergy < energy:
            psi = psi2/sqrt(sum(psi * psi * dx)) #replace old psi
        if index % 10000 == 0: #plot occasionally
            plt.scatter(xDomain,psi)
            plt.xlim(0,rBound+dx)
            #plt.ioff()
            plt.draw()
            plt.pause(1e-20)
            #plt.savefig('foo.png')
            #plt.show
            plt.clf()
            index += 1

#I was unable to get the code to properly plot. I had a hard time figuring out what the issue was because running the code does not call any errors. The plot box pops up, but it stays gray the whole time. However, by doing savefig() I was able to at least get a plot of the initial conditions (trial function). I also printed all of the energies and psi values to see if there was any issue, but was unsucessful in finding any cause of the gray plot.            
