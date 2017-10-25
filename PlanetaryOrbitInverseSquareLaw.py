#!/usr/bin/env python

###################################################################
# Ryan Flamm
# This program puts to question the inverse square law as it relates to the orbits of planets. By changing the exponent in the denominator we can see how the planets would act if they did not follow the inverse square law.

###################################################################
from numpy import arange,array,pi,sin,cos,exp,log,log2,abs
from numpy import linalg as LA
from pylab import plot,xlabel,ylabel,show,ylim
import matplotlib.pyplot as plt
###################################################################

# Variables and functions of orbit
def f(var,beta):
    G = 4*pi**2
    Ms = 1.
    x = var[0]
    vx = var[1]
    y = var[2]
    vy = var[3]
    r = array([x, y], float)
    xprime = vx
    yprime = vy
    vxprime = -(G*Ms*x)/(LA.norm(r)**(beta + 1))
    vyprime = -(G*Ms*y)/(LA.norm(r)**(beta + 1))
    return array([xprime,vxprime,yprime,vyprime],float)


# Runge Kutta method for planetary motion
def orbit(dt,t,initialConds,beta):
    x = [initialConds[0]]
    vx = [initialConds[1]]
    y = [initialConds[2]]
    vy = [initialConds[3]]
    time = [0]
    varvec = initialConds
    while time[-1] < maxtime:
    	  k1 = dt * f(varvec,beta)
     	  k2 = dt * f(varvec + 0.5 * k1,beta)
          k3 = dt * f(varvec + 0.5 * k2,beta)
	  k4 = dt * f(varvec + k3,beta)
	  varvec +=  (k1 + 2 * k2 + 2 * k3 + k4)/6.
	  x.append(varvec[0])
	  vx.append(varvec[1])
          y.append(varvec[2])
          vy.append(varvec[3])
	  time.append(time[-1] + dt)  
    return time, x, vx, y, vy


# Variables of time
dt = .02
t = 0
maxtime = 4


# Initial conditions
initialx = 1.
initialvx = 0.
initialy = 0.
initialvy = 2*pi
initialbeta = input("Enter value of beta: ")
initials = array([initialx,initialvx,initialy,initialvy])
time,x,vx,y,vy = orbit(dt,t,initials,initialbeta)


# Plot of orbit as time progresses
for i,k in enumerate(x):
    plt.scatter(k,y[i])
    plt.xlim(-3,3)
    plt.ylim(-3,3)
    plt.axes().set_aspect('equal')
    plt.draw()
    plt.pause(.00000000000000000000001)
plt.show()
