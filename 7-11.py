#!/usr/bin/env python

##############################################################
#Ryan Flamm
#6.1
# This code is meant to simulate the diffussion of particles
#in space via probabilities.
##############################################################

from numpy import copy,array,meshgrid,arange
from pylab import show
import numpy as np
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D

#############################################################

#initial conditions
def f(x,y):
    if -4 < x < 4 and -4 < y < 4:
        return    2
   # if -10 < x < -6 and -10 < y < -6:
   #     return 10
   # if 5 < x < 9 and 5 < y < 9:
   #     return 10
    else:
        return    0




    
#step sizes for time and space
dt = 1
dx = 1
dy = 1
Dx = (1/6.)*(dx**2)
Dy = (1/6.)*(dy**2)

#grid boundaries and points
ax = -10
bx = 10
ay = -10
by = 10
nx = int((bx-ax)/dx) 
ny = int((by-ay)/dy)

#set x and y as well as the meshgrid fro the plot
x = arange(ax,bx,dx)
y = arange(ay,by,dy)
xP,yP = meshgrid(x,y) 

#defined variables that go into probability function
Cx = (Dx*dt)/(dx**2)
Cy = (Dy*dt)/(dy**2)


p = array([[f(r,l) for r in x] for l in y],float)
p2 = copy(p)

fig = plot.figure()

    #used to plot the initial conditions
#ax = fig.add_subplot(111, projection='3d')   
#ax.plot_wireframe(xP,yP,p)
#show()
#exit()


tMax = 100
t = 0

for t in range(1,tMax):
    p[1:nx-1,1:ny-1] = p2[1:nx-1,1:nx-1]+Cx*(p2[2:nx,1:ny-1]\
	+p2[0:nx-2,1:ny-1]-2*p2[1:nx-1,1:ny-1])+Cy*(p2[1:nx-1,2:ny]\
	+p2[1:nx-1,0:ny-2]-2*p2[1:nx-1,1:ny-1])

     #This form would only plot in one dimension.
    #for i in range(1,nx-1):                                         
        #p[i] = p2[i] + Cx * (p2[i+1] + p2[i-1] - 2*p2[i])
        #for j in range(1,ny-1):
            #p[j] = p2[j] + Cy * (p2[j+1] + p2[j-1] - 2*p2[j])
            
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(xP,yP,p)
    ax.set_zlim(0,2)
    plot.draw()
    plot.pause(1e-30)
    plot.clf()
    p2 = copy(p)
