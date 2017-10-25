#!/usr/bin/env python

##############################################################
#Ryan Flamm
#5.10

#This code explores the over-relaxation method with regards the a
#point charge in empty space
##############################################################

from mpl_toolkits.mplot3d import Axes3D
from numpy import empty,zeros,max,pi,meshgrid,arange
from numpy.linalg import norm as LA
from pylab import imshow,gray,show,figure


def rho(x,y,xprime,yprime):
   qx1 = xprime + 2
   qx2 = xprime - 2
   qy1 = yprime + 2
   qy2 = yprime - 2
#   print x, y, 'current point'
#   print xprime, yprime, 'location of point charge'
   if qx2 < x < qx1 and qy2 < y < qy1:
#      print 'found point charge'
#      print x, y
      Q = 2
   else:
      Q = 0
   
   return Q

V = 10. #voltage
M = 100 # number of grid squares on a side
target = 10**(-6)
omega = .9
xprime = 60
yprime = 50



#arrays
phi = zeros([M+1,M+1],float)
phi[xprime,yprime] = V #position of point charge
#phi[0,:]=10  #voltage on wall


# Main Loop
keepGoing = True
while keepGoing:
   keepGoing = False
   for x in range(1,M):
      for y in range(1,M):
         before = phi[x,y]
         phi[x,y] = (1+omega)*(phi[x+1,y] + phi[x-1,y] + phi[x,y+1] + phi[x,y-1] + rho(x,y,xprime,yprime) )/4. - omega * phi[x,y]
         after = phi[x,y]

         if abs(after - before) > target:
            keepGoing = True
              #calc max difference from old values



x = arange(0,M+1,1)
y = arange(0,M+1,1)
x,y = meshgrid(x,y)

fig = figure(1)
ax = fig.gca(projection='3d')
ax.plot_surface(x,y,phi)

figure(2)
imshow(phi)
#gray()
show()
