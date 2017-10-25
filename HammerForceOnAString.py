#!/usr/bin/env python

##################################################################

#Ryan Flamm
#This program studies the strike of a hammer on a piano string. The effects are studied through the force on the soundboard as well as the frequencies produced.

################################################################## 
import matplotlib.pyplot as plt
from random import randint,random,sample,uniform,choice
from numpy import array,linalg,where,sum,exp,sin,abs,pi,arctan,arange,sqrt,cos,linspace,sign,copy,zeros,empty
from numpy.linalg import norm
from numpy.fft import fft
##################################################################

# Bounds
a = 0 		#Left bound
b = 1.06	#Right bound
dx = .00065
r = 1.
c = 340.
dt = r/c * dx
N = int((b-a)/dx)
x = linspace(a,b,N)
y = zeros(len(x),float)
yp = copy(y) 
ypp = copy(yp)


# Initial conditions
tMax = 0.05
t = 0
T = 349		# Tension in string
mh = .0043	# Mass of the hammer
f = 65		# Frequency of the string
vh = 4.	        # Speed of the hammer. Large value means a stronger hammer strike
zh = 0
m = .008	# Linear density of the string
mew =.008
Q = b/8.	# Location of the hammer strick
QIndex = int(Q/dx)
Fh = zeros(len(x),float)


# Initialize arrays for finding the force on the bridge
Fb = [] 
time = arange(0,tMax,dt)
index = 0
for t in time:        
        p = 3
        K = 1e11
        zf = zh - y[QIndex]
        if zf < 0:
                Fh[QIndex-2:QIndex+2] = 0
        else:
                Fh[QIndex-2:QIndex+2] = K * zf**p
        y[1:N-1] = 2*(1-r**2)*yp[1:N-1]-ypp[1:N-1]+r**2 * (yp[2:N]+yp[0:N-2]) + ((dt**2)/(mew*dx))*Fh[1:N-1]
        vh -= Fh[QIndex]/mh *dt
        zh +=  vh * dt 
	Fb.append(T *((y[1]-y[0])/dx))

	
# This section will run the animation
        if index % 100 == 0:
		plt.figure(1)
	        plt.scatter(x,y)
		plt.xlim(-.2,1.2)
		plt.ylim(-.01,.01)
	        plt.draw()
	        plt.pause(1e-10)
	        plt.clf()
	ypp = copy(yp)
	yp = copy(y)
	index += 1

	
# plot the force
plt.figure(3)
plt.ylim(-10,20)
plt.scatter(time,Fb)


# Fourier Transform
d = fft(Fb)
fN = 1/(2* dt) 
freq = linspace(0,fN,len(d)/2)
plt.figure(2)
plt.scatter(freq,abs(d[:len(d)/2]))
plt.xlim(0,2000)
plt.show()
