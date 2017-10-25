#!/usr/bin/env python

##############################################################
#Ryan Flamm
#6.1
# This code is meant to simulate the motion of a wave on a string as well as a fourier transform of it
##############################################################

from matplotlib.pyplot import figure,scatter,show,draw,pause,clf,ylim,xlim
from numpy import exp,pi,linspace,copy,arange,real,imag,empty
from numpy.linalg import norm as LA
from pylab import imshow,gray,show,figure
from numpy.fft import rfft, irfft
 
#################
#functions
def position(dt,dx,c,n):
	r = c*dt/dx
	return r

def wave(x):
	k = 30.0
	x0 = .35
	ynot = exp(-k * (x-x0)**2)
	return ynot

def dft(samples):
	N = len(samples)
	c = zeros(N//2+1,complex)
	for k in range(N//2+1):
		for n in range(N):
			c[k] += y[n]*exp(-2j*pi*k*n/N)
	return c

##################
a = 0
b = 1
dx = .01
dt = .01
c = 1.0
N = int((b-a)/dx)
r = c*dt/dx
samplingTime = 10.5
samples = int(samplingTime/dt)
freq = linspace(0,N,samples/2+1)

x = linspace(a,b,N)
y = [wave(n) for n in x]
yp = copy(y) 
ypp = copy(yp)

tMax = 30
t = 0
while t < tMax:
	for i in range(1,N-2):
		y[i] = 2*(1-r**2)*yp[i]-ypp[i]+r**2*(yp[i+1]+yp[i-1])
	y[0] = y[1]
	y[N-2] = y[N-3]
	t += dt
	ypp = copy(yp)
	yp = copy(y)
	scatter(x,y)
	ylim(-3,3)
	draw()
	pause(0.0001)
	clf()

#figure(1)
#scatter(freq,abs(c))
#xlim(0,5)
#show()
