#!/usr/bin/env python

##################################################################
#Ryan Flamm
#9-13

#This code uses Molecular Dynamics to study the melting point of
#a 3D system.
##################################################################

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from random import uniform
from numpy import array,abs,arange,mod,copy,sign
from numpy.linalg import norm

##################################################################

#calculates net force due to all the particle interactions
def getForce(positions,particle,sizeX,sizeY,sizeZ):
    fX = 0
    fY = 0
    fZ = 0
    #find components of r vector
    for i in positions:
        xDiff = (i - particle)[0]
        yDiff = (i - particle)[1]
        zDiff = (i - particle)[2]

        #establishes periodic boundaries
        if abs(xDiff) > sizeX/2. and abs(yDiff) > sizeY/2. and abs(zDiff) > sizeZ/2.:
            modifiedPos = array([i[0] - sign(xDiff) * sizeX, i[1] - sign(yDiff) * sizeY, i[2] - sign(zDiff) * sizeZ])
        elif abs(xDiff) > sizeX/2. and abs(yDiff) > sizeY/2.:
            modifiedPos = array([i[0] - sign(xDiff) * sizeX, i[1] - sign(yDiff) * sizeY, i[2]])
        elif abs(xDiff) > sizeX/2. and abs(zDiff) > sizeZ/2.:
            modifiedPos = array([i[0] - sign(xDiff) * sizeX, i[1], i[2] - sign(zDiff) * sizeZ])
        elif abs(yDiff) > sizeY/2. and abs(zDiff) > sizeZ/2.:
            modifiedPos = array([i[0], i[1] - sign(yDiff) * sizeY, i[2] - sign(zDiff) * sizeZ])
        elif abs(xDiff) > sizeX/2.:
            modifiedPos = array([i[0] - sign(xDiff) * sizeX, i[1], i[2]])
        elif abs(yDiff) > sizeY/2.:
            modifiedPos = array([i[0], i[1] - sign(yDiff) * sizeY, i[2]])
        elif abs(zDiff) > sizeZ/2.:
            modifiedPos = array([i[0], i[1], i[2] - sign(zDiff) * sizeZ])
        else:
            modifiedPos = copy(i)

        #redefines components using the shortest interparticle distance
        xDiff = (particle - modifiedPos)[0]
        yDiff = (particle - modifiedPos)[1]
        zDiff = (particle - modifiedPos)[2]
        distance = norm(particle - modifiedPos)

        #establishes the strength of the force using Lennard-Jones potential
        #these bounds are used to cut unnecessary calculations above r=3
        if 0.1 <= distance <= 3.:
            fX += 24. * (3./distance**13 - 1./distance**7) * xDiff/distance
            fY += 24. * (3./distance**13 - 1./distance**7) * yDiff/distance
            fZ += 24. * (3./distance**13 - 1./distance**7) * zDiff/distance

    return fX,fY,fZ

#sets up the initial conditions of the particles
def initializeParticles(lVecs,sizeX,sizeY,sizeZ):
    positions = []
    velocities = []
    deltaR = 0.1
    v0 = 10            #initial velocity

    for i in range(0,sizeX):
        for j in range(0,sizeY):
            for k in range(0,sizeZ):
                positions.append(i * lVecs[0] + j * lVecs[1] + k * lVecs[2] + 0.5 + [2 * (uniform(0,1) - 0.5) * deltaR, 2 * (uniform(0,1) - 0.5) * deltaR, 2* (uniform(0,1) - 0.5) * deltaR])
                velocities.append([2 * (uniform(0,1) - 0.5) * v0, 2 * (uniform(0,1) - 0.5) * v0, 2 * (uniform(0,1) - 0.5) * v0])
    return array(positions),array(velocities)

sizeX = 4                             #size of box in x direction
sizeY = 4                             #size of box in y direction
sizeZ = 4                             #size of box in z direction
nParticles = sizeX * sizeY * sizeZ    #number of particles found in the box
time = 0
tMax = 16
dt = 0.001
index = 0
lVecs = array([[1.,0,0],[0,1.,0],[0,0,1.]])


positions,velocities = initializeParticles(lVecs,sizeX,sizeY,sizeZ)

#sets current conditions to that found in InitialParticles
positionsC = copy(positions)

#previous positions
positionsP = positions - velocities * dt


while time < tMax:
    if index % 1 == 0:
        plt.scatter([n[0] for n in positions],[n[1] for n in positions],s=46)
        fig = plt.figure(1)
        ax = fig.add_subplot(111,projection='3d')
        ax.set_xlim(0,sizeX)
        ax.set_ylim(0,sizeY)
        ax.set_zlim(0,sizeZ)
        ax.scatter([n[0] for n in positions],[n[1] for n in positions],[n[2] for n in positions])
        plt.draw()
        plt.pause(0.000001)
#        if time < 13:
#            plt.clf()

    forces = []

    for i in range(len(positions)):
        fX,fY,fZ = getForce(positions,positions[i],sizeX,sizeY,sizeZ)
        forces.append([fX,fY,fZ])

        #verlet method
    for i in range(len(positions)):
        prevX = positionsP[i][0]
        currentX = positionsC[i][0]
        prevY = positionsP[i][1]
        currentY = positionsC[i][1]
        prevZ = positionsP[i][2]
        currentZ = positionsC[i][2]
        fX = forces[i][0]
        fY = forces[i][1]
        fZ = forces[i][2]
        #update positions
        positions[i][0] = 2. * currentX - prevX + fX * dt**2
        positions[i][1] = 2. * currentY - prevY + fY * dt**2
        positions[i][2] = 2. * currentZ - prevZ + fZ * dt**2
        #update velocities
        velocities[i][0] = (positions[i][0] - prevX)/(2. * dt)
        velocities[i][1] = (positions[i][1] - prevY)/(2. * dt)
        velocities[i][2] = (positions[i][2] - prevZ)/(2. * dt)
        #periodic boundary conditions
        if positions[i][0] > sizeX:
            positions[i][0] = positions[i][0] - sizeX
        if positions[i][0] < 0:
            positions[i][0] = positions[i][0] + sizeX
        if positions[i][1] > sizeY:
            positions[i][1] = positions[i][1] - sizeY
        if positions[i][1] < 0:
            positions[i][1] = positions[i][1] + sizeY
        if positions[i][2] > sizeZ:
            positions[i][2] = positions[i][2] - sizeZ
        if positions[i][2] < 0:
            positions[i][2] = positions[i][2] + sizeZ

    time += dt
    index += 1
    #turn current positions into the new previous positions
    positionsP = copy(positionsC)
    #update current positions
    positionsC = copy(positions)
plt.show() 
