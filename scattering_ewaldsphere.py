#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 10:30:45 2017

@author: ipoudyal
"""

import numpy as np
import center_mol as cm
import timeit
import scipy.constants
import sys

start = timeit.default_timer()

CMSF = np.empty((20,9))
CMSF[:] = np.NAN


#Cromer Mann Scattering factors 
#Uses nine parameter for each element
#Table 6.1.1.4. International Table of Crystallography Vol C. page 578
CMSF[5,:]  = [2.310,1.020,1.589,0.865,20.844,10.208,0.569,51.651,0.216]  #C a1,a2,a3,a4,b1,b2,b3,b4,c
CMSF[6,:]  = [12.213,3.132,2.013,1.166,0.006,9.893,28.997,0.583,-11.529] # N
CMSF[7,:]  = [3.049,2.287,1.546,0.867,13.277,5.701,0.324,32.909,0.251] # O
CMSF[12,:] = [6.420,1.900,1.594,1.965,3.039,0.7426,31.547,85.088,1.115] # Al
CMSF[14,:] = [6.435,4.179,1.780,1.491,1.907,27.157,0.526,68.165,1.115] # P
CMSF[15,:] = [6.905,5.203,1.438,1.586,1.468,22.215,0.254,56.172,0.867] # S

def scattering_factor(s,z,B):
    #Calculates electron scattering factor of an atom as a function of q
    # z is the Atomic number of the atom
    # q is the length of the scattering vector in A^-1.
    #Cromer-Mann Scattering factors is
    # f(s) = a1*exp(-b1*s^2/4) + a2*exp(-b2*s^2/4) + a3*exp(-b3*s^2/4) + a4*exp(-b4*s^2/4) + c 
    #s^2/4 = (sin(theta)/lambda)^2
    #2dsin(theta) = lambd; 1/d = q = 2Sin(theta)/lambda ; s=q/2=sin(theta)/lambda
    # Bernhard Rupp page 259
    #the coefficient above in CMSF are a1,a2,a3,a4,b1,b2,b3,b4 and c.
    
    a = CMSF[z-1,:]
    f = a[8] * np.ones((s.shape)) #include the constant coefficient c during initialization, see equation above
    for i in range(4):
        f = f + a[i]*np.exp(-a[i+4]*s**2/4.)  
    if B != 0 :
        f = f*np.exp(-B*s**2/4.)
    return f
    
def Flux_per_m2(photons,beamsize):
    ##photons = Number of photons
    ##beamsize = beam-size in microns
    beam_size = beamsize*1e-6 # convert microns to meters
    beam = beam_size*beam_size
    Flux = photons/beam
    return Flux
    

D = cm.length #Dimension of Box is D*oversampling ratio
d = 10.0 #Resolution is 2.0 Angstroms
lambd_A = 5.0#lambd*1e10 #Lambda in Angstroms
theta_max = 2*np.arcsin(lambd_A/(2*d)) #scattering_angle
qmax = np.tan(theta_max)/lambd_A  # if lambd_A=0 then qmax = 1./d
o = 2 #oversampling ratio
N = int(round(D*qmax*o))
M = N
L = N
#H,K,L = np.mgrid[-N:N+1,-M:M+1,-L:L+1]

qX, qY = np.mgrid[-N:N+1,-M:M+1]  # pixel coordinates of the detector

a = D*o     #Box size along first dimension
b = D*o
c = D*o

X = cm.X
Y = cm.Y 
Z = cm.Z
ATOM = cm.Type

qX = qX/float(a)   # pixel coordinates expressed in terms of q-space.
qY = qY/float(b) 
   

#Flux = 10^12 photons /(0.1 micron * 0.1 microns) #10^12 photons focused to 100*100nm^2.
photons = 1e12
beamsize = 0.1 # in microns
Flux = Flux_per_m2(photons,beamsize) # number per m**2
r_e = 2.8179e-15 #classical electron radius in meters

qabs = (2./lambd_A)*np.sin(0.5*np.arctan(lambd_A*np.sqrt((qX)**2 + (qY)**2)))
qx = qabs*np.sqrt(1-lambd_A**2*qabs**2/4.)*np.sin(np.arctan2(qX,qY)) 
qy = qabs*np.sqrt(1-lambd_A**2*qabs**2/4.)*np.cos(np.arctan2(qX,qY))   
qz = -0.5*lambd_A*qabs**2    
  
THETA = 2*np.arcsin(lambd_A*qabs/2.)

print("qX-shape:",qX.shape)
print("Size of molecule:",D)
print('lambd_A: Angstroms',lambd_A)
print("Flux per m^2:%.2E" %Flux)
print(" ")
print("scattering_angle:",theta_max*180./np.pi)
print("qabs_max,qabs_min",qabs.max(),qabs.min())   

qabs = qabs.reshape(qabs.size,1) 
qx = qx.reshape(qx.size,1)   
qy = qy.reshape(qy.size,1)
qz = qz.reshape(qz.size,1)
THETA = THETA.reshape(THETA.size,1)

SFC = scattering_factor(qabs,6,0.0)
SFN = scattering_factor(qabs,7,0.0)
SFO = scattering_factor(qabs,8,0.0)
SFS = scattering_factor(qabs,16,0.0)

print('SFC_shape',SFC.shape)

#Initialization to calculate phases
PhC = np.zeros((qabs.shape))
PhN = np.zeros((qabs.shape))
PhO = np.zeros((qabs.shape))
PhS = np.zeros((qabs.shape))


num_Atoms = 0
num_Electrons = 0
nC = 0
nN = 0
nO = 0
nS = 0

for j in range(len(ATOM)):
    if ATOM[j] in ['C']:
        PhC = PhC + np.exp(2*np.pi*1j*(qx*X[j] + qy*Y[j] + qz*Z[j]))
        num_Atoms = num_Atoms + 1
        num_Electrons = num_Electrons + 6
        nC = nC+1
    elif ATOM[j] in ['N']:
        PhN = PhN + np.exp(2*np.pi*1j*(qx*X[j] + qy*Y[j] + qz*Z[j]))
        num_Atoms = num_Atoms + 1
        num_Electrons = num_Electrons + 7  
        nN = nN + 1
    elif ATOM[j] in ['O']:
        PhO = PhO + np.exp(2*np.pi*1j*(qx*X[j] + qy*Y[j] + qz*Z[j]))
        num_Atoms = num_Atoms + 1
        num_Electrons = num_Electrons + 8
        nO = nO + 1
    elif ATOM[j] in ['S']:
        PhS = PhS + np.exp(2*np.pi*1j*(qx*X[j] + qy*Y[j] + qz*Z[j]))
        num_Atoms = num_Atoms + 1
        num_Electrons = num_Electrons + 16
        nS = nS +1


print("Total number of Atoms:",num_Atoms)
print("Total number of Electrons:",num_Electrons) 
print("Number of C,N,O,S:",(nC,nN,nO,nS))       
F = SFC*PhC + SFN*PhN + SFO*PhO + SFS*PhS
F2 = abs(F)**2
I = F2*Flux*(r_e**2)*(lambd_A/(o*D))**2  #Multiply intensity by prefactors
#I  = I *(np.cos(THETA))**2 #this is used for large scattering angle, for large scattering angle, photons recieved by pixels in the edge reduced. 
m = int((I.shape[0])**0.5)
I = I.reshape(m,m)
qabs = qabs.reshape(m,m)
print("I-shape",I.shape)
#indices = np.where(np.logical_and(qabs>=1./11,qabs<=1./10))
#indices = np.where(np.logical_and(qabs>=1./3.5,qabs<=1./3))
#nphotons = np.mean(I[indices])
#print("Number of photons in highest resolution shell:",nphotons)
stop = timeit.default_timer()

print('time-taken=%.3f mins'%((stop-start)/60.))

#import matplotlib.pyplot as plt
#I = np.minimum(I,I.max()/100.)
##I = np.minimum(I,50.0) # for noise plot
#plt.imshow(I)
##plt.imshow(qabs)
##plt.imshow(I,alpha=10.0)
#plt.colorbar()
#plt.show()
