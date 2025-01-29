#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:39:39 2024

@author: lorrainemarcelin
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp,pi,sqrt
from matplotlib.animation import FFMpegWriter


def BrownianParticles(N,R,L,TotalTime,dt,Skip):
    #Function modeling dynamics Brownian particles given the parameters:
    #N = number of steps
    # R = radius of particle
    # L = size of box determining periodic boundary
    #parameters
    kBT = 0.004 #Thermal energy in pN*um
    eta = 0.001 #viscocity of water
    zeta = 6*pi*eta*R #drag force
    c = sqrt(zeta*kBT/dt)# magnitude of stochastic force
    e = 0.1 #depth of potential well
    sigma = 2 * R #effective size of particles
    
    Steps = round(TotalTime/dt/Skip)
    print(Steps)
    
    #initializing lists to stor positions and time
    x = np.zeros((2,N,Steps))
    Time = [Skip *dt*i for i in range(Steps)]
    
    #Randomly scattering particles in box for initial positions
    #box of size L
    x[:,:,0]= L *np.random.rand(2,N)
    
    
    
    #initialize matrix to store the particle-particle interactions
    Fi= np.zeros((2,N))
    
    metadata=dict(title='BrownianParticleHW',artist='Matplotlib')
    writer=FFMpegWriter(fps=15,codec='h264',bitrate=3000,metadata=metadata)
    fig1 = plt.figure()
    plt.xlim(0,L)
    plt.ylim(0,L)
    
    data1,= plt.plot([],[],'b')

    with writer.saving(fig1,'BrownianParticleHW.mp4',100):
        for i in range(1,Steps):
            x[:,:,i] = x[:,:,i-1]
            for j in range(Skip):
                # compute the random force
                xi = c*np.random.randn(2,N)
                Dx = np.meshgrid(x[0,:,i])
                Dx = Dx - np.transpose(Dx)
                Dy = np.meshgrid(x[1,:,i])
                Dy = Dy - np.transpose(Dy)
                Dx[Dx<-L/2] = Dx[Dx<-L/2] + L
                Dx[Dx>L/2] = Dx[Dx>L/2] - L
                Dy[Dy<-L/2] = Dy[Dy<-L/2] + L
                Dy[Dy>L/2] = Dy[Dy>L/2] - L
                r = np.sqrt(Dx**2 + Dy**2)
                for k in range(N):
                    r[k,k] = 10
                    FiMag = (12*e/sigma**2)*( (sigma/r)**(14) - (sigma/r)**8 )
                    FiMag[FiMag>100] = 100
                    
                    Fi[0,:] = -np.sum(FiMag*Dx,axis=1)
                    Fi[1,:] = -np.sum(FiMag*Dy,axis=1)
                    x[:,:,i] = x[:,:,i] + (dt/zeta)*(xi + Fi)
            
                # impose periodic boundary condition
                Mask = x[0,:,i]>L
                x[0,Mask,i] = x[0,Mask,i] - L
                Mask = x[0,:,i]<0
                x[0,Mask,i] = x[0,Mask,i] + L
                Mask = x[1,:,i]>L
                x[1,Mask,i] = x[1,Mask,i] - L
                Mask = x[1,:,i]<0
                x[1,Mask,i] = x[1,Mask,i] + L
            
            #plotmovie frame
            data1.set_data(x[0, :, i],x[1, :, i])
    
            
            writer.grab_frame()
            plt.pause(0.02)
            
    plt.show()
    plt.scatter(x[0, :, -1], x[1, :, -1], s=2)

    plt.show()
    return x,Time

#Test run
BrownianParticles(500, 1,1000, 2, 0.01, 20)





            

