#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:22:12 2024

@author: lorrainemarcelin
"""

from matplotlib.animation import FFMpegWriter
def Pendulum_Oscillations(beta,wD,TotalTime,dt,Skip):
    #Function modeling motion of a driven, damped harmonic pendulum 
    #with large displacements
    #beta is the dampening  constant

    # wD is the magnitude of the angular displacement
    #N is total number of steps
    
    #defining parameters
    g = 9.81 #gravitational constant g
    L=2 #length of rod connecting pendulum to pivot point
    w0 = np.sqrt(g/L) #natural frequency of pendulm
    A = 3 # intial amplitude of driving force
    
    
    # define step size
    Steps = round(TotalTime/dt/Skip)
    
    #initialize angular posiiton and angular velocity
    theta = np.zeros((Steps),'float')
    v= np.zeros((Steps),'float')
    
    #initialize time array
    Time = np.zeros((Steps),'float')

    
    # setup movie plotting stuff
    metadata=dict(title='Final_Project',artist='Matplotlib')
    writer=FFMpegWriter(fps=15,codec='h264',bitrate=3000,metadata=metadata)
    fig1 = plt.figure()
    l, = plt.plot([],[],'b-')
    plt.xlim(0,TotalTime)
    plt.ylim(-2*np.pi,2*np.pi)
    
    
    with writer.saving(fig1,"Final_Project.mp4",100):
    
        for i in range(1,Steps):
        
            theta[i] = theta[i-1]
            v[i] = v[i-1]
        
            Time[i] = Time[i-1]
        
            for j in range(Skip):
            
                #computing theta at the half time step
                thetamid = theta[i] +  dt*v[i]/2
            
                # Computing magnitude of driving force at half timestep
                Fmag = A*np.cos(wD*(Time[i]+dt/2))
            
                #computing velocity at half time step   
                vmid = v[i] + dt*(-2*beta*v[i] - w0**2*np.sin(theta[i])
                                         + Fmag)
                # Updating full velocity using the half-step values
                v[i] = v[i] + dt * (-2 * beta * vmid - w0**2 * np.sin(thetamid) + Fmag)

                # Updating full angular position using the half-step velocity
                theta[i] = theta[i] + dt * vmid
                
                Time[i] = Time[i] + dt 
            
            
            l.set_data(Time[:i+1],theta[:i+1])
         
            writer.grab_frame()
            plt.pause(0.01)
        plt.ylabel('Amplitude')
        plt.xlabel('Time')
        
    plt.savefig('Final_Project_wD.png',bbox_inches='tight')  
    
        
        
    plot(Time,theta,'r')   
    plt.show()
    
    
  
    return Time,theta
'''
'''
beta_values=np.linspace(0.025, 0.5, 50)
wD=2.21
TotalTime=350
dt=0.01

def beta_Variation(wD,beta_beta_values):
    
    #Plots the maximum steady-state amplitude vs the driving frequency
    #takes a list of driving frequency values and a single beta value
    beta_Amplitudes=[]
    for i in beta_values: #Computing the function for each wD value
    
        beta = i          
        
        Time, theta,w0 = Pendulum_Oscillations_RK4(beta=beta, wD=wD, TotalTime=TotalTime, dt=dt)
        
        
        steady_state_theta = theta[int(250/dt):] #getting the steady state amplitude for this wD
        #amplitude = np.max(steady_state_theta)
        rms_amplitude = np.sqrt(np.mean(steady_state_theta**2)) #getting the root mean square of the amplitude
        beta_Amplitudes.append(rms_amplitude)  #finding maximum amplitude and appending it to list
        
    plt.figure()
    
    plt.xlabel('Dampening Constant')
    plt.ylabel('RMS Amplitude (radians)')
    plt.legend()
    
    plt.plot(beta_values,beta_Amplitudes,'g')
    plt.savefig('Final_Project_Beta_RMSresonance.png',bbox_inches='tight')    
    plt.show()
    
    print(beta_Amplitudes)
    return beta_Amplitudes

beta_Variation(wD,beta_values)
'''