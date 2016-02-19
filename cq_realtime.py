from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq

def cqRealTime(real_time_file,dipole_direction,kick_strength,damp_const):    
    '''
        (C) Joshua Goings 2016
        
        CQ_RealTime.py: a post-processing script for computing the absorption spectrum of
         Real Time Time Dependent SCF jobs in Chronus Quantum

        Computes the energy range, w (in eV) and dipole strength function S(w) for
         a given real time TD-SCF run. 

        real_time_file   ... type:string ; the RealTime_Dipole.csv file from a ChronusQ run
        dipole_direction ... type:char   ; which dipole moment contribution is computed (e.g. 'x','y', or 'z')
        kick_strength    ... type:float  ; in a.u., what was the applied field strength (e.g. 0.0001 au)
        damp_const       ... type:float  ; in a.u. of time, gives FWHM of 2/damp_const
    '''
   
    # chronusq file is CSV, also skip the header (first row)
    rt = np.genfromtxt(real_time_file,skip_header=1,delimiter=',')

    # choose which dipole axis you want
    if dipole_direction.lower() == 'x':
        direction = 2
    elif dipole_direction.lower() == 'y':
        direction = 3
    elif dipole_direction.lower() == 'z':
        direction = 4
    else:
        print "Not a valid direction for the dipole! Try: x,y,z "
        sys.exit(0)
    
    t      = rt[:,0] 

    # note 'z' is just generic dipole direction, converted from debye to au
    z      = rt[:,direction]*0.393456 

    # scale dipole signal  
    z0 = z[0]
    z = z - z0
    
    # add damping to give Lorenztian lineshape with FWHM of (2/damp_const)
    damp = np.exp(-(t-t[0])/damp_const)
    z = z * damp
    
    # pad signal with zeros. also not necessary, but makes spectra prettier
    #zero = np.linspace(0,0,10000)
    #z = np.hstack((z,zero))
   
    # do the fourier transform 
    fw = fft(z)

    # determine frequency range
    n = len(fw)                         # number samples, including padding
    timestep = t[1] - t[0]              # spacing between time samples; assumes constant time step
    w = fftfreq(n,d=timestep)*2.0*np.pi # frequency list
   
    fw_re = np.real(fw)                 # the real FFT frequencies
    fw_im = (np.imag(fw))               # the imaginary FFT frequencies
    fw_abs = abs(fw)                    # absolute value of frequencies
    
    # 'correct' equation for dipole strength function assuming you did SCF in static field
    #S = (2.0*w*w*fw_re)/(3.0*np.pi*137*kick_strength)

    # 'correct' equation for dipole strength function assuming you did a delta kick
    S = -(4.0*w*np.pi*fw_im)/(3.0*137*kick_strength)
    
    # you can print the integrated dipole strength function ... should equal number of electrons    
    #from scipy.integrate import simps
    #idx = np.where((w >= 0.0) & (w <=10000.0))
    #integral = simps(S[idx],w[idx]) 

    w = (w*27.2114)    # give frequencies in eV
    return w, S

if __name__ == '__main__':

    xFilename   = 'h2o_x_RealTime_Dipole.csv'
    yFilename   = 'h2o_y_RealTime_Dipole.csv'
    zFilename   = 'h2o_z_RealTime_Dipole.csv'
    
    kick        = 0.0001 # depends on system
    damping     = 150.0  # anywhere between 50-250 usually works well
    
    w, Sxx      = cqRealTime(xFilename,'x',kick,damping)
    w, Syy      = cqRealTime(yFilename,'y',kick,damping)
    w, Szz      = cqRealTime(zFilename,'z',kick,damping)
    
    plt.plot(w,Sxx+Syy+Szz,label='S')
    plt.ylim(0.0,2.0)  # y range
    plt.xlim(0,25)     # X range
    plt.legend()
    plt.show()
    #plt.savefig('myfile.pdf')

