# -*- coding: utf-8 -*-
"""GCC implementation

Author: Thejasvi Beleyur 
License: MIT License
"""
import matplotlib.pyplot as plt
import numpy as np 
import scipy.signal as signal 

def estimate_gcc(signal1, ref_ch):
    
    fft1 = np.fft.rfft(signal1)
    fft2 = np.fft.rfft(ref_ch)
    cross_spectrum = fft1*np.conjugate(fft2)
    gcc = cross_spectrum/(np.abs(cross_spectrum))
    
    ifft_gcc = np.fft.irfft(gcc)
    ifft_gcc = np.roll(ifft_gcc, int(ifft_gcc.size*0.5))
    return ifft_gcc
    

if __name__ == '__main__':
    
    toa_1 = 0.000
    toa_2 = 0.001
    fs = 192000
    signal1 = np.random.normal(0,1e-5,20000)
    durn = 0.003
    samples = int(durn*fs)
    t = np.linspace(0,durn,samples)
    
    
    signal1[int(toa_1*fs):int(toa_1*fs)+samples] += signal.chirp(t, 25000, t[-1], 90000)*0.25
    signal2 = np.random.normal(0,1e-5,20000)
    signal2[int(toa_2*fs):int(toa_2*fs)+samples] += signal.chirp(t, 25000, t[-1], 90000)*0.15
    
    gcc_raw = estimate_gcc(signal1, signal2)
    
        
    plt.figure()
    plt.plot(gcc_raw)
