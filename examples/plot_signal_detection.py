#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inbuilt Signal Detection Routines 
=================================

`INCOMPLETE EXAMPLE`

`batracker` has two inbuilt signal detection routines:
    
#. Envelope detector
    
    The envelope detector generates the Hilbert envelope of the input audio. 
    Any region that is :math:`\geq` a threshold in dB peak is considered a 
    valid signal.

#. dB rms profile detector
    
    The moving dBrms profile of the input audio is generated.  Regions :math:`\geq`
    a threshold are considered valid signal.

This example will go on to showcase how signal detection works using these two
routines.


The synthetic audio recording
-----------------------------
First let's make a bat-like echolocation call to test this with. The detection
routines should work regardless of the actual signal type as they only really
on how 'loud' the signal is recorded in the audio file. 

"""

import matplotlib.pyplot as plt 
plt.rcParams['agg.path.chunksize'] = 10000
import numpy as np 
import scipy.signal as signal

# Function which generates a bat call

def generate_5ms_bat_call():
    '''
    Returns
    -------
    bat_call
    fs
    duration
    
    '''
    duration = 0.005
    fs = 250000
    t = np.linspace(0, duration, int(fs*duration))
    bat_call = signal.chirp(t,90000, 25000, t[-1])
    bat_call *= 0.5
    return bat_call, fs, duration







