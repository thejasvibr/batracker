#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Start to end prototype example
==============================




"""
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize']=10000
import numpy as np 
import scipy.io.wavfile as wavfile
from threshold_detector_prototyping import cross_channel_threshold_detector

# %% 
# Load the simulated audio file 
fs, audio = wavfile.read('simulated_audio/batracker_simple.wav')

# %% 
# Detect the calls in each channel 
detections = cross_channel_threshold_detector(audio, fs)

plt.figure()
ax= plt.subplot(411)
plt.specgram(audio[:,0], Fs=fs)
for each in detections[0]:
    plt.vlines(each, 0, fs*0.5, linewidth=0.2)

for i in range(2,5):
    plt.subplot(410+i, sharex=ax)
    plt.specgram(audio[:,i-1], Fs=fs)
    for each in detections[i-1]:
        plt.vlines(each, 0, fs*0.5, linewidth=0.2)


# %% 
# Perform correspondence matching and generate the common audio boundaries



# %% 
# Estimate time-difference-of-arrival across different channels and sounds



# %% 
# Use the TDOAs to calculate positions of sound sources






