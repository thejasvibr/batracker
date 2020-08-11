#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Start to end prototype example
==============================




"""
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize']=10000
import numpy as np 
import pandas as pd
import scipy.io.wavfile as wavfile
from threshold_detector_prototyping import cross_channel_threshold_detector
from correspondence_match_prototyping import generate_crosscor_boundaries
# %% 
# Load the simulated audio file 
fs, full_audio = wavfile.read('simulated_audio/batracker_simple.wav')
audio = full_audio[:fs,:] # keep the whole example short first!


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
# Perform correspondence matching and generate the common boundaries
# across channels for cross=correlation. Also, load the mic array geometry 
# as the max- inter-mic distances are required  for the calculation of max
# inter-mic delays

array_geom = pd.read_csv('simulated_audio/array_geom.csv')

crosscor_boundaries = generate_crosscor_boundaries(detections,array_geom)

num_channels = audio.shape[1]
for each in range(1,num_channels):
        plt.subplot(num_channels*100  + 10+each+1, sharex = ax)    
        for every in crosscor_boundaries:
            plt.vlines(every, 0, fs*0.5,'r', linewidth=0.5, )

# %% 
# Estimate time-difference-of-arrival across different channels and sounds



# %% 
# Use the TDOAs to calculate positions of sound sources






