#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correspondence matching prototyping
===================================
Implements the most basic type of correspondence matching -  which matches
signals across all channels to a single source emission. 

In this prototyping example, the simplest algorithm is used to find correspondences. 
The maximum inter-channel delay is calculated based on the longest dimension of the 
array, from the reference mic to the focal mic. 



Possible approaches
"""
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize'] = 10000
import numpy as np 
import scipy.io.wavfile as wavfile
import scipy.spatial.distance as distance

from threshold_detector_prototyping import cross_channel_threshold_detector

# %%
# Load the simulated audio file that was previously made

fs, audio = wavfile.read('simulated_audio/batracker_simple.wav')
short_audio = audio[:fs,:]
multi_detections = cross_channel_threshold_detector(short_audio, fs)

# %% 
# Also need some kind of visualiser to see the output of the all the detections 
# and diagnose problems.

num_channels = audio.shape[1]
all_axes = []
plt.figure(figsize=(10,8))
ax0 = plt.subplot(num_channels*100  + 11)
plt.specgram(short_audio[:, 0], Fs=fs, NFFT=256, noverlap=255)

for each in range(1,num_channels):
    plt.subplot(num_channels*100  + 10+each+1, sharex = ax0)
    plt.specgram(short_audio[:, each], Fs=fs, NFFT=256, noverlap=255)
    
    for every in multi_detections[each]:
        plt.vlines(every, 0, fs*0.5, linewidth=0.5)

