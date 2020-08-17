#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Start to end prototype example
==============================




"""
import batracker
from batracker.localisation import schau_robinson_1987 as sr1987
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize']=10000
import numpy as np 
import pandas as pd
import scipy.io.wavfile as wavfile
from plot_threshold_detector import cross_channel_threshold_detector
from plot_correspondence_matching import generate_crosscor_boundaries
from plot_tdoa_prototyper import measure_tdoa

# %% 
# Load the simulated audio file 
# THE TACOST WORKFLOW THROUGH THE CLI IS BROKEN -- VERIFY !!! 
# fs, full_audio = wavfile.read('simulated_audio/batracker_simple.wav')
# audio = full_audio[:int(fs*0.75),:] # keep the whole example short first!

import tacost
from tacost import calculate_toa as ctoa

# %% 
# Generate an array geometry and use previously made  source positions

mic_positions = np.array([[0,0,1],
                       [1,0,0],
                       [-1,0,0],
                       [0,1,0]
                       ])


test_pos = np.random.choice(np.arange(-2,5,0.5),24).reshape(-1,3)
toa = ctoa.calculate_mic_arrival_times(test_pos,
                                       array_geometry=mic_positions)

test_audio, fs = tacost.make_sim_audio.create_audio_for_mic_arrival_times(toa,
                                                                     call_snr=80)
audio = test_audio.copy()

# %% 
# Detect the calls in each channel 
detections = cross_channel_threshold_detector(audio, fs,
                                              dbrms_window=0.5*10**-3,
                                              dbrms_threshold=-60)


# Spectrogram of the cross-corr boundaries
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


# Waveformsof the detection 

plt.figure()
ax= plt.subplot(411)
plt.plot(audio[:,0])
for each in detections[0]:
    plt.vlines(np.array(each)*fs, np.min(audio), np.max(audio), 'k',linewidth=0.5)

for i in range(2,5):
    plt.subplot(410+i, sharex=ax)
    plt.plot(audio[:,i-1])
    for each in detections[i-1]:
        plt.vlines(np.array(each)*fs,
                   np.min(audio), np.max(audio), 'k',linewidth=0.5)

# %% 
# The cross-cor boundary needs to be calculated keeping the whole signal duration 
# in mind too!!

# %% 
# Perform correspondence matching and generate the common boundaries
# across channels for cross=correlation. Also, load the mic array geometry 
# as the max- inter-mic distances are required  for the calculation of max
# inter-mic delays

ag = pd.DataFrame(mic_positions)
ag.columns  = ['x','y','z']

crosscor_boundaries = generate_crosscor_boundaries(detections, ag)

num_channels = audio.shape[1]


# %% 
# Estimate time-difference-of-arrival across different channels and sounds
all_tdoas = {}
for i,each_common in enumerate(crosscor_boundaries):
    start, stop = each_common
    start_sample, stop_sample = int(start*fs), int(stop*fs)
    
    tdoas = measure_tdoa(audio[start_sample:stop_sample,:], fs)
    all_tdoas[i] = tdoas

# %% 
# Use the TDOAs to calculate positions of sound sources
vsound = 338.0
all_positions = []
for det_number, tdoas in all_tdoas.items():
    try:
        d = vsound*tdoas.reshape(3,1)
        pos = sr1987.schau_robinson_solution(mic_positions, d)
        all_positions.append(pos)
    except:
        print(f'Couldnt calculate for detection {det_number}')

# %% 
# The point I"m stuck at here is that the SR1987 gives two solutions for any
# set of of TDOAs, how to choose one of them without prior information. It seems
# like you'll always need some kind of range estimate to choose the one that 
# seems 'relevant' -- but this is exactly what you don't want to do in a situation
# where you'd like complete automation. In their paper, Schau and Robinson
# themselves seem to suggest this : `'In these situations, the two locations are
# usually far enough apart that the correct solution can be discerned by other 
# physical reasoning such as one solution lying outside the domain of interest.'`


print(all_positions)

print(test_pos)

