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
import pandas as pd
import scipy.io.wavfile as wavfile
import scipy.spatial as spatial

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

# %% 
# Matching detections to each other. 
# The sounds are matched together according to the furthest microphone distances. 
# eg. if the mics 1 and 4 are about 2 m apart, then we expect all sounds
# to arrive within :math:`\pm`6 ms of each other from on
array_geometry = pd.read_csv('simulated_audio/array_geom.csv')

def calc_intermic_distances(array_geom):
    '''
    Parameters
    ----------
    array_geom: pd.DataFrame
        With compulsory column names, 'x', 'y', and 'z'
        Each mic's coordinates must be placed in one row, 
        in the order matching the recorded channels of
        the audio file. 
    Returns
    -------
    mic2mic_distance : np.array
        Nmic x Nmix array, a distance matrix
    
    '''
    mics_xyz = array_geom[['x','y','z']].to_numpy()
    mic2mic_distance = spatial.distance_matrix(mics_xyz, mics_xyz)
    return mic2mic_distance

def largest_mic2mic_distance(array_geom):
    '''
    Provides the largest mic-2-mic distance in the array. 
    In the presence of multiple equal mic-2-mic distances, 
    provides the first one. 
    
    Parameters
    ----------
    array_geom: pd.DataFrame
        With compulsory column names, 'x', 'y', and 'z'
        Each mic's coordinates must be placed in one row, 
        in the order matching the recorded channels of
        the audio file. 
    
    Returns
    -------
    maxval_location : np.array
        Array with 2 ints, with the row and column indices. 
        This corresponds to the mic pair with the largest mic-2-mic distance.
    maxval : float
        The largest mic-2-mic distance found. 
    
    See Also
    --------
    calc_intermic_distances
    
    '''
    mic2mic_distances = calc_intermic_distances(array_geom)
    maxval = np.max(mic2mic_distances)
    maxval_location = np.argwhere(maxval==mic2mic_distances)[0]
    return maxval_location, maxval
    
loc, max_dist  = largest_mic2mic_distance(array_geometry)

# %% 
# Now, find all sounds within +/- max delay distance of that mic pair, and
# put them into one boundary region for cross-correlation

