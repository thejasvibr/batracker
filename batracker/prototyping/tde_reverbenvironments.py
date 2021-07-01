# -*- coding: utf-8 -*-
"""
Trying out TDE matching in an array in reverberant environments
---------------------------------------------------------------
This module will describe my efforts at working out how to perform TDE (time delay estimation)
in reverb environments using a simulated recording made with pyroomacoustics [1]. 

The 120cm tristar array is placed a little ahead of the wall and 10ms chirps are 
played from 3 locations. 

I'd like to see if the TDE pairs form 'curves' as described in Zhayida et al. 2016 [2].

Estimating TDEs
~~~~~~~~~~~~~~~
Instead of using cross-correlation as I've done in the past, I'll now use the generalised cross-correlation. 
In general, it's supposed to be more reverb resistant - but I ahve the feeling it also shows the 'peaks' for the 
delays a bit better?

References
~~~~~~~~~~
[1] Scheibler, R., Bezzam, E., & Dokmanić, I. (2018, April). Pyroomacoustics: A python package for audio room simulation and array processing algorithms. In 2018 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) (pp. 351-355). IEEE.
[2] Zhayida, S., Rex, S. S., Kuang, Y., Andersson, F., & Åström, K. (2016). An automatic system for acoustic microphone geometry calibration based on minimal solvers. arXiv preprint arXiv:1610.02392.
"""

import soundfile as sf 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np 
import scipy.signal as signal 
import scipy.spatial as spatial
from gen_cross_corr import estimate_gcc
from room_tdoa_simulation import make_source_positions,make_tristar_geom
from batracker.localisation import spiesberger_wahlberg_2002 as sw02


def pick_top_peaks(X, n_peaks=4):
    '''
    Parameters
    ----------
    X: np.array
    n_peaks : int>0, optional
        Number of peaks to be chosen. Defaults to 4. 
        
    Returns
    -------
    peak_inds : np.array
    '''
    # thanks https://stackoverflow.com/a/23734295/4955732
    peak_inds = np.argpartition(X, -n_peaks)[-n_peaks:]
    return peak_inds
    

dB = lambda X: 20*np.log10(np.abs(X))

# load the simulated audio
audio, fs = sf.read('tristar120_roomsimulation.wav')
n_pbks = 20
# select the signal portions for which the TDE must be performed. 
seg_samples = int(0.07*fs)
seg_inds = np.array_split(np.arange(audio.shape[0]), n_pbks)
segment_inds = [ slice(each[0], each[-1]) for each in seg_inds]


plt.figure()
a0 = plt.subplot(411)
plt.specgram(audio[:,0], Fs=fs)

for each in range(1,4):
    plt.subplot(411+each, sharex=a0, sharey=a0)
    plt.specgram(audio[:,each], Fs=fs)


#%% Trying for Call1
audio_segments = [audio[segment_inds[each],:] for each in range(n_pbks)]
tdes_m01 = []
tdes_m02 = []
tdes_m03 = []

for i,each in enumerate(audio_segments):
    tdes_m01.append(estimate_gcc(each[:,1], each[:,0]))
    tdes_m02.append(estimate_gcc(each[:,2], each[:,0]))
    tdes_m03.append(estimate_gcc(each[:,3], each[:,0]))

num_peaks = 10

tde_peaks_m01 = [pick_top_peaks(tdes_m01[each], num_peaks) for each in range(len(tdes_m01))]
tde_peaks_m02 = [pick_top_peaks(tdes_m02[each], num_peaks) for each in range(len(tdes_m02))]
tde_peaks_m03 = [pick_top_peaks(tdes_m03[each], num_peaks) for each in range(len(tdes_m03))]


#%% 
vsound = 338.0
source_positions = make_source_positions()
mic_positions = make_tristar_geom()
source_m0 = spatial.distance.euclidean(source_positions[0,:], mic_positions[:,0])
source_mics = []
for each in [1,2,3]:
    source_mics.append(spatial.distance.euclidean(source_positions[0,:], mic_positions[:,each]))

delta_R_m01 = np.array(source_mics) - source_m0
tau_m01 = delta_R_m01/vsound

tde_peaks_m01_time = [(each- tdes_m01[i].size/2.0)/fs for i,each in enumerate(tde_peaks_m01)]
tde_peaks_m02_time = [(each- tdes_m02[i].size/2.0)/fs for i,each in enumerate(tde_peaks_m02)]
tde_peaks_m03_time = [(each- tdes_m03[i].size/2.0)/fs for i,each in enumerate(tde_peaks_m03)]

#%% 
# What happens if we explore the positions resulting from all possible TDOAs between 
# N channels? Is there some way we can figure out which ones are 'true' positions? Here we have 64 potential
# localisations!!! 

position_number = 5
all_potential_tdes = [ each[position_number] for each in [tde_peaks_m01_time,tde_peaks_m02_time,tde_peaks_m03_time]]

potential_tdoas = []
for each_m01 in all_potential_tdes[0]:
    for each_m02 in all_potential_tdes[1]:
        for each_m03 in all_potential_tdes[2]:
            potential_tdoas.append([each_m01, each_m02, each_m03])
    
#%% 
# Now generate the positions that each tdoa vector implies. 
all_possible_positions = []
for each_tdoa in potential_tdoas:
    d_potential = np.array(each_tdoa)*vsound
    pos1, pos2 = sw02.spiesberger_wahlberg_solution(mic_positions.T, d_potential)
    # keep only those with +ve y solution 
    if not np.isnan(np.sum(pos1)):
        if np.logical_and(pos1[1]>=0,pos1[2]>=0):
            if pos1[0]>=0:
                all_possible_positions.append(pos1)
    if not np.isnan(np.sum(pos1)):
        if np.logical_and(pos2[1]>=0,pos2[2]>=0):
            if pos2[0]>=0:
                all_possible_positions.append(pos2)

candidate_positions = np.array(all_possible_positions).reshape(-1,3)

#%% 

f1 = plt.figure()
af1 = f1.add_subplot(111, projection='3d')
af1.plot(candidate_positions[:,0], candidate_positions[:,1], candidate_positions[:,2],'*')
af1.plot(source_positions[position_number,0], 
         source_positions[position_number,1], source_positions[position_number,2],'g*')

tristar_geom = mic_positions.copy()
m0 = tristar_geom[:,0]
for each in range(1,4):    
    line = np.column_stack((m0, tristar_geom[:,each]))
    x,y,z = tristar_geom[:,each]
    af1.text3D(x,y,z,'ch '+str(each+1))
    af1.plot(line[0,:], line[1,:], line[2,:], '-k')

af1.set_ylim(0,5)
af1.set_xlim(0,7)
af1.set_zlim(0,5)

# is the 'real' localisation actually even there at all?
abs_diff = np.abs(candidate_positions-source_positions[position_number,:])
abs_euc = np.sqrt(np.sum(abs_diff**2, 1))
index = np.argmin(abs_euc)
closest_match = candidate_positions[index,:]
af1.plot(closest_match[0], closest_match[1], closest_match[2],'r^')


print(f'calculated {candidate_positions[index,:]}, \n actual:{source_positions[position_number,:]} \n error:{abs_euc[index]}')

#%% How to choose the correct candidate position?
# My thoughts right now are around seeing if each candidate position predicts the 
# TOA for each channel properly. My intuition tells me the 'correct' position will
# minimise the TOA for all channels - as it tells us the most direct path taken.
# So the idea would be to calculate the TOAs for all candidate positions - and see which one 
# creates the lowest TOAs for all channels. 