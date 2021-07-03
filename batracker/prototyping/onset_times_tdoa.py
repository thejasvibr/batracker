#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Another approach to estimating TDOAs is by looking at the onset times 
of a sound across channels. 

Let's try this out for botht he simulated and real audio. 
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
from tde_reverbenvironments import pick_top_peaks
#%% 
real_audio, fs = sf.read('SPKRPLAYBACK_multichirp_2018-08-18_09-15-06_ch9-12_first2s.wav')
sim_audio, fs = sf.read('tristar120_roomsimulation.wav')

tristar_geom = make_tristar_geom()

# process the real audio a bit - 
sync = real_audio[:,0]
sync[sync<0] = 0
cam_first_frame = np.argwhere(sync>=np.percentile(sync, 95)).flatten()[0]

array_audio = real_audio[cam_first_frame:,1:]
num_samples = 2*fs - array_audio.shape[0] 
array_audio = np.row_stack((array_audio, np.zeros((num_samples,4))))
sample_ind = np.arange(array_audio.shape[0])
segment_inds = np.array_split(sample_ind, 9)

audio_segments = [ array_audio[each,:] for each in segment_inds]

plt.figure()
plt.specgram(audio_segments[-1][:,0])
#%%
# Now try to see how the time of onset measurements help in good TDE measurement

def filtfilt_multichannel_audio(audio, b,a):
    '''
    '''
    bp_audio = np.apply_along_axis(lambda X: signal.filtfilt(b,a,X), 0, audio)
    return bp_audio

def multi_band_pass(audio,fs,**kwargs):
    '''
    Creates versions of the same audio filtered at multiple 
    'bands of interest'. By default splits the 0-Nyquist 
    into 4 equal bands. 
    '''
    band_centres = np.linspace(0,fs*0.5,5)
    bands = [ (each, band_centres[i+1]) for i,each in enumerate(band_centres[:-1])]
    bands_of_interest = kwargs.get('bands_of_interest', bands)
    
    bandpass_versions = []
    
    for each_band in bands_of_interest:
        peak_freq = np.mean(each_band)
        bw = np.abs(np.diff(each_band))
        Q = peak_freq/bw
        b,a = signal.iirpeak(peak_freq, Q, fs)
        bp_audio = filtfilt_multichannel_audio(audio, b,a)
        bandpass_versions.append(bp_audio)
        
    return bandpass_versions
        
#%%
boi = [(20000,40000),
       (50000,70000)]
multi_bp = multi_band_pass(audio_segments[1],fs, bands_of_interest=boi)

# Now calculate the onset time with the hilbert envelope for each version of the 
# bp audio 

dB = lambda X: 20*np.log10(np.abs(X))

def detect_onset_time(X,fs,**kwargs):
    envelope = np.abs(signal.hilbert(X))
    baseline = np.percentile(envelope, kwargs.get('baseline', 5)) # what is the 'silence'
    db_baseline = dB(baseline)
    db_envelope = dB(envelope)
    threshold = kwargs.get('threshold_dBre_baseline', 20)
    db_threshold = db_baseline + threshold
    geq_threshold = db_envelope >= db_threshold
    
    onset_sample = np.argwhere(geq_threshold).flatten()[0]
    onset_time = onset_sample/fs
    return onset_time
#%%
# The onset times are pretty blah actually -- not all of them match up. 
vsound = 338.0
for i,_ in enumerate(boi):
    onset_times = np.apply_along_axis(detect_onset_time, 0, multi_bp[i],
                                      fs, **{'threshold_dBre_baseline':36})
    tdoas = onset_times[1:] - onset_times[0]
    d_potential = np.array(tdoas)*vsound
    pos1, pos2 = sw02.spiesberger_wahlberg_solution(tristar_geom.T, d_potential)
    if pos1[1]>=0:
        print(pos1)
    elif pos2[1]>=0:
        print(pos2)
#%% 
def choose_positive_y(pos1,pos2):
    if pos1[1]>=0:
        return pos1
    elif pos2[1]>=0:
        return pos2
    else:
        return [np.nan, np.nan, np.nan]
#%% What if we do the gcc style estimation for each bp segment? 
band_tde_combinations = {}
for i,each in enumerate(multi_bp):
    gcc_m10 = estimate_gcc(each[:,1], each[:,0])
    gcc_m20 = estimate_gcc(each[:,2], each[:,0])
    gcc_m30 = estimate_gcc(each[:,3], each[:,0])
    
    tde_m10 = (pick_top_peaks(gcc_m10,3) - gcc_m10.size/2.0)/fs
    tde_m20 = (pick_top_peaks(gcc_m20,3) - gcc_m20.size/2.0)/fs
    tde_m30 = (pick_top_peaks(gcc_m30,3) - gcc_m30.size/2.0)/fs

    # pick top 2 peaks from ech gcc
    tdoa_combinations = []
    for td10 in tde_m10:
        for td20 in tde_m20:
            for td30 in tde_m30:
                tdoa_combinations.append(np.array([td10, td20, td30]))
    band_tde_combinations[i] = tdoa_combinations
#%%
# calculate implied positions from the tdoa combinations per band

band_positions = {}
for band, bp_tdoas in band_tde_combinations.items():
    this_band_positions = []
    for each in bp_tdoas:
        d_potential = vsound* each
        pos1, pos2 = sw02.spiesberger_wahlberg_solution(tristar_geom.T, d_potential)
        this_band_positions.append(choose_positive_y(pos1, pos2))
    band_positions[band] = this_band_positions


#%% Using
    


#%% Take the ACC of each channel:
small_audioseg = audio_segments[0][int(0.18*fs):,:]

acc0 = signal.correlate(small_audioseg[:,0],small_audioseg[:,0],'full')
acc1 = signal.correlate(small_audioseg[:,1],small_audioseg[:,1],'full')
cc = signal.correlate(small_audioseg[:,1],small_audioseg[:,0],'full')

plt.figure()
a0 = plt.subplot(311)
plt.plot(cc)
plt.subplot(312,sharey=a0,sharex=a0)
plt.plot(acc0, label='channel 1')
plt.legend()
plt.subplot(313,sharey=a0,sharex=a0)
plt.plot(acc1, label='channel 2')
plt.legend()

plt.figure()
a1 = plt.subplot(211)
plt.specgram(small_audioseg[:,0],Fs=fs)
plt.subplot(212, sharey=a1, sharex=a1)
plt.specgram(small_audioseg[:,1],Fs=fs)