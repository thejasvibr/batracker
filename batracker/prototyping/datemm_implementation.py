#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DATEMM implementation attempt
=============================
Scheuing and Yang 2008 [1] present the DATEMM algorithm which implements many of
the ideas I was only vaguely intuiting. Most importantly they look at two criteria:
    1) The auto-correlation of each channel tells us the 'lags' at which copies of the signal arrived. 
    2) The 'correct' set of TDOA's can be figured out by using the zero-cyclic sum condition between 
        any set of 3 microphones. That is :math:`\tau_{20} +\tau_{21}+\tau_{10} = 0 ` and so on, 
        for each 'sub-triangle'


References
~~~~~~~~~~
[1] Scheuing, J., & Yang, B. (2008). Disambiguation of TDOA estimation for multiple sources in reverberant environments.
    IEEE transactions on audio, speech, and language processing, 16(8), 1479-1489.


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

#%% 


small_audioseg = audio_segments[3][int(0.12*fs):int(0.140*fs):,:]
ch0 = small_audioseg[:,0]
chx = small_audioseg[:,3]

plt.figure()
plt.specgram(chx, Fs=fs)

#%%
acc0 = signal.correlate(ch0,ch0,'full')
acc1 = signal.correlate(chx,chx,'full')

min_timefluc = 1e-4
highest_freq = 1/min_timefluc
nyquist = fs*0.5
b,a = signal.butter(1, highest_freq/nyquist, 'lowpass')

acc0_lp = signal.filtfilt(b,a,acc0**2)
acc1_lp = signal.filtfilt(b,a,acc1**2)

dB_above = 60
thresh_ac0 = np.percentile(acc0_lp,10)*10**(dB_above/20.0)
thresh_ac1 = np.percentile(acc1_lp,10)*10**(dB_above/20.0)

cc = signal.correlate(chx,ch0,'full')
#%%
plt.figure()
#plt.plot(acc0_lp)
plt.plot(acc1_lp)
plt.hlines(thresh_ac1, 0, acc1_lp.size)


#%%
plt.figure()
a1 = plt.subplot(221)
plt.specgram(ch0,Fs=fs)
a2 = plt.subplot(222)
t = np.linspace(0,acc1_lp.size/fs, acc1_lp.size)
t -= np.mean(t)
plt.plot(t, acc0_lp)
plt.hlines(thresh_ac1, t[0], t[-1])
plt.ylim(0, thresh_ac1*2)
plt.subplot(223, sharey=a1, sharex=a1)
plt.specgram(chx,Fs=fs)
a4 = plt.subplot(224, sharex=a2)
t = np.linspace(0,acc1_lp.size/fs, acc1_lp.size)
t -= np.mean(t)
plt.plot(t, acc1_lp)
plt.hlines(thresh_ac1, t[0], t[-1])
plt.ylim(0, thresh_ac1*2)
