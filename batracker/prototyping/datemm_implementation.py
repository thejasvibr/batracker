#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DATEMM implementation attempt
=============================
Scheuing and Yang 2008 [1] present the DATEMM algorithm which implements many of
the ideas I was only vaguely intuiting. Most importantly they look at two criteria:
    1) The auto-correlation of each channel tells us the 'lags' at which copies of the signal arrived. 
        * The last peak in the auto-correlation tells us the first time of arrival of the signal 
        * These 'extremum' positions can then inform sensible TDE's between microphone pairs
        * A general set of potential TDOA's can thus be generated
    2) The 'valid' set of TDOA's can be figured out by using the zero-cyclic sum condition between 
        any set of 3 microphones. That is :math:`\tau_{20} +\tau_{21}+\tau_{10} = 0 ` and so on, 
        for each 'sub-triangle'
    3) 
    


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


#%% 

def generate_multich_gcc(audio):
    '''
    '''
    multich_gcc = {}
    nchannels = audio.shape[1]
    
    for i in range(nchannels):
        for j in range(nchannels):
            multich_gcc[i,j] = estimate_gcc(audio[:,i], audio[:,j])
    return multich_gcc


#%%
# Pick up TDE's across all channel pairs

def generate_multich_tdes(multichannel_gcc, fs,**kwargs):
    '''
    Parameters
    ----------
    multichannel_gcc : dictionary
        With channel pair (i,j) as key and GCC np.array as entry.
    fs : float>0
        Sampling rate
    n_peaks : int>0, optional 
        Defaults to 4 peaks. 
    
    Returns
    -------
    multich_tdes : dictionary
        channel pair (i,j) as key and the top `N` time delay peaks as
        np.array entry.

    See Also
    --------
    pick_top_peaks
    '''
    multich_tdes = {}
    for channel_pair, gcc in multichannel_gcc.items():
        peaks = pick_top_peaks(gcc, **kwargs) # in samples
        peaks_seconds = (peaks - gcc.size/2.0)/fs
        multich_tdes[channel_pair] = peaks_seconds
    return multich_tdes
    

#%%
# Generate all autocorrrelations

def generate_multich_autocorrelation(audio, fs):
    '''
    Parameters
    ----------
    audio : np.array
        Multichannel audio 
    fs : float>0
        sampling rate
    Returns 
    -------
    acc_matrix: np.array
        Array with autocorrelation for each channel (Msamples x Nchannels)
    '''
    acc_matrix = np.zeros(audio.shape)
    for i in range(audio.shape[1]):
        acc_matrix[:,i] = signal.correlate(audio[:,i], audio[:,i],'same')
    return acc_matrix


#%% Refine autocorrelations and detect peaks. 

def detect_autocorrelation_peaks(acc_mat, fs, **kwargs):
    '''
    Peak detection is done by squaring the autocorrelation signal 
    and lowpass filtering it to remove small time-scale fluctuations. 
    
    Peaks are detected for all parts of the 

    Parameters
    ----------
    acc_mat : np.array
        Msamples x Nchannels with autocorrelation of each audio channel 
    fs : float>)
    time_resolution : float>0, optional
        The time resolution of the autocorrelation peak detection. 
        Defaults to 0.1 milliseconds
    min_duration : float>0, optional 
        The minimum duration of a sound. Defaults to 3 milliseconds
    
    Returns 
    -------
    
    '''


if __name__ == '__main__':
    #%%
    fname = ['SPKRPLAYBACK_multichirp_2018-08-18_09-15-06_ch9-12_first2s.wav', 
             'tristar120_roomsimulation.wav']
    sim_audio, fs = sf.read(fname[1])
    real_audio, fs = sf.read(fname[0])
    
    
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
    
    small_audioseg = audio_segments[3][int(0.12*fs):int(0.140*fs):,:]
    
    #%%
    
    # generate multichannel GCC
    mch_gcc = generate_multich_gcc(small_audioseg)
    # choose top X peaks of time delay     
    mch_tdes = generate_multich_tdes(mch_gcc, fs, n_peaks=10) 

    # Autocorrelation for each channel     
    acmat = generate_multich_autocorrelation(small_audioseg,fs)
    
    b,a = signal.butter(1, )
    acc_lp_ch = signa
    channel = 1
    thresh = np.percentile(np.abs(acmat[:,channel]), 25)
    rel_thresh = thresh/np.max(np.abs(acmat[:,channel]))
    peaks_00, _ = signal.find_peaks(np.abs(acmat[:,channel]),height=rel_thresh,distance=(fs*0.5e-3))

    plt.figure()
    plt.plot(np.abs(acmat[:,channel]))
    plt.plot(peaks_00, np.abs(acmat[peaks_00,channel]),'*')    