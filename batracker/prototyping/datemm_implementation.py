#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DATEMM implementation attempt
=============================
Scheuing and Yang 2008 [1] present the DATEMM algorithm which implements many of
the ideas I was only vaguely intuiting. Most importantly the:
    * The auto-correlation of each channel tells us the 'lags' at which indirect paths of the signal arrived
        eg. Signal A at channel 0, has lags of 10, 25, 41, and channel 1 has lags of 25, 32, 5
    * The indirect path delays show up in combinations in the channel pair cross-correlations:
        eg. the cc has peaks at 15 (result of ch1:lag25 - ch0:lag10) or 16 (result of ch0:lag41-ch1:lag25)
    * Figuring out which of the channel-pair time delays can be explained by 

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

def detect_peaks(acc, fs, **kwargs):
    '''
    Peak detection is done by taking the input auto/cross-correlation, and
    performing peak detection. 



    Parameters
    ----------
    acc : np.array
        Auto/crosss-correlation of one audio channel 
    fs: float>0
        Sampling rate in Hz
    peak_threshold: 0<float<100, optional
        The percentile threshold of the autocorrelation signal 
        that sets the absolute value threshold for a peak. Defaults to 10..
    min_interpeak_distance: float>0, optional 
        expected interpeak distance., defaults to 0.2 milliseconds

    Returns 
    -------
    peaks : np.array
        With indices of peak locations. 
    
    See Also
    --------
    scipy.signal.find_peaks
    
    '''
    interpeak = int(fs*kwargs.get('min_interpeak_distance', 5e-4))
    thresh = np.percentile(acc, kwargs.get('peak_threshold',10))
    rel_thresh = thresh/np.max(acc)
    peaks, _ = signal.find_peaks(acc,
                                    height=rel_thresh,distance=interpeak)
    return peaks

#%%

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
    mch_gcc_peaks = {}
    
    for channel_pair, gcc in mch_gcc.items():
        mch_gcc_peaks[channel_pair] = detect_peaks(gcc,fs,peak_threshold=15)

    plt.figure()
    plt.plot(mch_gcc[0,1])
    plt.plot(mch_gcc_peaks[0,1], mch_gcc[0,1][mch_gcc_peaks[0,1]], '*')
    
    # choose top X peaks of time delay     
    mch_tdes = generate_multich_tdes(mch_gcc, fs, n_peaks=10) 
    #%%
    # Autocorrelation for each channel     
    acmat = generate_multich_autocorrelation(small_audioseg,fs, )
    multich_accpeaks = {}
    for each in range(4):
        multich_accpeaks[each] = detect_peaks(acmat[:,each],fs,
                                                    speak_threshold=15)

#%% Now check which of the cross-correlation peaks can be explained by 
# the result of a non-direct path. 

