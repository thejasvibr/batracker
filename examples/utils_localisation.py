# -*- coding: utf-8 -*-
"""
Utils to handle long audio files
================================
Created on Wed Oct  1 21:09:40 2025

@author: theja
"""
import numpy as np 
import sys
sys.path.append('C:/Users/theja/Documents/research_repos/batracker')
from batracker import localisation
from batracker.localisation import spiesberger_wahlberg_2002 as sw2002
from batracker.localisation import friedlander_1987 as f1987 
from batracker.localisation import schau_robinson_1987 as sr87
from batracker.tdoa_estimation import tdoa_estimators as tdoa_est
friedlander = f1987.solve_friedlander1987
wahlbergspiesberger = sw2002.spiesberger_wahlberg_solution
from scipy import spatial
def dB(X):
    return 20*np.log10(X)

def dB_to_linear(X):
    return 10**(X/20)

def break_into_chunks(multich_audio, num_chunks):
    '''
    '''
    return np.array_split(multich_audio, num_chunks)

    

def channels_above_threshold(multich_audio, peak_thresholds):
    '''
    '''
    peak_values = np.max(abs(multich_audio),0)
    above_threshold_channels = []
    for i,(peakval, threshold) in enumerate(zip(peak_values, peak_thresholds)):
        if peakval >= threshold:
            above_threshold_channels.append(i)
    return above_threshold_channels



def toa_from_known_signal(audio, knownsignal, **kwargs):
    '''
    Gives the TOA of the first peak - WARNING - many hard-coded parameters
    '''
    cc_position = signal.correlate(audio,
                                   knownsignal, mode='same', method='fft')
    cc_envelop = abs(signal.hilbert(cc_position))
    
    smoothed_cc = signal.convolve(cc_envelop, np.ones(96)/96, mode='same')
    min_ht = np.percentile(smoothed_cc, 99)
    
    det_peaks, _ = signal.find_peaks(smoothed_cc, width=int(fs*50e-6),
                                      height=min_ht)
    return det_peaks[0]

def tdoa_from_known_signal_multich(multich_audio, **kwargs):
    '''
    Parameters
    ----------
    multich_audio : (M,N) np.array 
        Multichannel audio with M samples, and N channels
    
    Keyword Argument
    ----------------
    knownsignal : np.array
        A known signal to cross-correlate with 
    
    Returns 
    -------
    tdoa : (N,) np.array 
        Time-difference-of-arrival in samples with regards to the 
        first channel. 
    '''
    all_toas = []
    for ch in range(multich_audio.shape[1]):
        toa = toa_from_known_signal(multich_audio[:,ch], kwargs['knownsignal'])
        all_toas.append(toa)
    tdoa = np.array(all_toas)
    tdoa -= tdoa[0]
    return tdoa

def localise_sounds_friedlander(audio, fs, mic_geom, **kwargs):
    '''
    Parameters
    ----------
    audio 
    fs
    mic_geom 
    
    Keyword Arguments
    -----------------
    tdoa_function : function 
        Function with the first argument of tdoa_function(audio, **kwargs)
        and all other parameters referred to with keyword arguments. 
        The function must output time-difference-of-arrivals in sample units, 
        and with reference to channel 1 always. 
    
    '''
    tdoa_function = kwargs.get('tdoa_function', None)
    tdoas = tdoa_function(audio, **kwargs)
    rangediff = (tdoas/fs)*343 # metres
    source = friedlander(mic_geom, rangediff[1:], j=0)
    return source
    

def check_loop_residual_from_toas(t1, t2, t3):
    '''
    '''
    t21 = t2 - t1 
    t31 = t3 - t1
    t32 = t3 - t2 
    return calc_loop_sum(t21, t31, t32)

def calc_loop_sum(t21, t31, t32):
     return t32 + t21 - t31

if __name__ == '__main__':
    import soundfile as sf
    import matplotlib.pyplot as plt
    import pandas as pd
    import scipy.signal as signal 
    #%%
        
    totalstation_yxz = pd.read_csv('sibbesept2025.csv', header=None)
    totalstation_xyz = totalstation_yxz.loc[:,[0,2,1,3,4]]
    totalstation_xyz.columns=['pointname','x','y','z','pt-type']
    
    ch_order_ttstation = ['ch'+str(chnum) for chnum in range(1,9)]
    
    ttstation_ctrl = pd.concat([totalstation_xyz[totalstation_xyz.loc[:,'pointname']==ch] for ch in ch_order_ttstation]).reset_index(drop=True)
    ttstation_ctrl_xyz = ttstation_ctrl.loc[:,'x':'z'].to_numpy()
        
    #%%
    
    calib_file = 'audio/multichannel_2025-09-27_19-30-15.wav'
    bat_file = 'audio/multichannel_2025-09-27_23-48-00.wav'
    input_file = bat_file
    

    start, stop = 0,10 #3.35, 3.7#8.9, 10.06#4, 8#2.235, 2.25 # 2.485, 2.504# 2.723, 2.748 
    fs = sf.info(input_file).samplerate
    minfreq, maxfreq = 20e3, 95e3
    audio, fs = sf.read(input_file, start=int(fs*start), stop=int(fs*stop))
    b,a = signal.butter(2, np.array([minfreq,maxfreq])/(fs*0.5), 'band')
    fn_bandpass = lambda X: signal.filtfilt(b,a,X)
    audio_bp = np.apply_along_axis(fn_bandpass, 0, audio)
    
    silent_chunk, _ = sf.read(calib_file, start=int(fs*1.55), stop=int(fs*1.7))
    silent_bp = np.apply_along_axis(fn_bandpass, 0, silent_chunk)
    silent_peak = np.max(abs(silent_bp), 0)
    silent_abs_dBpeak = dB(silent_peak)
    #%%
    
    
    
    durn = 10e-3
    t = np.linspace(0,durn,int(fs*durn))
    start_f, end_f = 45e3, 15e3
    linear_sweep = signal.chirp(t, start_f, t[-1], end_f)
    linear_sweep *=signal.windows.tukey(linear_sweep.size, 0.95)
    #%%
    audio_chunks = break_into_chunks(audio_bp, int((audio_bp.shape[0]/fs)/25e-3))
    all_sources = []
    t = 0
    all_t = []
    chunk_ids = []
    for i, chunk in enumerate(audio_chunks):
        t += chunk.shape[0]/fs
        ch_above = channels_above_threshold(chunk, dB_to_linear(silent_abs_dBpeak+10))
        if len(ch_above)>4:
            valid_mics = ttstation_ctrl_xyz[ch_above,:]
            valid_audio = chunk[:, ch_above]
            chunk_ids.append(i)
            all_t.append(t)
            try:
                tdoas = tdoa_from_known_signal_multich(valid_audio, knownsignal=linear_sweep)
                xyz = localise_sounds_friedlander(valid_audio, fs, mic_geom=valid_mics,
                                            tdoa_function=tdoa_from_known_signal_multich, 
                                            knownsignal=linear_sweep)
                print(t+start, f'chunk {i}', xyz, ch_above)
                
                all_sources.append(xyz.flatten())
            except:
                pass
    #%%
    consecutive_distances = []
    all_sources = np.array(all_sources)
    big_distmatrix = spatial.distance_matrix(all_sources, all_sources)
    for j in range(1,all_sources.shape[0]):
        consecutive_distances.append(big_distmatrix[j,j-1])
    #%%
    all_sources = []
    for i,chunkid in enumerate(chunk_ids):
        chunk = audio_chunks[chunkid]
        ch_above = channels_above_threshold(chunk, dB_to_linear(silent_abs_dBpeak+10))
        valid_mics = ttstation_ctrl_xyz[ch_above,:]
        valid_audio = chunk[:, ch_above]
        chwise_peaks = {}
        num_subplots = len(ch_above)
        f, axs = plt.subplots(num_subplots, 1, sharey=True, sharex=True)
        f1, axs1 = plt.subplots(num_subplots, 1, sharey=True, sharex=True)
        tdoas = []
        plt.sca(axs1[0])
        plt.plot(valid_audio[:,0])
        for ch in range(1,valid_audio.shape[1]):
            # cc_position = signal.correlate(valid_audio[:,ch],
            #                                 linear_sweep, mode='same', method='fft')
            norm_ch0 = valid_audio[:,0]/np.max(abs(valid_audio[:,0]))
            norm_ch = valid_audio[:,ch]/np.max(abs(valid_audio[:,ch]))
            # cc_position = signal.correlate(valid_audio[:,ch], valid_audio[:,0], mode='same', 
            #                                 method='fft')
            cc_position = signal.correlate(norm_ch, norm_ch0, mode='same', 
                                             method='fft')
            cc_envelop = abs(signal.hilbert(cc_position))
            smoothed_cc = signal.convolve(cc_envelop, np.ones(1)/1, mode='same')
            min_ht = np.percentile(smoothed_cc, 99)
            
            det_peaks, _ = signal.find_peaks(smoothed_cc,
                                              height=min_ht)
            taller_peak = det_peaks[np.argmax(smoothed_cc[det_peaks])]
            tdoa = taller_peak - cc_envelop.size*0.5
            tdoas.append(tdoa)
            # plt.sca(axs[ch])
            # plt.title(f'Global Channel#: {ch_above[ch]}, time:{all_t[i]+start}')
            # plt.plot(smoothed_cc)
            # plt.plot(cc_envelop)
            # plt.plot(det_peaks, smoothed_cc[det_peaks], '*')
            # plt.plot(taller_peak, smoothed_cc[taller_peak], 'r*')
            # plt.sca(axs1[ch])
            # plt.plot(valid_audio[:,ch])
            # plt.title(f'Global Channel#: {ch_above[ch]}, time:{all_t[i]+start}')
            
        source_df = pd.DataFrame(data={}, index=[0],
                                 columns=['t_chunk', 'chunk_id', 'x', 'y','z', 'mic_nums'])
        if len(tdoas)>3:
            tdoas = np.array(tdoas)
            rangediff = (tdoas/fs)*343
            source = friedlander(valid_mics, rangediff, j=0).flatten()
            source_df.loc[0,'chunk_id'] = chunkid
            source_df.loc[0,'x':'z'] = source
            source_df.loc[0,'mic_nums'] = str(ch_above)
            
            all_sources.append(source_df)
                    
    all_t = np.array(all_t) + start
    all_sources = pd.concat(all_sources)
    all_sources_xyz = all_sources.loc[:,'x':'z'].to_numpy(dtype=np.float64)
    all_sources.to_csv('bat_localisations.csv')
    
    #%%
    plt.figure()
    a0 = plt.subplot(111, projection='3d')
    plt.plot(ttstation_ctrl_xyz[:,0], ttstation_ctrl_xyz[:,1], ttstation_ctrl_xyz[:,2], '*')
    for i, each in enumerate(ttstation_ctrl_xyz):
        a0.text(*each, 'channel '+str(i+1))
    plt.plot(all_sources_xyz[:,0], all_sources_xyz[:,1], all_sources_xyz[:,2],'-*')
    for i, each in enumerate(all_sources_xyz):
        a0.text(*each, f't: {np.round(all_t[i]+start,3)}')
    a0.set_aspect('equal')
    #%%
    accam_points = totalstation_xyz.loc[18:,'x':'z'].to_numpy()
    plt.figure()
    a0 = plt.subplot(111, projection='3d')
    plt.plot(ttstation_ctrl_xyz[:,0], ttstation_ctrl_xyz[:,1], ttstation_ctrl_xyz[:,2], '*')
    for i, each in enumerate(ttstation_ctrl_xyz):
        a0.text(*each, 'channel '+str(i+1))
    plt.plot(accam_points[:4,0], accam_points[:4,1], accam_points[:4,2], 'r*')
    plt.plot(accam_points[4:,0], accam_points[4:,1], accam_points[4:,2], 'g*')
    idx = -3
    plt.plot(all_sources_xyz[idx,0], all_sources_xyz[idx,1], all_sources_xyz[idx,2],'-*')
    a0.set_aspect('equal')
    
    
    
    #%%
    tdoas = np.array(tdoas)
    rangediff = (tdoas/fs)*343
    print(friedlander(valid_mics, rangediff, j=0))
        
        
        
    
