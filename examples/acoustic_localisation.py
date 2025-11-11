# -*- coding: utf-8 -*-
"""
Acoustic localisation of the Sibbe Sept 2025
============================================

Created on Wed Oct  1 16:32:50 2025

@author: theja
"""
import soundfile as sf
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd
import sys 
import scipy
import scipy.signal as signal
sys.path.append('C:/Users/theja/Documents/research_repos/batracker')
sys.path.append('C:/Users/theja/Documents/software/gccestimating')
from gccestimating import GCC
from batracker.localisation import friedlander_1987 as f1987
from batracker.localisation import spiesberger_wahlberg_2002 as sw2002
from batracker.localisation import schau_robinson_1987 as sr87
 # from mini_ccg.timediffestim_copy import *
# from mini_ccg.graph_manip import make_consistent_fls, make_ccg_matrix
# from mini_ccg.combineall import combine_all

euclidean = scipy.spatial.distance.euclidean
import pandas as pd
np.random.seed(560085)

#spiesberger_wahlberg2002
friedlander = f1987.solve_friedlander1987
spieswahl_berg = sw2002.spiesberger_wahlberg_solution
schaurob = sr87.schau_robinson_solution
dB = lambda X: 20*np.log10(X)
#%%
totalstation_yxz = pd.read_csv('sibbesept2025.csv', header=None)
totalstation_xyz = totalstation_yxz.loc[:,[0,2,1,3,4]]
totalstation_xyz.columns=['pointname','x','y','z','pt-type']

ch_order_ttstation = ['ch'+str(chnum) for chnum in range(1,9)]

ttstation_ctrl = pd.concat([totalstation_xyz[totalstation_xyz.loc[:,'pointname']==ch] for ch in ch_order_ttstation]).reset_index(drop=True)
ttstation_ctrl.to_csv('sibbegroeve_sept2025_micxyz.csv')
ttstation_ctrl_xyz = ttstation_ctrl.loc[:,'x':'z'].to_numpy()

#%%
calib_file = 'audio/multichannel_2025-09-27_19-30-15.wav'
bat_file = 'audio/multichannel_2025-09-27_23-48-00.wav'
start, stop = 10.486, 10.496 # 2.235, 2.25 # 2.485, 2.504# 2.723, 2.748 
fs = sf.info(calib_file).samplerate
audio, fs = sf.read(calib_file, start=int(fs*start), stop=int(fs*stop))
b,a = signal.butter(2, np.array([10e3,60e3])/(fs*0.5), 'band')
fn_bandpass = lambda X: signal.filtfilt(b,a,X)
audio_bp = np.apply_along_axis(fn_bandpass, 0, audio)

silent_chunk, _ = sf.read(calib_file, start=int(fs*0.134), stop=int(fs*0.199))
silent_bp = np.apply_along_axis(fn_bandpass, 0, silent_chunk)
silent_peak = np.max(abs(silent_bp), 0)
silent_abs_dBpeak = dB(silent_peak)
#%%

durn = 10e-3
t = np.linspace(0,durn,int(fs*durn))
start_f, end_f = 45e3, 15e3
linear_sweep = signal.chirp(t, start_f, t[-1], end_f)
linear_sweep *=signal.windows.tukey(linear_sweep.size, 0.95)
sf.write('10ms_sweep_45-15_kHz_sibbegroevesept2025.wav', linear_sweep, fs)
#%%


good_chs = []
good_toas = []
for ch in range(audio_bp.shape[1]):
    channel_peak = np.max(abs(audio_bp[:,ch]))
    if dB(channel_peak) >= silent_abs_dBpeak[ch] + 6:
        cc_position = signal.correlate(audio_bp[:,ch],
                                       linear_sweep, mode='same', method='fft')
        cc_envelop = abs(signal.hilbert(cc_position))
        min_ht = np.percentile(cc_envelop, 95)
        smoothed_cc = signal.convolve(cc_envelop, np.ones(96)/96, mode='same')
        det_peaks, _ = signal.find_peaks(cc_envelop, width=int(fs*30e-6),
                                          height=min_ht)
        #det_peaks = np.array([np.argmax(smoothed_cc)])
        first_peak = det_peaks[0]
        good_chs.append(ch)
        good_toas.append(first_peak/fs)
        
        plt.figure()
        aa = plt.subplot(412)
        t_ch = np.linspace(0, cc_position.size/fs, cc_position.size)
        plt.plot(t_ch,smoothed_cc)
        plt.plot(t_ch, cc_envelop)
        plt.plot(det_peaks/fs, cc_envelop[det_peaks],'*')
        plt.hlines(min_ht, 0, cc_envelop.size/fs)
        plt.subplot(411, sharex=aa)
        t_ch = np.linspace(0, cc_position.size/fs, cc_position.size)
        plt.plot(t_ch,cc_position)
        
        detection = first_peak/fs
        plt.subplot(413, sharex=aa)
        plt.specgram(audio_bp[:,ch], Fs=fs, NFFT=96, noverlap=90)
        plt.vlines([detection, detection+5e-3], 0, 15e3)
        plt.ylim(200, 95e3)
        
        plt.subplot(414, sharex=aa)
        null_audio = np.zeros(audio_bp[:,ch].size)
        try:
            ind_start = first_peak-int(fs*5e-3)
            null_audio[ind_start:ind_start+linear_sweep.size] += linear_sweep
            plt.specgram(null_audio, Fs=fs, NFFT=96, noverlap=90)
        except:
            pass
            
        plt.savefig(f'channel_{ch}_crosscorr.png')
        plt.close()
#%%
tdos_refch1 = np.array(good_toas)
tdos_refch1 -= tdos_refch1[0]

good_mics_xyz = ttstation_ctrl.loc[good_chs, 'x':'z'].to_numpy()


rangediff = tdos_refch1[1:]*343

fr_out = friedlander(good_mics_xyz, rangediff, j=0)
print(fr_out)
#%%
