#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Understanding some of the things in the Scheuing &  Yang 2008 paper. 

Created on Mon Jul  5 10:48:01 2021

@author: thejasvi
"""
#%%
import numpy as np 
import scipy.signal  as signal 
import scipy.spatial as spatial 
import matplotlib.pyplot as plt 

#%%
# microphone array geometry
R = 1.2
theta = np.pi/3
tristar = np.row_stack(([0,0,0],
                        [-R*np.sin(theta), 0, -R*np.cos(theta)],
                        [R*np.sin(theta), 0, -R*np.cos(theta)],
                        [0,0, R]))

sound_pos = np.array([3,2,1])
reflection_source = np.array([4,2,0])
direct_indirect_sources = np.row_stack((sound_pos, reflection_source))
# direct path propagation:

dist_mat = spatial.distance_matrix(direct_indirect_sources, tristar)
#add the distance of propagation from source to reflection point
source_to_reflectionpoint = spatial.distance.euclidean(sound_pos, reflection_source)
dist_mat[1,:] += source_to_reflectionpoint

# make the direct

chirp_durn = 0.003
fs = 192000
t = np.linspace(0,chirp_durn,int(fs*chirp_durn))
chirp = signal.chirp(t,80000,t[-1],25000)
chirp *= signal.hann(chirp.size)*0.5


vsound = 340.0
audio = np.zeros((int(fs*0.03),4))
toa_sounds = dist_mat/vsound
toa_samples = np.int64(toa_sounds*fs)
for channel in range(4):
    random_atten = np.random.choice(np.linspace(0.2,0.9,20),2)
    start_direct, start_indirect = toa_samples[0,channel], toa_samples[1,channel]
    audio[start_direct:start_direct+chirp.size,channel] += chirp*random_atten[0]
    audio[start_indirect:start_indirect+chirp.size,channel] += chirp*random_atten[1]
audio += np.random.normal(0,1e-3,audio.size).reshape(audio.shape)

#random_toa = int(fs*0.019)
#audio[random_toa:random_toa+chirp.size,0] += chirp*0.1
#audio[random_toa+20:random_toa+20+chirp.size,0] += chirp*0.05

#%% 
# Now let's look at the simple case where there's 1 direct and 1 indirect path.s
plt.figure()
a0 = plt.subplot(411)
plt.specgram(audio[:,0],Fs=fs)
for i in range(1,4):
    plt.subplot(411+i, sharex=a0, sharey=a0)
    plt.specgram(audio[:,i],Fs=fs)

#%% 
# Generate the auto-corr for each channel 
multich_acc = np.apply_along_axis(lambda X: signal.correlate(X,X,'same'),0, audio)

plt.figure()
plt.plot(multich_acc[:,1])

#%% Generate the cross-corr for each channel pair
multich_cc = {}
for i in range(4):
    for j in range(1,4):
        if i!=j:
            multich_cc[i,j] = signal.correlate(audio[:,j],audio[:,i],'same')
#%%
plt.figure()
cc0 = plt.subplot(311)
plt.plot(multich_cc[(0,1)])
plt.subplot(312, sharex=cc0)
plt.plot(multich_acc[:,0])
plt.subplot(313, sharex=cc0)
plt.plot(multich_acc[:,1])

#%% extract auto-cor and cross-cor peaks 
cc_and_acc_peaks = lambda X:  signal.find_peaks(X, 0.11,distance=int(fs*1e-4))[0] 

acc_peaks = []
for each in range(4):
    raw_peaks = cc_and_acc_peaks(multich_acc[:,each])
    raw_peaks_re_centre = raw_peaks-multich_acc[:,each].size/2.0
    acc_peaks.append(raw_peaks_re_centre)

#%%
# now get the cc peaks 
cc_peaks = {}
for i in range(4) :
    for j in range(1,4):
        if i!=j:
            cc_peaks[(i,j)] = cc_and_acc_peaks(multich_cc[(i,j)])
##%% 
plt.figure()
plt.plot(multich_cc[(0,1)])
peakinds = cc_peaks[(0,1)]
plt.plot(peakinds, multich_cc[(0,1)][peakinds],'*')
#plt.subplot(312, sharex=cc1)
#plt.plot(multich_acc[:,0])
#plt.plot(acc_peaks[0], multich_acc[acc_peaks[0],0],'*')
#plt.subplot(313, sharex=cc1)
#plt.plot(multich_acc[:,1])
#plt.plot(acc_peaks[1], multich_acc[acc_peaks[1],1],'*')


#%% Now, perform the direct/echo path matching. Which of the cross-cor peaks
# can be explained by the presence of indirect path delays? 

cc_peak01 = cc_peaks[(0,1)]
peak_quality = multich_cc[(0,1)][cc_peak01 ] # initial values
gamma_twrm = 20
direct_path_hits = np.zeros(len(cc_peak01 ))

#%%
# also make the quality vector now
q_array = multich_cc[(0,1)][cc_peak01]

#%%

def calculate_TDE_qualityfactor(crosscor_ba, acc_a, acc_b):
    '''
    Implements direct path matching and calculates quality factor. 
    The relevant portions are described in eqns. 11-15
    
    Parameters
    ----------
    crosscor_ba : tuple
        Tuple containing the crosscorrelation and the peaks: (crosscorrelation, peaks).
        The crosscorrelation is of channel b with ref. to channel a.
    acc_a, acc_b : tuple
        Tuple with autocorrelation and peaks: (autocorreation, peaks_autocorr)
    twrm: int>0, optional
        The 'tolerance width of raster match' in samples. Defaults to 10 samples. 
    
    Returns 
    -------
    quality_factor : np.array
        Of same size as crosscorrelation peaks.
    
    Notes
    -----
    * All peak values should in units of sample indices - and not in time
    * The peaks for autocorrelation must be 'centred' (those that are on the left
    of the centre must have -ve values and on the right have +ve values)

    '''
    cc_ba, cc_peaks = crosscor_ba
    autocc_a, accpeaks_a = acc_a
    autocc_b, accpeaks_b = acc_b

    quality_factor = cc_ba[cc_peaks] #set initial values to r_kls
    
    for i, each_ccpeak in enumerate(cc_peaks[:-1]):
        # check the left-to-right matching for chA
        for other_ccpeak in enumerate(cc_peaks[i+1:]):
            ccpeak_gap = other_ccpeak - each_ccpeak
            if ccpeak_gap
        
        
        # check the right-to-left matching for chB
    
        #store matching ACC peaks into P_kkprime and P_llprime
        
        # 
def TFRM(autocorr_delay, cc_peak1, cc_peak2, tfrm):
    difference = np.abs(autocorr_delay) - np.abs(cc_peak2-cc_peak1)
    
    if difference < 0.5*tfrm:
        return 1 - np.abs(difference)/(0.5*tfrm)
    else:
        return 0

for i, each_ccpeak in enumerate(cc_peak01[:-1]):
    # search from L--R for ch0
    for j in range(i+1,len(cc_peak01)):
        peak_diff = cc_peak01[i] - cc_peak01[j]
        within_twrm = np.abs(np.abs(peak_diff)-acc_peaks[0])<= 0.5*gamma_twrm
        if np.sum(within_twrm)>0:
            direct_path_hits[i] += 1            

n_ccs = len(cc_peak01)-1
for i, each_ccpeak in enumerate(cc_peak01[:-1]):
    # search from L--R for ch0
    for j in range(i+1,len(cc_peak01)):
        peak_diff = cc_peak01[n_ccs-i] - cc_peak01[n_ccs-j]
        within_twrm = np.abs(np.abs(peak_diff)-acc_peaks[1])<= 0.5*gamma_twrm
        if np.sum(within_twrm)>0:
            direct_path_hits[n_ccs-i] += 1        
