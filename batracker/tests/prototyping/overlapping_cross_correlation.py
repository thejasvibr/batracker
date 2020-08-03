#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" What to do when there are two calls which arrive very close to one another?

Solution 1:
    
    #. Do the cross-correlation anyway - the multiple peaks might be picked up, and 
    the relative channel delays suggested by both the peaks could be tried out anyway
    #. Do the cross-correlation after 'intra-audio normalisation' - that is, make sure
    the RMS profile within the sound region to be cross correlated is uniform. This 
    prevents false positive peaks from forming during the cross correlation. 

Broad findings:
    
    #. Hilbert envelope based audio-normalisation could work, at least for the 
    simple case of overlapping calls (without reverb etc.)

"""
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize'] = 100000
import scipy.io.wavfile as wavfile
import scipy.signal as signal 
import numpy as np

def make_wav_and_spec(X, fs=192000):
    time = np.linspace(0, X.size/fs, X.size)
    plt.figure()
    a0  = plt.subplot(211)
    plt.specgram(X, Fs=fs)
    plt.subplot(212, sharex=a0)
    plt.plot(time, X)


fs = 192000
start_f, end_f =  90000, 30000
durn = 0.004
num_samples = int(fs*durn)
t = np.linspace(0,durn,num_samples)
call = signal.chirp(t, start_f, t[-1], end_f)
call *= signal.hanning(call.size)
call *= 0.5


make_wav_and_spec(call)

# %%
# Now make a channel with no overlaps, and a channel with two calls arriving very
# close to  each other. 

ch1 = np.zeros(1920)
ch1[100:100+num_samples] += call
ch2 = np.zeros(1920)
ch2[300:300+num_samples] += call
ch2[700:700+num_samples] += call


make_wav_and_spec(ch2)


cc = signal.correlate(ch1, ch2, 'same')
plt.figure()
plt.plot(cc)
plt.vlines(cc.size/2.0, 0, np.max(cc))



# %% 
# What happens if one call is louder than the other in the channel with overlaps

mixed_ampch2 = np.zeros(1920)
mixed_ampch2[300:300+num_samples] += call*0.3
mixed_ampch2[700:700+num_samples] += call

make_wav_and_spec(mixed_ampch2)

cc_mixedamp = signal.correlate(ch1, mixed_ampch2, 'same')
plt.figure()
plt.plot(cc_mixedamp)
plt.vlines(cc_mixedamp.size/2.0, 0, np.max(cc_mixedamp))

# %% 
# The peaks are of very different heights which means that the delays will of course 
# be interpreted badly. Now, let's equalise the waveform somehow. 

hilbert_tr = signal.hilbert(mixed_ampch2)
hilbert_envelope = np.abs(hilbert_tr)
plt.figure()
plt.plot(mixed_ampch2)
plt.plot(hilbert_envelope)

max_env = np.max(hilbert_envelope)
amp_factor = hilbert_envelope/max_env

uniform_amp = mixed_ampch2*(1/amp_factor)

plt.figure()
plt.plot(hilbert_envelope)
plt.plot(uniform_amp)


processed_cc = signal.correlate(ch1, uniform_amp, 'same')
plt.figure()
plt.plot(processed_cc)
plt.vlines(processed_cc.size/2.0, 0, np.max(processed_cc))


# %% 
# Let's compare the raw overlapping calls CC with the intra-audio normlised
# The peaks from the intra-audio normalised CC is definitely much larger than 
# just the raw CC - which shows that the intra-audio normalisation definitely helps!

plt.figure()
plt.subplot(211)
plt.plot(cc_mixedamp)
plt.vlines(cc_mixedamp.size/2.0, 0, np.max(cc_mixedamp))
plt.title('Raw mixed amplitude CC')
plt.xticks([])
plt.subplot(212)
plt.plot(processed_cc)
plt.vlines(processed_cc.size/2.0, 0, np.max(processed_cc))
plt.title('Amplitude adjusted CC')


# %% 
# One thing to pay attention to - the SNR before actually doing the envelope 
# based normalisation. If the SNR is already low - then it'd mean the noise would 
# be amplified -- and mess up the CC. 











