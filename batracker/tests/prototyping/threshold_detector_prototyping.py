#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Threshold detector protoyping
=============================

TODO
----

#. Write the actual detector function :p

"""

import os
import matplotlib.pyplot as plt
plt.rcParams['agg.path.chunksize']=10000
import scipy.io.wavfile as wav

# % Load the raw audio 

audio_path = os.path.join('simulated_audio','batracker_simple.wav')
fs, audio = wav.read(audio_path)

## % threshold detector
plt.figure()
plt.plot(audio[:,0])

## % Detect audio within one channel 
