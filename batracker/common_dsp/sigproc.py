#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common Signal Processing Functions 
==================================
This module is a collection of all the common signal processing related 
functions like dB, rms, etc. 
"""

import numpy as np 

dB = lambda X : 20*np.log10(np.abs(X))

def rms(X):
    '''Root mean square of a signal '''
    return np.sqrt(np.mean(X**2.0))
