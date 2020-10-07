#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Measuring TDOA
==============


"""

import numpy as np 
import scipy.signal as signal 

def measure_tdoa(audio, fs,**kwargs):
    '''
    Parameters
    ----------
    audio : np.array
        Nchannels x Msamples
    fs : float>0
    vsound: float, optional
        Speed of sound in m/s, defaults  to 338m/s
    tdoa_method: function, optional 
        The function used to measure the TDOA. Defaults
        to cross_cor_w_ref_channel.  See Notes

    Returns 
    -------
    tdoa: np.array
        The TDOA between channels in seconds. 

    Notes
    -----
    The TDOAs between channels can be measured in various ways ie. cross-corre-
    lation, generalised cross-correlation, or a user-defined custom algorithm. 
    
    The tdoa_method function in general  should give out an array of inter-channel 
    TDOAs in units of samples.
    '''
    tdoa_method = kwargs.get('tdoa_method', tdoa_by_cc_w_ref_channel)
    tdoa_estimate_samples = tdoa_method(audio, **kwargs)
    tdoa = tdoa_estimate_samples/fs
    return tdoa

def tdoa_by_cc_w_ref_channel(audio, **kwargs):
    '''
    Cross correlates a common set of samples with a common reference channel.
    The distance of the cross-correlation peak from the centre is considered to
    the delay between each channel pair considered.
    
    Parameters
    ----------
    audio : np.array 
        Multichannel audio  
    ref_channel : int < Nchannels, optional 
        The reference channel used for all comparisons. Defaults to 
        the last channel. 
    
    Returns
    -------
    cc_tdoa : np.array
        All the estimated TDOAs in units of samples. The order of the 
        outputs is as follows. If :math:`N_{ref}` is the  reference
        channel,  where :math:`N_{ref}` can be any one of :math:`1-N_{channels}`, 
        the order of  the TDOAs are :math:`\Delta \tau_{ref-i}, where i=1:N_{channels}
        and i \neq N_{ref}`
    
    '''
    num_channels = audio.shape[1]
    ref_channel = kwargs.get('ref_channel', num_channels-1)
    

    all_delays = []
    for i in range(num_channels):
        if i != ref_channel:
            cc = signal.correlate(audio[:,i], audio[:,ref_channel], 'same')
            delay = np.argmax(cc)-(cc.size*0.5)
            all_delays.append(delay)
            
    cc_tdoa = np.array(all_delays)
    return cc_tdoa