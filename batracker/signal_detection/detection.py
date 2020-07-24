'''
Detection API
=============

Deals with the actual detection of signals in multichannel audio files. 
There are two problems that need to solved while detecting a signal of interest.
    #. within-channel signal detection
    #. across-channel correspondence matching

Within-channel signal detection
-------------------------------
This task involves `locally` checking if there are any signals of interest in one channel at a time. The exact methods used for 
the within-channel can be set by the user, though the simplest is of course a basic threshold-type detector. Whenever the 
signal goes beyond a particular threshold, a signal is considered to be in that region.


Across-channel correspondence matching
--------------------------------------
After having detected the signals across all channels, they now need to be matched to each other. For sensible time-difference-of-arrivals, 
the 'correct' signal detections across channels need to be matched. In single-animal cases, where the animal emits discrete sounds with 
large enough inter-sound intervals, this can be achieved simple by assuming all sounds within :math:`\Delta\ T` of each other are from the
same source. However, with larger arrays, and multiple animals things can get tricky. Moreover, sometimes it may not be possible to decide the
correspondence solely from the audio data itself, and may require the assistance of another data channel, eg. video. 

The simplest case `batracker` deals with is the single animal case, where the signal correspondence problem is easily solved by checking for 
all signals within :math:`\pm\ \DeltaT` of the focal signal.

'''

import scipy.signal as signal 
import numpy as np 


def detect_signal(audio,fs, **kwargs):
    '''
    Parameters
    ----------
    audio  : np.array
        Msamples x Nchannels array with the multicchannel audio.
    fs : float>0
        sample rate
    signal_detector : function, optional
        The function used to detect a signal. Defaults to the click_detector function.

    Returns
    -------
    multichannel_candidates : List with sublists.
        The start and stop times where the signal is in each channel. 
        Each sublist contains the candidate signal regions for a channel.
        Each candidate signal region is a tuple with the start and end times
        of the signal detected in that channel.

    Note
    ----
    The `signal_times` list can look like so 
    [ [(0.1,0.3),(0.5,0.8),(0.9,1.2)],
      [],
      [(0.2,0.9), (.9,1.2), (1.23,1.26)]
    ]
    This implies the input audio was a 3 channel recording. 
    The first channel showed three detections, the second channel had zero detections, 
    and the third channel had three detections too. Note that all channels do not
    have to have the same number of detections, and they may be different for a 
    variety of reasons. 

    By default, the simplest signal detection algorithm is chosen, the 'click-detection'
    algorithm, where the sound is considered to start and stop as long as the 
    waveform is above a given threshold. 
    
    See Also
    --------
    click_detector
    '''
    signal_detector = kwargs.get('signal_detector', click_detector)
    
    # apply click detection on each channel
    n_channels = audio.shape[1]
    multichannel_candidates = []
    for each_channel in range(n_channels):
        channel_candidates = signal_detector(audio[:,each_channel], fs, **kwargs)
        multichannel_candidates.append(channel_candidates)
    return multichannel_candidates
        

def click_detector(one_channel, fs, **kwargs):
    '''
    Any region above a given threshold is considered  valid sound, and 
    its start and stop times are output.

    Parameters
    ----------
    one_channel : np.array
        Audio for a single channel
    fs : float>0
        Samplerate
    threshold: float, optional
        The threshold above which a sound is considered a signal 
    
    Returns
    -------
    candidate_regions : list with tuples
        Each tuple consists of two floats with the 
        start and stop time of the candidate signal region

    '''
    threshold = kwargs.get('threshold', -20) # dB rms / PEAKS -DECIDE
    
    # Generate a moving rms or some kind of other sound profile ? -- 
    # Perhaps also a Hilbert transform envelope? 
    
    
    # Use scipy.ndimage.find_objects to get all continuous stretches of 
    # signal above the threshold 
    
    # reformat the ndimage tuples and return them.
    
    
    






