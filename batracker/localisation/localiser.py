'''Localiser module

Deals with assigning signals across channels to the same sound, 
and localising the sound to its source
'''

import numpy as np 
#from schau_robinson_1987 import schau_robinson_solution


def localise(multichannel_candidates, microphone_positions, **kwargs):
    '''
    Parameters
    ----------
    multichannel_candidates : list with sublists.
        Each sublist contains candidate regions where the signal of interest
        has been detected
    microphone_positions : np.array
        Mmics x 3 np.array with the xyz positions of the mics. 
        Each mic is in a new row
    v_sound : float, optional
        Speed of sound in m/s. Defaults to 331 m/s. 
    formulation : str
        The mathematical formulation used to localise sounds
        Available options include 
        #. Schau & Robinson 1987, IEEE Trans. Acoust. Speech & Sig. Proc.
            

    Returns
    -------
    source_positions : np.array

    
    Notes
    -----
    Error in the input data can contribute to errors in the actual source localisation.
    
    Error in microphone position estimation can contribute to errors in source position
    estimation, especially when dealing with sources that are a bit further away [1]. 
    At least in air, variation in the speed of sound is rather unlikely to contribute to big 
    errors, and the error will also be homogenous in space (cause a uniform contraction
    or expansion of all position estimates).
    
    References
    ----------
    [1] Wahlberg, M., Møhl, B., & Teglberg Madsen, P. (2001).
        Estimating source position accuracy of a large-aperture hydrophone array 
        for bioacoustics. The Journal of the Acoustical Society of America, 109(1),
        397-406.
    
    See Also
    --------
    signal_detection.detection
    '''

    # generate correspondence map for all sounds 
    correspondence_map = match_sounds_to_each_other(multichannel_candidates,
                                                    microphone_positions,
                                                    **kwargs)

    # generate 3D positions for the correspondence mapped sounds
    source_positions = localise_source(correspondence_map, 
                                               multichannel_candidates,
                                                   microphone_positions, 
                                                               **kwargs)
    return source_positions


def match_sounds_to_each_other(multichannel_candidates,
                               microphone_positions,
                                                       **kwargs):
    '''
    Multiple sounds are detected within a channel and across channel. 
    Each detecting in a channel needs to be grouped together with 
    other detections across channels to proceed further with localisations. 
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    '''
    
    
def localise_source(correspondence_map, multichannel_candidates,
                                                       microphone_positions, 
                                                       **kwargs):
    '''
    '''
    
    
    
    
    
    
    
    
    
    
