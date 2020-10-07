import unittest
import numpy as np
import scipy.signal as signal
from batracker.signal_detection.detection import *


def check_num_detections(detections, exp_num_detections):
    '''Internal test function to check the number of detections.
    '''
    obtained_detections = len(detections)
    if obtained_detections!=exp_num_detections:
        raise ValueError(f'Expected detections ({exp_num_detections}) dont match obtained ({obtained_detections})')

def check_duration(one_detection, expected_duration, margin=0.1):
    '''
    Parameters
    ----------
    one_detection : tuple with 2 entries
    expected_duration: float
    margin : float>0
        The tolerated +/- factor in percentage. ie. if obtained duration 
        is 0.0045 and the actual is 0.005 s. This corresponds to 0.9 accuracy. 
        If the margin in 0.1, then this passes as it is within 0.9-1.1 of the
        actual. Defaults to 0.1. 
    
    Raises
    ------
    ValueError
        If the accuracy of the duration estimate is beyond the acceptable margin of 
        error. 
    
    '''
    start, stop = one_detection
    obtained_duration = stop-start
    accuracy = obtained_duration/expected_duration
    
    lower, higher = 1-margin, 1+margin
    if np.logical_or(accuracy<lower, accuracy>higher):
        raise ValueError(f'Obtained duration accuracy is {accuracy}, and not within {lower}-{higher}')

def generate_5ms_bat_call():
    '''
    Returns
    -------
    bat_call
    fs
    duration
    
    '''
    duration = 0.005
    fs = 250000
    t = np.linspace(0, duration, int(fs*duration))
    bat_call = signal.chirp(t,90000, 25000, t[-1])
    bat_call *= 0.5
    return bat_call, fs, duration


class dBrmsDetection(unittest.TestCase):
    
    def setUp(self):
        
        bat_call, self.fs, self.duration = generate_5ms_bat_call()
        
        background = -60 # dB rms
        self.audio = np.random.normal(0, 10**(background/20), self.fs)
        
        self.sound_start = 0.05

    
    
    def test_simple(self):
        '''
        pure signal without any kind of windowing.
        Checks that the start and stop times are detected within 1 ms error.
        '''

        t = np.linspace(0, self.duration, int(self.fs*self.duration))
        bat_call = signal.chirp(t,90000, 25000, t[-1])
        bat_call *= 0.5
        sound_stop = self.sound_start+self.duration
        
        start, end = np.int32(np.array([self.sound_start,
                                        sound_stop])*self.fs)
        self.audio[start:end] += bat_call
        
        detections = dBrms_detector(self.audio, self.fs, dbrms_threshold = -12)
        check_num_detections(detections, 1)
        check_duration(detections[0],  self.duration)

class EnvelopeDetection(unittest.TestCase):
    
    def setUp(self):
        
        self.bat_call, self.fs, self.duration = generate_5ms_bat_call()        
        background = -60 # dB rms
        self.audio = np.random.normal(0, 10**(background/20), self.fs)

        
    
    def test_simple_envelope(self):
        ''' 
        A box window bat call which needs to be detected. 
        '''
        sound_start = 0.5
        sound_stop = sound_start + self.duration
        start, stop = np.int32(np.array([sound_start, sound_stop])*self.fs)
        
        self.audio[start:stop] += self.bat_call
        
        envelope_dets = envelope_detector(self.audio, self.fs, threshold_dbpeak=-20)
        print(envelope_dets)
        check_num_detections(envelope_dets, 1)
        check_duration(envelope_dets[0],  self.duration)


if __name__ == '__main__':
    unittest.main()
