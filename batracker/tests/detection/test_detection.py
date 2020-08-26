import unittest
import scipy.signal as signal
import batracker.signal_detection
from batracker.signal_detection import *

def check_num_detections(detections, exp_num_detections):
    '''
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
    


class dBrmsDetection(unittest.TestCase):
    
    def setUp(self):
        self.fs = 250000
        background = -60 # dB rms
        self.audio = np.random.normal(0, 10**(background/20), fs)
        
        self.sound_start = 0.05

    
    
    def test_simple(self):
        '''
        pure signal without any kind of windowing.
        Checks that the start and stop times are detected within 1 ms error.
        '''
        duration = 0.005
        t = np.linspace(0, duration, int(self.fs*duration))
        bat_call = signal.chirp(t,90000, 25000, t[-1])
        bat_call *= 0.5
        sound_stop = self.sound_start+duration
        
        start, end = np.int32(np.array([self.sound_start,
                                        sound_stop])*self.fs)
        self.audio[start:end] += bat_call
        
        detections = dBrms_detector(self.audio, self.fs, dbrms_threshold = -12)
        print(detections)
        check_num_detections(detections, 1)
        check_duration(detections[0],  duration)






if __name__ == '__main__':
    unittest.main()