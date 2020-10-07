Signal detection
================
In this package, signal detection refers to the general act of 1)'recognising' the signal of interest within a given channel. 

Expected Inputs
----------------
The input data to run the different signal detection algorithms are the raw multichannel audio data. The most important requirement is 
that all the channels are `synchronised`!

Outputs to expect
-----------------
Each algorithm will at least provide the detected signals list, and may optionally 
produce additional information too (eg. confidence intervals, prediction probabilities., etc.).
All detected signals are output in a list with sublists holding candidate signal regions for each channel. The sublists are ordered in the order 
of the audio files. For example:

.. code-block:: python

    [ [(0.1,0.3),   (0.5,0.8),   (0.9,1.2)],
      [(0.15,0.35), (0.23, 0.45)          ],
      [(0.2,0.9),   (.9,1.2),  (1.23,1.26)],
      [(0.1,0.7),   (1.,1.1),  (1.29,1.38)]  
    ]

This corresponds to an input of a 4 channel audio. There are 3,2,3,3 candidate signal regions detected in 1st-4th channels respectively.


Each candidate signal region consists of a tuple with the start and end time in the audio channel. For example:

.. code-block:: python 

    # let's take the same signal detections above, but only look at the 1st channel:
    [ [(0.1,0.3),   (0.5,0.8),   (0.9,1.2)],
      ....... 
    ]

The candidate region :code:`(0.1,0.3)`, has a start time of 0.1s and and end time of 0.3s. This means the signal of interest 
is located between these two time stamps in the first channel. This means there are three candidate regions in the first channel
starting at 0.1, 0.5 and 0.9 seconds. 


The simplest method: threshold based detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Whenever the audio goes beyond a particular level (in RMS), a signal is considered detected.
This is a common (and relatively robust) algorithm that's used in many programs. 

- One thing to consider with the threshold method is that if the signal isn't received with the same intensity on all channels (and also in case of mic directionality  etc.) -- then the 'wrong parts' of a signal will be cross-correlated - leading to poor localisations. How to overcome this problem? 


Less simple methods
~~~~~~~~~~~~~~~~~~~
Threshold based methods could lead to false positives, where there are similar non-target sounds in the same 
frequency range (eg. two species calling at the same time). Here there are a bunch of options



Signal Detection API
--------------------

.. automodule:: batracker.signal_detection.detection
   :members:



