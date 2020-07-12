# batracker
flexible acoustic tracking software - mainly built with bats and such in mind

*User Warning: This is a fictional spec sheet, I'm writing this out only to understand what all needs to be done to make a working prototype of a program like this!*

### What batracker does 

- acoustic tracking using pre-recorded multichannel audio data.
- workflows to handle relatively complex audio scenes with the possibility of handling dense audio recordings (with many sounds close together in time)


### The main parts of the batracker system

- Like any tracking system there are two parts to it:
    - Signal/call detection
    - Localisation

### Signal detection 
This part of the package allows the detection of sounds according to various algorithms already in place, and those that can be custom written too. 
For instance, a very common algorithm used to detect bat calls is the 'click-detection' method, where sounds are detected based on when the waveform 
crosses a certain threshold level. Of course such methods may not work, especially in dense audio recordings with many calls in them. Signal detection in 
dense recordings can be done by the use of neural-network based methods. ```batrack``` allows you to train a network, let it do the signal detection, and keep
the rest of the workflow intact.

The following signal detection algorithms are implemented
    - click-detection : the signal of interest is understood to be a short and easily identifiable sound (eg. bat call). Signals are detected whenever they exceed a 
threshold value eg. X dB peak amplitude

You can also add your own custom signal detection algorithms - all custom detection must take the raw multichannel audio as input and output candidate time regions
for each channel. For example, you could use an already published bat call detection network eg. [batdetective](linktotheoriginalpaper here), use to it detect bat calls,  and then output the start and end times of the bat calls in across all channels. The only constraints are the input and the output, the rest is upto you. 

### Localisation 
Localisation can be achieved through multiple algorithms and methods. One of the most standard localisation methods uses time-difference-of-arrivals (TDOA) of a single sound across channels. The TDOAs can be calculated through many methods (cross-correlation, GC-PHAT, etc). TDOAs

Acoustic tracking in 2D can be performed with two synchronised channels, and in 3D with >= 4 channels. ```batrack``` handles acoustic tracking with >= 2 channels and along with inbuilt tracking algorithms, is built to be extensible. The user can also define custom localisation algorithms for tracking. For instance, even among the methods using TDOAs, there are various formulations, the 

## Program Inputs and Outputs

### Inputs
The most basic information ```batracker``` requires is the geometry of the microphone array with which the recordings were made. If you only have recordings and no direct measurements of the microphone position, try out other packages like [StructureFromSound](linkhere) and come back with the microphone position estimates. 

### Outputs
The outputs produced by ```batracker``` are the 3D positions of the emitted sound. Especially when dealing with recordings from the field - these position estimates need to be taken with a pinch of salt! Remember that under low signal-to-noise ratios the position estimation will be off, and moreover, in the case of tracking animals, there are many variables that can't be controlled (eg. call directionality, reflections) which will invariably lead to poor fixes. The raw positions need to be further processed by a separate package such as [position_inspector](fictionalpagehere).

