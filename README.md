# batracker
This is a (not-so) stable package to perform acoustic localisation. Especially the ```localisation``` module is heavily in use by members of the [Active Sensing Collectives](activesensingcollectives.com) lab. 

### What batracker does 

- acoustic tracking using pre-recorded multichannel audio data.
and then output the start and end times of the bat calls in across all channels. The only constraints are the input and the output, the rest is upto you. 

### Localisation 
Localisation can be achieved through multiple algorithms and methods. One of the most standard localisation methods uses time-difference-of-arrivals (TDOA) of a single sound across channels. The TDOAs can be calculated through many methods (cross-correlation, GC-PHAT, etc). TDOAs

Acoustic tracking in 2D can be performed with two synchronised channels, and in 3D with >= 4 channels. ```batrack``` handles acoustic tracking with >= 2 channels and along with inbuilt tracking algorithms, is built to be extensible. The user can also define custom localisation algorithms for tracking. For instance, even among the methods using TDOAs, there are various formulations, the 

## Program Inputs and Outputs

### Inputs
The most basic information ```batracker``` requires are 1) a multichannel audio recording and 3) the xyz coordinates of the microphone array with which the recordings were made. If you only have recordings and no direct measurements of the microphone position, try out other packages like [StructureFromSound](linkhere) and come back with the microphone position estimates. 

### Outputs
The outputs produced by ```batracker``` are the 3D positions of the emitted sound. Especially when dealing with recordings from the field - these position estimates need to be taken with a pinch of salt! Remember that under low signal-to-noise ratios the position estimation will be off, and moreover, in the case of tracking animals, there are many variables that can't be controlled (eg. call directionality, reflections) which will invariably lead to poor fixes. 
