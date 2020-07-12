# batracker
flexible acoustic tracking software - mainly built with bats and such in mind


### What batracker does 

- acoustic tracking using pre-recorded multichannel audio data.
- workflows to handle relatively complex audio scenes with the possibility of handling dense audio recordings (with many sounds close together in time)


### The main parts of the batracker system

- Like any tracking system there are three parts to it:
    - Signal detection
    - Localisation
    - Point inspection (*Is this really required?..how about making it a separate and optional sub-package*)

### Signal detection 
This part of the package allows the detection of sounds according to various algorithms already in place, and those that can be custom written too. 
For instance, a very common algorithm used to detect bat calls is the 'click-detection' method, where sounds are detected based on when the waveform 
crosses a certain threshold level. Of course such methods may not work, especially in dense audio recordings with many calls in them. Signal detection in 
dense recordings can be done by the use of neural-network based methods. ```batrack``` allows you to train a network, let it do the signal detection, and keep
the rest of the workflow intact. 

### Localisation 
Acoustic tracking in 2D can be performed with two synchronised channels, and in 3D with >= 4 channels. ```batrack``` handles acoustic tracking with >= 2 channels and along with inbuilt tracking algorithms, is built to be extensible. The user can also define custom localisation algorithms for tracking.

### Point inspection (*TBD*)
Despite having the best automated tracking
