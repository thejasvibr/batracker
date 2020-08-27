#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inbuilt Signal Detection Routines
=================================

`batracker` has two inbuilt signal detection routines:
    
#. Envelope detector
    
    The envelope detector generates the Hilbert envelope of the input audio. 
    Any region that is :math:`\geq` a threshold in dB peak is considered a 
    valid signal.

#. dB rms profile detector
    
    The moving dBrms profile of the input audio is generated.  Regions :math:`\geq`
    a threshold are considered valid signal.

This example will go on to showcase how signal detection works using these two
routines.

"""



