#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Correspondence matching prototyping
===================================
Implements the most basic type of correspondence matching -  which matches
signals across all channels to a single source emission. 

In this prototyping example, the simplest algorithm is used to find correspondences. 
The maximum inter-channel delay is calculated based on the longest dimension of the 
array, from the reference mic to the focal mic. 


Important variables
-------------------

#. Maximum distance from reference mic to other mics
#. Inter-sound interval
#. Sound duration 

Scenarios to handle
-------------------
The most important thing to check is if correspondence matching based on the 
max array distance works. It `will` work when 

What happens when the inter-sound interval drops below that of the max-array distance, 
and which rules work better then?

Possible approaches
"""

import numpy as np 
import scipy.distance as distance

# %%
# Example source emissions with mic array
# Make different sets with different inter-sound delays 

mic_array =  np.array([[0,0,1],
                       [0,1,0],
                       [1,0,0],
                       [0,0,0]])

sources = np.array([[10,0,0],
                    [0,10,10],
                    [0,10,10]
                    ])

intersound_interval = [0.004, 0.010, 0.02]

