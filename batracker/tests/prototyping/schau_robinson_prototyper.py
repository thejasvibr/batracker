#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prototyping module: Schau & Robinson 1987
=========================================
This module contains test examples to check if things are working and if they are 
accurate.
"""

import numpy as np 
import scipy.spatial as spatial

# this needs to be fixed at some point...
import sys 
sys.path.append('../../')
from localisation import schau_robinson_1987
# %% 
# Let's imagine a 4 channel array somewhere with a couple of sound sources around it.
# To check if my Schau & Robinson implementation is working I need a couple of things
# 1) range diffferences between the mics to sources
# 2) array geometry
 


# M matrix
raw_array_geom = np.random.rand(12).reshape(4,3) # each row is one mic
array_geom_ref_m4 = raw_array_geom-raw_array_geom[-1,:]

# %% Making the simulated data
source_position = np.array([10, 10, 10])
source_position_ref_m4 = source_position-raw_array_geom[-1,:]
M = source_position_ref_m4[:-1,:]

#%% Range to source
D_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom_ref_m4,
                               source_position_ref_m4)
# %% Range difference , ref. Range of mic 4
d_matrix = D_matrix - D_matrix[-1]
d = d_matrix[:-1]

#%% Mic to origin distances
R_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom_ref_m4,
                               np.zeros(3))

R = R_matrix[:-1]

# %% Delta matrix
delta = R**2-d**2 


# %% Solutions!!

# %% Equation 10
schau_robinson_1987.equation_10()

