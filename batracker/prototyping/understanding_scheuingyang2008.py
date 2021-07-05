#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Understanding some of the things in the Scheuing &  Yang 2008 paper. 

Created on Mon Jul  5 10:48:01 2021

@author: thejasvi
"""
import numpy as np 
import scipy.spatial as spatial 
import matplotlib.pyplot as plt 

# microphone array geometry
R = 1.2
theta = np.pi/3
tristar = np.row_stack(([0,0,0],
                        [-R*np.sin(theta), 0, -R*np.cos(theta)],
                        [R*np.sin(theta), 0, -R*np.cos(theta)],
                        [0,0, R]))

sound_pos = np.array([3,2,1])

reflection_source = np.array([4,2,1])

