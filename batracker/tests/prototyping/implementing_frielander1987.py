#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementing Friedlander 1987
=============================

@author: autumn
"""
import numpy as np 
# %% 
# Notation
# :math:`x_{i} = [x,y,z]^{T}` : the position of the ith sensor
# :math:`R_{is}` : distance between source and sensor i 
# :math:`R_{io}` : distance between origin and sensor i 
# :math:`r_{ij}= R_{is}-R_{js}`

# %%
# Equation 7b in action
# `j` is the reference sensor

def make_Sj(j, mic_posns):
    '''
    Parameters
    ----------
    j : int >=0
        The row number of the reference microphone
    mic_posns: np.array
        A Nmicx3 np.array
    
    Returns 
    -------
    Sj : np.array
        An Nchannels-1 x 3 array with :math:`\Delta` x,y,z descirbing the variable 
        :math:`S_{j}` in equation 7b. 
    '''
    num_channels, coods = mic_posns.shape
    if coods !=3:
        raise ValueError(f'Expected 3 coordinates, but got {coods} coordinates')
    
    Sj = np.zeros((num_channels-1, 3))
    
    valid_i_array = np.arange(num_channels)
    value_to_delete = int(np.argwhere(valid_i_array==j))
    valid_i_array = np.delete(valid_i_array, value_to_delete)
    
    for rownum, i in enumerate(valid_i_array):
        Sj[rownum,:] = mic_posns[i,:] - mic_posns[j,:]
    return Sj

def make_muj():
    '''

    
    Parameters
    ----------
    
    
    Returns
    -------
    
    Notes
    -----
    I think there may be a typo in eq. 7c because the first row entry should bee
    :math:`{R_{1o}}^2 - -{R_{jo}}^2 - {r_{1j}^2}` rather than just 
    :math:`{R_{1o}} - -{R_{jo}}^2 - {r_{1j}^2}`. This is because :math:`\mu_{j}`
    is obtained from eq. 6, where the term is :math:`({R_{io}^2}-{R_{jo}^2})-{R_{ij}^2}`
    '''
    
 
     
     

