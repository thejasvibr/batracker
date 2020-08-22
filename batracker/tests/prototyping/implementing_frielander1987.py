#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementing Friedlander 1987
=============================


@author: autumn
"""
import numpy as np 
import scipy.spatial.distance as distance
# %% 
# Notation
# :math:`x_{i} = [x,y,z]^{T}` : the position of the ith sensor
# :math:`R_{is}` : distance between source and sensor i 
# :math:`R_{io}` : distance between origin and sensor i 
# :math:`r_{ij}= R_{is}-R_{js}`

# %%
# Equation 7b in action
# `j` is the reference sensor

def solve_friedlander1987(array_geometry, d, **kwargs):
    '''
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    pass


def friedlander1987(mic_posns, j, rij):
    '''
    '''
    # make Sj
    Sj = make_Sj(ref_mic, mic_posns)    
    # make muj
    muj = make_muj(j,rij,  mic_posns).reshape(-1,1)    
    # make rhoj
    rhoj = rij.copy().reshape(-1,1)
    # make Pj
    Pj = make_Pj(Sj)

    # find solution from eq. 16
    frac_term = rhoj.T.dot(Pj).dot(muj)/rhoj.T.dot(Pj).dot(rhoj)
    rhs = muj - frac_term.dot(rhoj)
    
    xs = np.linalg.solve(Sj, rhs)
    return xs
    


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
    
    valid_i_array = make_valid_i_array(j, num_channels)
    
    for rownum, i in enumerate(valid_i_array):
        Sj[rownum,:] = mic_posns[i,:] - mic_posns[j,:]
    return Sj

def make_muj(j, rij, mic_posns):
    '''

    Parameters
    ----------
    j : int >=0
        The row number of the reference microphone
    rij : np.array
        Difference in range distances (:math:`R_{is}-R_{js}`)
    mic_posns: np.array
        A Nmicx3 np.array    

    Returns
    -------
    mu_j : np.array
        N-1 x 1 np.array with the :math:`mu_{j}`
    
    Notes
    -----
    I think there may be a typo in eq. 7c because the first row entry should bee
    :math:`{R_{1o}}^2 - -{R_{jo}}^2 - {r_{1j}^2}` rather than just 
    :math:`{R_{1o}} - -{R_{jo}}^2 - {r_{1j}^2}`. This is because :math:`\mu_{j}`
    is obtained from eq. 6, where the term is :math:`({R_{io}^2}-{R_{jo}^2})-{R_{ij}^2}`
    Here I"m going with my own interpretation of what is correct. 
    
    The positions given in mic_posns are assumed to be centred around the origin 
    (0,0,0).
    
    '''
    num_channels = mic_posns.shape[0]                           
    valid_i = make_valid_i_array(j, num_channels)
    
    R_sq_io = np.sum(mic_posns[valid_i,:]**2.0,1)
    R_sq_jo = np.sum(mic_posns[j,:]**2.0)
    r_sq_ij = rij**2
    
    mu_j = 0.5 *( R_sq_io - R_sq_jo - r_sq_ij)
    
    return mu_j
    
    
def make_Pj(Sj):
    '''
    Pj term in eq. 13
    '''
    sjtsj = Sj.T.dot(Sj)
    Pj = np.eye(Sj.shape[0]) - Sj.dot(np.linalg.inv(sjtsj)).dot(Sj)
    return Pj


def distance_to_point(x,y):
    return distance.euclidean(x,y)

    
def make_valid_i_array(j, num_channels):
    '''
    Parameters
    ----------
    j:  0>int
        Index of reference sensor
    num_channels: int
        Number of sensors
    
    Returns
    -------
    valid_i_array: np.array
        num_channels-1 array with all channel numbers apart from j.        
    '''
    valid_i_array = np.arange(num_channels)
    value_to_delete = int(np.argwhere(valid_i_array==j))
    valid_i_array = np.delete(valid_i_array, value_to_delete)
    return valid_i_array
    
if __name__=='__main__':
        
    import tacost
    from tacost import calculate_toa as ctoa

    
    mic_posns = np.array([[0,0,1],
                          [1,0,0],
                          [-1,0,0],
                          [0,1,0],
                          [1,1,0]
                          ])
    
    
    test_pos = np.array([5.5, 100, 15.2])
    toa = ctoa.calculate_mic_arrival_times(test_pos,
                                           array_geometry=mic_posns)
    j = 2
    vsound = 338.0 
    valid_i = make_valid_i_array(j, mic_posns.shape[0])
    rij = (toa[valid_i]-toa[j])*vsound
    
    Sj = make_Sj(j, mic_posns)    
    # make muj
    muj = make_muj(j,rij,  mic_posns).reshape(-1,1)
    # make rhoj
    rhoj = rij.copy().reshape(-1,1)

    # %% 
    from scipy.linalg import circulant


    Dj = np.linalg.inv(np.diag(rhoj.flatten()))
    
    array_to_shift = np.concatenate((np.zeros(rhoj.size-1),
                                     np.array([1])                                     
                                   )).flatten()
    
    Z = circulant(array_to_shift)
    Mj = (np.eye(rhoj.size)-Z).dot(Dj) # equation 8a

    
    MjSj = Mj.dot(Sj)
    Mjmuj = Mj.dot(muj)
    # %%
    #     
    xs,resid, _,_ = np.linalg.lstsq(MjSj, Mjmuj)
    
    left_portion = np.linalg.inv(Sj.T.dot(Mj.T.dot(Mj.dot(Sj))))
    right_portion = Sj.T.dot(Mj.T.dot(Mj.dot(muj)))
    xs_2 = left_portion.dot(right_portion)
    print( xs,test_pos, xs_2)
