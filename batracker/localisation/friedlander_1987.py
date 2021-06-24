#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementing Friedlander 1987
=============================
Friedlander 1987 [1] formulates the localisation problem as the source being the  
intersection point of multiple major axis lines. Each set of three mics
lies on an ellipsoid, where the major axis passes  through the source. 
The final source position is calculated as the point where all the major axes
of the ellipsoids intersect. See Figure 2. of the paper for a better idea. 



References
----------
[1] Friedlander, B. (1987). A passive localization algorithm and its accuracy analysis.
    IEEE Journal of Oceanic engineering, 12(1), 234-245.

TODO:
    
    #. Run tests with simulated data and check accuracy

"""
import numpy as np 
import scipy.spatial.distance as distance
from scipy.linalg import circulant


def solve_friedlander1987(array_geometry, d, **kwargs):
    '''
    Parameters
    ----------
    array_geometry: np.array
        A Nx3 array with xyz coordinates of N mics.
        The reference microphone channel number (j) used in the range difference
        estimation must be explicitly stated.
    d:  np.array
        A N-1x1 np.array with the range_differences to the source. 
        All range_differences (:math:`r_{ij}`), are calculated by
        taking :math:`R_{is}-R_{js}`, where :math:`R_{is}` is the direct
        range from mic i to the source, and :math:`j` is the reference
        microphone.
    j : int >=0
        The channel number of the reference microphone.
        The first channel is the 0th.
    use_analytical : bool, optional
        Whether to use the analytical solution to find the source position
        or not. Defaults to True. If False, a least squares routine is 
        used to find the source position. The difference is likely to be 
        negligible for most cases.
        
    Returns
    -------
    x_real_world : list
        If the array_geometry consists of >4 mics, then the output  is 
        a list with a 3x1 np.array with the x,y,z positions of the source.
        If the array_geometry consists of 4 mics, then the solution is 
        a list with 2 candidate xyz source positions. 
    '''
    
    j = kwargs['j']
    
    num_channels, coods = array_geometry.shape
    if num_channels <= 4:
        raise IndexError(f'Friedlander 1987 supports only >4 mics for a 3d localisation.\
                         The input audio has {num_channels} mics')

    
    if coods != 3:
        raise IndexError(f'The current array geometry is {coods} dimensional.Only 3d coordinates accepted')
    
    x_real_world = friedlander1987(array_geometry, j, d, 
                                           kwargs.get('use_analytical', True))
    
    return x_real_world
    
    

def friedlander1987(mic_posns, j, rij, use_analytical=True):
    '''
    
    
    '''
    Sj = make_Sj(j, mic_posns)    
    muj = make_muj(j,rij,  mic_posns).reshape(-1,1)
    rhoj = rij.copy().reshape(-1,1)
    Dj = np.linalg.inv(np.diag(rhoj.flatten()))
    
    array_to_shift = np.concatenate((np.zeros(rhoj.size-1),
                                     np.array([1])                                     
                                   )).flatten()
    Z = circulant(array_to_shift)
    Mj = (np.eye(rhoj.size)-Z).dot(Dj) # equation 8a
   
    MjSj = Mj.dot(Sj)
    Mjmuj = Mj.dot(muj)
    if use_analytical:
        left_portion = np.linalg.inv(Sj.T.dot(Mj.T.dot(Mj.dot(Sj))))
        right_portion = Sj.T.dot(Mj.T.dot(Mj.dot(muj)))
        xs = left_portion.dot(right_portion)
    else:
        xs,resid, _,_ = np.linalg.lstsq(MjSj, Mjmuj)
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
    

