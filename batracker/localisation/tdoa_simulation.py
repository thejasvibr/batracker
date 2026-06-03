# -*- coding: utf-8 -*-
"""
TDOA prediction
===============
Created on Wed Jun  3 11:25:41 2026

@author: theja
"""
import scipy.spatial as spl
import numpy as np 
from itertools import combinations 

def mic_pair_combis(nchannels):
    '''
    '''
    return list(combinations(range(nchannels), 2))

def generate_tdoa_predictions(source, mic_geom, vsound=340):
    '''
    '''
    mic_pairs = mic_pair_combis(mic_geom.shape[0])
    source_to_array_distances = spl.distance_matrix(source.reshape(-1,3), mic_geom).flatten()
    time_of_arrivals = source_to_array_distances/vsound # seconds - assuming t_emission = 0
    
    tdoa_predictions = np.zeros(len(mic_pairs))
    for ii,(i,j) in enumerate(mic_pairs):
        tdoa_predictions[ii] = time_of_arrivals[i] - time_of_arrivals[j]
    return tdoa_predictions.reshape(-1,1)


def generate_tdoa_predictions_multisource(sources, mic_geom):
    '''
    Parameters
    ----------
    sources : (M,3) np.array
        XYZ coordinates of M sources 
    mic_geom : (N,3) np.array
        XYZ coordinates of N mics
    
    Returns
    -------
    tdoa_sources : (M, (NxN-1/2)) np.array
        Expected TDOAs of all the mic pairs (columns) for each source (row-wise)
   
    '''
    nsources = sources.shape[0]
    mic_pair_inds = len(mic_pair_combis(mic_geom.shape[0]))
    tdoa_sources = np.zeros((nsources, mic_pair_inds))
    for i in range(nsources):
        tdoa_sources[i,:] = generate_tdoa_predictions(sources[i,:], mic_geom).flatten()
    return tdoa_sources
