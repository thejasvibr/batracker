#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Let's do NRM-DATEMM first. (Zannini et al. 2010) (non-raster matching)

@author: Thejasvi
"""


#%%
import copy
import numpy as np 
np.random.seed(82319)
import pandas as pd
import scipy.signal  as signal 
import scipy.spatial as spatial 
import matplotlib.pyplot as plt 
from itertools import combinations

from understanding_scheuingyang2008 import simulate_sound_propagation, filter_cc_peaks_by_plausibility_and_minrkl
from understanding_scheuingyang2008 import generate_multich_crosscorr, gamma_tftm
from understanding_scheuingyang2008 import generate_multich_autocorr
from batracker.localisation import spiesberger_wahlberg_2002 as sw02

def cc_and_acc_peaks(X, fs=192000, **kwargs):
    return signal.find_peaks(X, 0.11,distance=int(fs*1e-4))[0] 

def get_tdoa_matrix(source_estimate, array):
    '''
    '''
    nmics = array.shape[0]
    tdoas = np.zeros((nmics, nmics))
    
    R_source_mics = [spatial.distance.euclidean(array[each,:],source_estimate) for each in range(nmics)]
    for i in range(nmics):
        for j in range(nmics):
            tdoas[j,i] = R_source_mics[j] - R_source_mics[i]
    return tdoas
#%% 



# make simulated audio and generate cross-cor + auto-cor
vsound = 338.0

sim_audio, dist_mat, array_geom, (source, reflector) = simulate_sound_propagation()
tdoas = np.row_stack((dist_mat[0,1:]-dist_mat[0,0],dist_mat[1,1:]-dist_mat[1,0]))/vsound
fs = 192000
multich_cc =  generate_multich_crosscorr(sim_audio)
multich_acc = generate_multich_autocorr(sim_audio)
twrm_samples = 10
#%% extract peaks for each cc pair

crosscor_peaks = {}
for each_pair, crosscor in multich_cc.items():
    ch_b, ch_a = each_pair
    cc_peaks = cc_and_acc_peaks(crosscor)
    autocorr_peaks_chb = cc_and_acc_peaks(multich_acc[:,ch_b])
    autocorr_peaks_cha = cc_and_acc_peaks(multich_acc[:,ch_a])
    acc_n_peaks = (multich_acc[:,ch_b], autocorr_peaks_chb,
                   multich_acc[:,ch_a], autocorr_peaks_cha,)

    q_vector = np.ones(cc_peaks.size)*20
    good_cc_peaks = filter_cc_peaks_by_plausibility_and_minrkl(each_pair,
                                                (crosscor, cc_peaks),
                                                         q_vector,
                                                         array_geom,
                                                         fs)
    # 'centre' the cross-correlation peaks to get negative
    # and positive values for values on the left and right 
    # of the 0 lag in the centre.
    crosscor_peaks[each_pair] = good_cc_peaks-multich_cc[(ch_b,ch_a)].size*0.5
#%% What are the actual cross-cor peaks I need to find if it's the real sound source? 
actual_tdes = {}
for each_pair, _ in crosscor_peaks.items():
    m_a, m_b = each_pair
    
    delta_R = dist_mat[0,m_a] - dist_mat[0,m_b]
    tde = delta_R/vsound
    tde_samples = int(tde*fs)
    actual_tdes[each_pair] = tde_samples

actual_tdoas = np.zeros((4,4))
for pair, tdes in actual_tdes.items():
    actual_tdoas[pair[0],pair[1]] = tdes 

#%% Now create *all* pair TDOA peaks, using the relation 
# TDE_ba = -TDE_ab
comprehensive_crosscorpeaks = {}
for pair, tdes in crosscor_peaks.items():
    mic_x, mic_y = pair
    comprehensive_crosscorpeaks[(mic_y, mic_x)] = -1*tdes

# and then fuse the dictionary with the old one 
comprehensive_crosscorpeaks.update(crosscor_peaks)


#%% Make all triples and check their consistency
all_triples = list(combinations(range(4), 3))
# for each triple, try out all possible tdoas
tftm_seconds = 1e-3 # seconds
tftm_samples = int(tftm_seconds*fs)
good_triples = {}
good_triple_scores = {}

for each_triple in all_triples:
    mic_a, mic_b, mic_c = each_triple
    # the equation to test 
    # 0 = delay_ba + delay_bc + delay_ac 
    for tde_ba in comprehensive_crosscorpeaks[(mic_b, mic_a)]:
        for tde_bc in comprehensive_crosscorpeaks[(mic_b, mic_c)]:
            for tde_ac in comprehensive_crosscorpeaks[(mic_a, mic_c)]:
                consistency = -tde_ba + tde_bc - tde_ac
                trip_serialnum = 0
                if np.abs(consistency)<=tftm_samples:
                    key = str(f'{mic_a},{mic_b},{mic_c};{trip_serialnum}')
                    while key in list(good_triples.keys()):
                        trip_serialnum += 1 
                        key = str(f'{mic_a},{mic_b},{mic_c};{trip_serialnum}')
                    graph = np.zeros((3,3))
                    graph[1,0] = tde_ba; graph[0,1] = -tde_ba;
                    graph[1,2] = tde_bc; graph[2,1] = -tde_bc;
                    graph[0,2] = tde_ac; graph[2,0] = -tde_ac
                    df_graph = pd.DataFrame(graph)
                    df_graph.columns = [mic_a, mic_b, mic_c]
                    df_graph.index = [mic_a, mic_b, mic_c]
                    
                    consistency_score = gamma_tftm(consistency, tftm_samples)
                    good_triple_scores[key] = consistency_score
                    good_triples[key] = df_graph

#%%
def does_it_make_a_quadruple(three_triples):
    '''
    '''
    
    # NOT IMPLEMENTED HERE -- CHECK IF THE THREE TRIPLE ACTUALLY
    # ARE UNIQUE, and that there are no repeated triples!
    
    nodes = [set(each.columns) for each in three_triples]
    edges1,edges2,edges3 = [set(combinations(each,2)) for each in nodes]
    # make all common pairs between the available nodes and 
    common12 = edges1.intersection(edges2)
    common23 = edges2.intersection(edges3)
    common13 = edges1.intersection(edges3)
    common_edges = list(common12.union(common23).union(common13))
    # now check if all the common edges have the same values in at least two
    # of these triples 
    common_score = 0
    for common_edge in common_edges:
        node_a, node_b = common_edge
        edge_values = []
        for each_graph in three_triples:
            try:
                edge_values.append(each_graph.loc[node_a,node_b])
            except:
                pass
        if len(np.unique(edge_values))==1:
            common_score += 1 
    if common_score == len(common_edges):
        return True
    else: 
        return False

def make_quadruple_from_3_triples(triples):
    '''
    '''
    quadruple = np.zeros((4,4))
    diag_rows, diag_cols = np.diag_indices(4)
    quadruple[diag_rows,diag_cols] = np.nan

    graph = pd.DataFrame(quadruple)
    graph.columns = np.unique(np.concatenate([each.columns for each in triples]))
    graph.index = graph.columns
    
    for triple in triples:
        nodes = triple.columns.tolist()
        for node_a in nodes:
            for node_b in nodes:
                if node_a!=node_b:
                    graph.loc[node_a,node_b] = triple.loc[node_a,node_b]
    return graph

def calculate_consistency_of_quadruple(graph, **kwargs):
    '''
    
    '''
    triplets = combinations(range(4),3)
    # relation is tde_ab +tde_bc - tde_ac = 0
    total_consistency = 0
    for each in triplets:
        mic_a, mic_b, mic_c = each
        deviation = graph.loc[mic_a,mic_b]+graph.loc[mic_b,mic_c]-graph.loc[mic_a,mic_c]
        total_consistency += gamma_tftm(np.abs(deviation), kwargs.get('tftm',10))
    return total_consistency

#%% Choose the most consistent triple and proceed to make a quadruplet 
# do this multiple times until no more quadruplets can be formed.

best_starting_triple = list(good_triples.keys())[0]#[9]
best_triple = good_triples[best_starting_triple]
starting_mics = set(best_starting_triple.split(';')[0].split(','))

candidate_pool =  copy.deepcopy(good_triples)
candidate_pool.pop(best_starting_triple)

potential_additions = []
for each_triple, graph in candidate_pool.items():
    # if the triple shares 2 nodes and and edge - keep  it
    abc = set(each_triple.split(';')[0].split(','))
    common_mics = starting_mics.intersection(abc)
    two_mics_common = len(common_mics)==2
    if two_mics_common:
        # look for a common edge
        common_mic1, common_mic2 = [int(each) for each in common_mics]
        if graph.loc[common_mic1, common_mic2] == best_triple.loc[common_mic1, common_mic2]:
            # put it  into the list only if the same mic combination hasn't been
            # appended already 
            hits  = 0 
            for each in potential_additions:
                if np.allclose(each.columns, graph.columns):
                    hits += 1 
            if hits == 0:
                potential_additions.append(graph)
        

valid_quadruple = does_it_make_a_quadruple(potential_additions)
if valid_quadruple:
    iop = make_quadruple_from_3_triples(potential_additions)
    
    d_matrix = (iop.loc[:,0]/fs)*vsound
    d = d_matrix.to_numpy()[1:]
    array_geom_tracking = array_geom.copy()
    
    pos1, pos2 = sw02.spiesberger_wahlberg_solution(array_geom_tracking, d)
    print(pos1, pos2)
    print(calculate_consistency_of_quadruple(iop))

#%% now perform the TDOA re-projection to check if there were some phantom mirrors
positions = {0:pos1, 1:pos2}
position_with_positive_y = int(np.argwhere([pos1[1]>=0, pos2[1]>=0])[0])
sourcepos_estimate = positions[position_with_positive_y]
reprojected_tdoas = get_tdoa_matrix(sourcepos_estimate, array_geom)
reprojected_tdoas_sampels = (reprojected_tdoas/vsound)*fs

np.nansum(np.tril(np.abs(reprojected_tdoas_sampels-iop)))
