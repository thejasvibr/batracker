#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Understanding some of the things in the Scheuing &  Yang 2008 paper. 

This is an *extremely* simplified situation where there is one direct
path and one indirect path that is recorded by all microphones. 

Created on Mon Jul  5 10:48:01 2021

@author: thejasvi
"""
#%%
import numpy as np 
np.random.seed(82319)
import scipy.signal  as signal 
import scipy.spatial as spatial 
import matplotlib.pyplot as plt 
from itertools import combinations

#%%

def simulate_sound_propagation(**kwargs):
    # microphone array geometry
    R = 1.2
    theta = np.pi/3
    tristar = np.row_stack(([0,0,0],
                            [-R*np.sin(theta), 0, -R*np.cos(theta)],
                            [R*np.sin(theta), 0, -R*np.cos(theta)],
                            [0,0, R]))
    
    sound_pos = kwargs.get('sound_pos',np.array([3,2,1]))
    reflection_source = kwargs.get('reflection_source',np.array([1,4,1]))
    direct_indirect_sources = np.row_stack((sound_pos, reflection_source))
    # direct path propagation:
    
    dist_mat = spatial.distance_matrix(direct_indirect_sources, tristar)
    #add the distance of propagation from source to reflection point
    source_to_reflectionpoint = spatial.distance.euclidean(sound_pos, reflection_source)
    dist_mat[1,:] += source_to_reflectionpoint
    
    # make the direct
    
    chirp_durn = 0.003
    fs = 192000
    t = np.linspace(0,chirp_durn,int(fs*chirp_durn))
    chirp = signal.chirp(t,80000,t[-1],25000)
    chirp *= signal.hann(chirp.size)*0.5
    
    
    vsound = 340.0
    audio = np.zeros((int(fs*0.03),4))
    toa_sounds = dist_mat/vsound
    toa_samples = np.int64(toa_sounds*fs)
    for channel in range(4):
        random_atten = np.random.choice(np.linspace(0.2,0.9,20),2)
        start_direct, start_indirect = toa_samples[0,channel], toa_samples[1,channel]
        audio[start_direct:start_direct+chirp.size,channel] += chirp*random_atten[0]
        audio[start_indirect:start_indirect+chirp.size,channel] += chirp*random_atten[1]
    audio += np.random.normal(0,1e-5,audio.size).reshape(audio.shape)
    return audio , dist_mat, tristar


#%% Generate the cross-corr for each channel pair

def generate_multich_crosscorr(input_audio):
    '''
    Generates all unique pair cross-correlations: (NxN-1)/2 pairs. Each pair is
    designated by a tuple where the second number is the reference channel, eg. (1,0)
    where channel 1 is cross-correlated with reference to channel 0. 

    Parameters
    ----------
    input_audio: np.array
        M samples x N channels

    Returns
    -------
    multichannel_cc : dictionary 
        Keys indicate the channel pair, and entries are the cross-correlation. 
        Each cross-correlation has M samples (same size as one audio channel).
    '''
    num_channels = input_audio.shape[1]
    unique_pairs = list(combinations(range(num_channels), 2))
    multichannel_cc = {}
    for cha, chb in unique_pairs:
        # make sure the lower channel number is the reference signal 
        signal_ch, ref_signal = sorted([cha, chb], reverse=True)
        multichannel_cc[(signal_ch, ref_signal)] = signal.correlate(input_audio[:,signal_ch],
                                                                 input_audio[:,ref_signal],'same')
    return multichannel_cc

def generate_multich_autocorr(input_audio):
    '''
    
    Parameters
    ----------
    input_audio : np.array
        M samples x Nchannels
    
    Returns 
    -------
    multichannel_autocor : np.array
        M samples x Nchannels
    '''
    return np.apply_along_axis(lambda X: signal.correlate(X,X,'same'),0, input_audio)

##%% extract auto-cor and cross-cor peaks 
cc_and_acc_peaks = lambda X:  signal.find_peaks(X, 0.11,distance=int(fs*1e-4))[0] 

def make_p_prime_kk_ll(crosscor_peaks, acc_peaks_ch1, acc_peaks_ch2, twrm):
    '''
    Identifies the peaks in the cross-correlation which are from echo paths. 
    
    Parameters
    ----------
    crosscor_peaks : np.array
        With indices of peak locations 
    acc_peaks_ch1, acc_peaks_ch2: np.array
        With centred indices of peak locations of the autocorrelation. 
        'Centred indices' mean that the 'centre' of the signal is the 0th 
        index, all those to the left of it have -ve indices and to the right 
        have +ve. 

    Returns 
    -------
    p_prime_kk, p_prime_ll: dict
        Dictionaries with the echo path delay as key and the corresponding cross-cor
        peaks as entries. 
    '''
    peak_pairs_that_match_acc = []
    p_prime_kk = {}
    p_prime_ll = {}
    eta_mu = np.concatenate((acc_peaks_ch1, acc_peaks_ch2))

    for focal_peak in crosscor_peaks:
        non_focalpeaks = np.argwhere(crosscor_peaks!=focal_peak)
        cc_peaks_wo_focal = crosscor_peaks[non_focalpeaks]
        for each in cc_peaks_wo_focal:
            difference01 = each-focal_peak
            acc_match = np.abs(eta_mu-np.abs(difference01))
    
            if np.any(acc_match<=twrm*0.5):
                peak_pairs_that_match_acc.append(focal_peak)
                peak_pairs_that_match_acc.append(each.tolist()[0])
                
                # save the acc peak into the P_prime list
                acc_delay = eta_mu[np.argmin(acc_match)]
                # figure out if it's from acc ch1/ch2
                correct_channel_acc =  [acc_delay in acc_peaks_ch2, acc_delay in acc_peaks_ch1]
                if np.logical_and(sum(correct_channel_acc)>0, sum(correct_channel_acc)<2):

                    if np.argmax(correct_channel_acc)==1:
                        p_prime_ll[acc_delay] = [focal_peak, each.tolist()[0]]
                    elif np.argmax(correct_channel_acc)==0:
                        p_prime_kk[acc_delay] = [each.tolist()[0], focal_peak, ]
                else:
                    ValueError('Problem with finding correct acc channel ')

    return p_prime_kk, p_prime_ll
 
def gamma_tfrm(autocorr_delay, cc_peak_diff, tfrm):
    difference = np.abs(autocorr_delay) - np.abs(cc_peak_diff)
    if difference < 0.5*tfrm:
        return 1 - (difference/(0.5*tfrm))
    else:
        return 0

def calculate_quality_value(crosscor_and_peaks, p_primekk_ll, acc_and_peaks, twrm):
    '''
    
    Parameters
    ----------
    crosscor_and_peaks : tuple
        With (cross_cor, crosscor_peaks)
    p_primekk_ll : tuple
        With (p_prime_kk, p_prime_ll)
    acc_and_peaks : tuple
        With (acc_channel1, acc1_peaks, acc_channel2, acc2_peaks).
        The autocorrelation peaks here are expected to be 'centred'. 
    twrm : int>0
        The tolerance width of raster matching in samples. 
    
    Returns 
    -------
    quality_values : np.array
        Same size as the number of crosscor_peaks, with the corresponding 
        quality score. 

    '''
    cross_cor, cc_peaks = crosscor_and_peaks
    p_primekk, p_primell = p_primekk_ll
    acc_ch1, acc1_peaks, acc_ch2, acc2_peaks = acc_and_peaks

    quality_values = np.zeros(cc_peaks.size)
    # where peak1 = eta_{gamma} and peak2 = eta_{mu}

    for i,each_ccpeak in enumerate(cc_peaks):
        rkl_mu = cross_cor[each_ccpeak]
        ch1_autocorr_term = 0
        ch2_autocorr_term = 0
        
        for acc_delay, (peak1, peak2) in p_primekk.items():
            
            uncentred = int(acc_ch1.size*0.5)+int(acc_delay)
            acc_peak_value = acc_ch1[uncentred]
            cc_delay = peak1-peak2
            thisdelay_autocorr_term = np.sign(peak1-peak2)*np.abs(acc_peak_value)
            thisdelay_autocorr_term *= gamma_tfrm(acc_delay, cc_delay, twrm)
            ch1_autocorr_term += thisdelay_autocorr_term

        for acc_delay, (peak1, peak2) in p_primell.items():
            uncentred = int(acc_ch2.size*0.5)+int(acc_delay)
            acc_peak_value = acc_ch2[uncentred]
            cc_delay = peak1-peak2
            thisdelay_autocorr_term = np.sign(peak1-peak2)*np.abs(acc_peak_value)
            thisdelay_autocorr_term *= gamma_tfrm(acc_delay, cc_delay, twrm)
            ch2_autocorr_term += thisdelay_autocorr_term

        quality_values[i] = rkl_mu + ch1_autocorr_term + ch2_autocorr_term

    return quality_values

def filter_cc_peaks_by_plausibility_and_minrkl(mic_pair, crosscor_and_peaks, quality_values,
                                               array_geom, fs, **kwargs):
    '''
    
    '''
    mic_1, mic_2 = mic_pair
    cross_cor, cc_peaks = crosscor_and_peaks
    min_rkl = np.min(cross_cor[cc_peaks])
    quality_peaks = cc_peaks[cc_peaks>=min_rkl]
    
    intermic_distance = spatial.distance.euclidean(array_geom[mic_1,:], 
                                                   array_geom[mic_2,:])
    vsound = kwargs.get('vsound', 338.0)
    max_intermic_delay = intermic_distance/vsound
    centred_peaks = quality_peaks - cross_cor.size*0.5
    peak_delays = centred_peaks/fs
    
    return quality_peaks[peak_delays<=max_intermic_delay]
#%%
    
def gamma_tftm(delay, tftm, **kwargs):
    '''
    '''
    if np.abs(delay)>=tftm:
        return 0
    else:
        
        if kwargs.get('weighting_function') is None:
            weighting_function = np.cos
            theta = 0.25*np.pi*(delay/(tftm*0.5))
            return weighting_function(theta)
        else:
            return weighting_function(delay, tftm)
#%% 
    def cyclic_tde_sum(tde_ba, tde_bc, tde_ac):
        return tde_ba+tde_bc+tde_ac


    def calculate_connectivity():
        '''
        the 'w' term described in equation 24, where broadly speaking:

            :math:`w = \Sigma_{all\:consistent\:triples} \Gamma_{TFTM}(cyclic \:sum \:of \:triple)
        
            
        
        '''
    
    def parse_triplet_graph_key(triplet_key):
        '''
        parses and splits the input string of the following format:
        
        
        'micA,micB,micC;tde_ba,tde_bc,tde_ac'
        Here micA,B,C are >=0 integers, while tde_ba,_bc,_ac are 
        also integers that can be <=0 or >0. 
        
        
        Parameters
        ----------
        triplet_key : str
            See description.
        
        Returns
        -------
        (mica,micb,micc),(tdeba,tdebc,tdeac): str
        '''
        mic_ids, tde_values = triplet_key.split(';')
        mica, micb, micc = mic_ids.split(',')
        tde_a, tde_b, tde_c = tde_values.split(',')
        return (mica, micb, micc), (tde_a, tde_b, tde_c)
#%%
                
    def find_quadruplet_from_triplet_set(starter_triplet, candidate_triplets, **kwargs):
        '''
        The candidate triplet list is first cut down to show only those whose 
        nodes have at least one different microphone. 
        
        Then all candidates are checked for two common nodes and the same edge length. 
        
        
        '''
        # remove all candidate triplets with the exact same mic numbers     
        subset_of_triplets = remove_exactly_same_triplets(starter_triplet,
                                                                candidate_triplet)
        if not len(subset_of_triplets)>0:
            return None
            
        for candidate_triplet in subset_of_triplets:
            commonality_found = check_for_common_edge_and_nodes(starter_triplet, candidate_triplet)
            if commonality_found: 
                quadruplet = fuse_two_triplets(starter_triplet, candidate_triplet)
                return quadruplet

#%% 
if __name__ == '__main__':
    #%%
    # make simulated audio and generate cross-cor + auto-cor
    vsound = 338.0
    
    sim_audio, dist_mat, array_geom = simulate_sound_propagation()
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
        pprime_kk_ll = make_p_prime_kk_ll(cc_peaks, 
                                          autocorr_peaks_chb,
                                          autocorr_peaks_cha,
                                          twrm_samples)
        acc_n_peaks = (multich_acc[:,ch_b], autocorr_peaks_chb,
                       multich_acc[:,ch_a], autocorr_peaks_cha,)
        q_vector = calculate_quality_value((crosscor,cc_peaks),
                                           pprime_kk_ll,
                                           acc_n_peaks,
                                           twrm_samples)
        good_cc_peaks = filter_cc_peaks_by_plausibility_and_minrkl(each_pair,
                                                    (crosscor, cc_peaks),
                                                             q_vector,
                                                             array_geom,
                                                             fs)
        # 'centre' the cross-correlation peaks to get negative
        # and positive values for values on the left and right 
        # of the 0 lag in the centre.
        crosscor_peaks[each_pair] = good_cc_peaks-sim_audio.shape[0]*0.5
    #%% Now create *all* pair TDOA peaks, using the relation 
    # TDE_ba = -TDE_ab
    comprehensive_crosscorpeaks = {}
    for pair, tdes in crosscor_peaks.items():
        mic_x, mic_y = pair
        comprehensive_crosscorpeaks[(mic_y, mic_x)] = -1*tdes
    
    # and then fuse the dictionary with the old one 
    comprehensive_crosscorpeaks.update(crosscor_peaks)
    
    #%% 
    all_triples = list(combinations(range(4), 3))
    # for each triple, try out all possible tdoas
    tftm_seconds = 1.5e-3
    tftm_samples = int(tftm_seconds*fs)
    consistent_triples = {}
    
    for each_triple in all_triples:
        mic_a, mic_b, mic_c = each_triple
        
        # the equation to test 
        # 0 = delay_ba + delay_bc + delay_ac 
        for tde_ba in comprehensive_crosscorpeaks[(mic_b, mic_a)]:
            for tde_bc in comprehensive_crosscorpeaks[(mic_b, mic_c)]:
                for tde_ac in comprehensive_crosscorpeaks[(mic_a, mic_c)]:
                    consistency = tde_ba + tde_bc + tde_ac
                    if np.abs(consistency)<=tftm_samples:
                        key = str(f'{mic_a},{mic_b},{mic_c};{tde_ba},{tde_bc},{tde_ac}')
                        consistency_score = gamma_tftm(consistency, tftm_samples)
                        consistent_triples[key] = consistency_score
    #%% Choose the most consistent triple and proceed to make a quadruplet 
    import copy
    values = sorted(list(consistent_triples.values()))
    sorted_values = sorted(values, reverse=True)
    keys = list(consistent_triples.keys())
    most_consistent_triplet = keys[values.index(sorted_values[0])]
    # now remove this triplet and all other triplets with the same mics in them 
    candidate_triplets = copy.deepcopy(consistent_triples)
    candidate_triplets.pop(most_consistent_triplet)
    
    candidate_triplet_graphs = list(candidate_triplets.keys())
    
    #%% Proceed iteratively and generate 
    
