# -*- coding: utf-8 -*-
"""
Spiesberger & Wahlberg 2002
===========================
A neat formulation to implement - the original paper also formulates 
a variable speed of sound propagation at each receiver that is NOT 
implemented here. 

This formulation handles >=4 microphones

Reference
---------
1. Spiesberger & Wahlberg 2002, Probability density functions for hyperbolic and isodiachronic locations, 
   JASA, 112, 3046 (2002); doi: 10.1121/1.1513648

Notes
-----
If the microphones are perfectly in plane (like in a tristar array), the 
algorithm fails. To make it work, a small 'jitter' should be added to the axis with 
common coordinate values. 

"""
import numpy as np 
import scipy.spatial as spatial
matmul = np.matmul
from batracker.localisation.tdoa_simulation import generate_tdoa_predictions

def spiesberger_wahlberg_solution(array_geometry, d, **kwargs):
    '''
    Parameters
    ----------
    array_geometry: np.array
        A 4x3 array with xyz coordinates of >= 4 mics.
        The first mic will be taken as the reference microphone.
    d:  np.array
        A 3x1 np.array with the range_differences to the source. 
        All range_differences (:math:`d_{0..N}`), are calculated by
        taking :math:`D_{i}-D_{0}`, where :math:`D` is the direct
        range from a mic to the source.

    Returns
    -------
    s : list
        A list with two 3x1 np.arrays with the x,y,z positions of the source.
        The two xyz positions describe two possible solutions to the given 
        array geometry and range differences.
    
    Notes
    -----
    The first mic in this formulation must be the origin (0,0,0). If array_geometry
    doesn't have the first mic's position as 0,0,0 - this is taken care of using 
    relative subtraction and addition. 

    '''
    c = kwargs.get('c', 338.0) # m/s
    # check that the 1st mic is origin - else set it to 0,0,0
    if not np.array_equal(array_geometry[0,:], np.array([0,0,0])):
        mic1_notorigin = True
        mic1 = array_geometry[0,:]
        array_geometry = array_geometry - mic1
    else:
        mic1_notorigin = False
        
    # the receiver matrix- excluding the first channel.
    R = array_geometry[1:,:]
    tau = d.copy()/c # to keep symbol conventions the same
    
    try:
        R_inv = np.linalg.inv(R)
    except:
        R_inv = np.linalg.pinv(R)
    
    Nrec_minus1 = R.shape[0]
    b = np.zeros(Nrec_minus1)
    f = np.zeros(Nrec_minus1)
    g = np.zeros(Nrec_minus1)
    for i in range(Nrec_minus1):
        b[i] = np.linalg.norm(R[i,:])**2 - (c*tau[i])**2
        f[i] = (c**2)*tau[i]
        g[i] = 0.5*(c**2-c**2)
    
    
    a1 = matmul(matmul(R_inv, b).T, matmul(R_inv,b))
    a2 = matmul(matmul(R_inv, b).T, matmul(R_inv,f))
    a3 = matmul(matmul(R_inv, f).T, matmul(R_inv,f))
    
    # quadratic equation is ax^2 + bx + c = 0
    # solution is x = (-b +/- sqrt(b^2 - 4ac))/2a
    # replace 
    
    a_quad = a3 - c**2
    b_quad = -a2
    c_quad = a1/4.0
    
    t_solution1 = (-b_quad + np.sqrt(b_quad**2 - 4*a_quad*c_quad))/(2*a_quad)
    t_solution2 = (-b_quad - np.sqrt(b_quad**2 - 4*a_quad*c_quad))/(2*a_quad)
    t1 = (t_solution1 , t_solution2)
    
    
    
    s = [matmul(R_inv,b*0.5) - matmul(R_inv,f)*t1[0],
         matmul(R_inv,b*0.5) - matmul(R_inv,f)*t1[1]]
    if mic1_notorigin:
        for each in s:
            each += mic1
    return s


def sw_tau51_better_solution(obs_rangediffs, sources, micgeom, vsound=340):
    '''
    Choosing the correct solution from SW2002 based on the TDOA of the last channel pair. 
    
    Parameters
    ----------
    obs_rangediffs : (Nmics-1,) np.array
        The rangedifferences with ref the first channel in meters. 
    sources : list with Nsources entries. Each source is a (3,) np.array
        List with candidate sources for a given range difference
    micgeom : (Nmics,3) np.array
    vsound : float>0, optional
        Defaults to 340 m/s
    
    Notes
    -----
    According to Spiesberger & Wahlberg 2002: 'For each ambiguous source location, one
    can generate a model for t 51 and choose the root for t1 that
    yields a model for t 51 that is closest to that measured' - where they talk about 
    the 5 channel case allowing the selection of the correct source location. 
    
    Here I've generalised it to include the last channel pair in general - and simulation tests
    indicate this does a good job of choosing the correct source. 
    
    '''
    minerror_in_prediction = []
    nmics = micgeom.shape[0]
    if nmics <=4:
        raise ValueError(f'{nmics} mics detected - SW2002 only gives reliable selection >=5 channels!')
    for i, source in enumerate(sources):
        pred_tdoas = generate_tdoa_predictions(source, micgeom)
        pred_tdoas = pred_tdoas[:micgeom.shape[0]-1]
        pred_rangediff = pred_tdoas*vsound
        error_index = pred_rangediff[-1] - obs_rangediffs[-1]
        minerror_in_prediction.append(abs(error_index))
    better_solution = sources[np.argmin(minerror_in_prediction)]
    return better_solution

def sw_choose_robust_better_solution(obs_rangediffs, sources, micgeom, vsound=343.0):
    '''
    The 'robustness' comes from the fact that the TDOAs are not just checked
    for overall error - but also for flipped 'polarities'...TDOA & -1*TDOA
    
    Parameters
    ----------
    obs_rangediffs : (N,) np.array
        Nmics-1 range differences w.r.t reference mic. 
    sources : list with (3,) np.arrays
        Candidate sources predicted for a given array geometry
    micgeom : (Nmics,3) np.array
        Microphone geometry xyz coordinates. One mic per row.
    vsound: float>0, optional
        Speed of sound. Defaults to 343 m/s

    Returns
    -------
    better_solution : (3,) np.array
        xyz coordinates of the best solution. 
    '''
    minerror_in_prediction = []
    for i, source in enumerate(sources):
        pred_tdoas = generate_tdoa_predictions(source, micgeom)
        pred_tdoas = pred_tdoas[:micgeom.shape[0]-1]
        pred_rangediff = pred_tdoas*vsound
        polarities = np.array([1,-1])
        tdoas_polarities = polarities*pred_rangediff
        tdoas_deviations = tdoas_polarities - obs_rangediffs.reshape(-1,1)
        error_index = np.max(abs(tdoas_deviations), axis=0)
        minerror_in_prediction.append(np.min(error_index))
    better_solution = sources[np.argmin(minerror_in_prediction)]
    return better_solution

if __name__ == '__main__':
    
    R = 1.2 # meters
    theta = np.pi/3
    other_x_position = 0.5
    theta2 = np.arctan(other_x_position/(R*np.cos(theta)))
    R_2 = np.sqrt(other_x_position**2 +  (R*np.cos(theta))**2)
    arbit_y = 10**-6
    Ro = np.array([[0,arbit_y,0],
                    [R_2*np.sin(theta2),  arbit_y, -R*np.cos(theta), ],
                    [-R*np.sin(theta), arbit_y, -R*np.cos(theta)],
                    [1,1.2,1],
                    [1.0,arbit_y,R]])
    
    Ro = Ro-Ro[0,:]
    Ro[:,1] += arbit_y

    source_pos = np.array([1,14.25,-15])
    d_matrix = np.zeros(Ro.shape[0])
    for i in range(d_matrix.size):
        d_matrix[i] = spatial.distance.euclidean(Ro[i,:], source_pos)
    d = (d_matrix[1:] - d_matrix[0])
    #d = np.array([0.00128646, 0.00266667, 0.00220833])*338
    outputs = spiesberger_wahlberg_solution(Ro, d)
    print(outputs)
    
    better_solution = sw_choose_robust_better_solution(d, outputs, Ro)
    print('\n',source_pos, better_solution)
    
    #%%
    sourcepos = np.array([7.79179179, -1.86586587  ,6.302302])
    micgeom = np.array([[-3.13057326,  3.34174603,  3.56756099],
                     [-9.73651937, -5.29910719,  3.19106483],
                     [-9.09807254,  9.86464726, -1.3601096 ],
                     [-1.50989502, -6.33514563, -4.63342018],
                     [ 0.6214693,   4.81233551, -1.56560442]])
    from batracker.localisation.tdoa_simulation import generate_tdoa_predictions
    tdoas = generate_tdoa_predictions(sourcepos, micgeom, vsound=340)
    nmics = micgeom.shape[0]
    rangediff = tdoas[:nmics-1]*340
    output_positions = spiesberger_wahlberg_solution(micgeom, rangediff)
    better_solution = sw_choose_robust_better_solution(rangediff,
                                                            output_positions,
                                                            micgeom)
    bet_soln = sw_tau51_better_solution(rangediff, output_positions, 
                                        micgeom)
    print('actual source:', sourcepos )
    for source in output_positions:
        tdoas_source = generate_tdoa_predictions(source, micgeom, vsound=340)
        print(source, tdoas[nmics-1], tdoas_source[nmics-1])
    
    
    
    