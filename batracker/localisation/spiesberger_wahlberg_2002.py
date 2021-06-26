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
    print(spiesberger_wahlberg_solution(Ro, d))
