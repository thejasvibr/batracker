#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Localising sounds with the  Schau & Robinson 1987 algorithm
===========================================================

The Schau & Robinson algorithm formulates the localisatin problem into finding
a series of intersecting spheres, instead of hyperbolas. This also means that 
there are two possible solutions to the localisation. The user will need to choose
the relevant option. 

The input variables
-------------------
The Schau & Robinson paper defines a bunch of input variables:
    #. :math:`M`
    #. :math:`d`
    #. :math:`x`
    #. :math:`\Delta`  , where :math:`\Delta` itself is a matrix with
       :math:`{R_{i}}^2 - {d_{i4}}^2`. :math:`R_{i}` is the distance from 
       mic `i` to the origin (mic 4), and :math:`d_{i4}` is the range
       difference to the source between mic `i` and mic 4 (Eqn. 1), :math:`d_{i4}=D_{i}-D_{4}`


References
----------

[1] Schau, H. C., & Robinson, A. Z. (1987). Passive source localization employing intersecting spherical surfaces from time-of-arrival differences.
IEEE Transactions on Acoustics, Speech, and Signal Processing, 35(8), 1223-1225.

TODO:
    #. run tests with a simulation and then check if the answers are correct or not!

"""
import numpy as np



def schau_robinson_solution():
    '''
    Parameters
    ----------
    array_geometry
    TDOA
    v_sound
    '''

def parse_for_equation13(M, d, Delta):
    '''Converts the simple M,d,Delta inputs into the quadratic
    type forms.
    
    Parameters
    ----------
    M
    d
    Delta
    
    Returns
    -------
    a
    b
    c
    '''
    
    


def equation_13(a, b, c):
    '''
    Provides two solutions as it solves a quadratic solution to the range from source,
    :math:`R_{s} = \\frac{-b \pm \sqrt{b^2 - 4ac}}{2a}`. Being a quadratic equation, there
    will be two valid solutions, and the user will need to  choose the range solution 
    that is relevant (or not).
    
    Parameters
    ----------
    a : float
        Defined as :math:`4 - 4 d^T {M^{-1}}^{T} M^{-1} d`
    b : float
        Defined as :math:`2d^T {M^{-1}}^{T} M^{-1} \Delta + 2{\Delta}^{T}{M^{-1}}^{T}M^{-1} d`
    c : float
        Defined as :math:`-[{\Delta}^T {M^{-1}}^{T}M^{-1} \Delta]`
        
    
    
    '''
    Rs_plus = (-b + np.sqrt(b**2-4*a*c))/(2*a)
    Rs_minus = (-b - np.sqrt(b**2-4*a*c))/(2*a) 
    return Rs_plus, Rs_minus


def equation_10(M, Delta, Rs, d):
    '''
    Parameters
    ----------
    M : np.array
        3 x 3 array with the xyz coordinates of the microphones.
        All positions are relative to mic 4, as the last mic is considered
        to be the origin (0,0,0)
    Delta  : np.array
        Matrix, :math:`\Delta` in the paper, with :math:`{R_{i}}^2 - {d_{i4}}^2`, where is is from 1-3,
        and refers to mic number. 
        :math:`R_{i}` is the distance from mic `i` to the origin. 
        :math:`d_{i4}` is the difference in source to mic travel distance
        between mic `i` and mic 4 (See Eqn 1&2).
    Rs : tuple with 2 entries
        The two values of source range arising from the quadratic equation in Eq.13. 
    d : np.array
        The distances between each of the mics to the reference fourth mic. 
        Each distance is notated as :math:`d_{i4}` in the paper. 
    
    
    Returns
    -------
    X_1, X_2 : np.arrays
        The x,y,z coordinates of the source, with the origin of the coordinate 
        system set at mic 4. 


    Notes
    -----
    The :math:`R_{s}` will have two range solutions, which also means that there
    will be two solutions of X. The user must make a prudent choice about which values
    are relevant for the situation.   
    
    '''
    
    X_1 = 0.5*np.invert(M)*(Delta-2*Rs[0]*d)
    X_2 = 0.5*np.invert(M)*(Delta-2*Rs[1]*d)
    return X_1, X_2
