#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prototyping module: Schau & Robinson 1987
=========================================
This module contains test examples to check if things are working and if they are 
accurate.
"""

import numpy as np 
import scipy.spatial as spatial

# this needs to be fixed at some point...
import sys 
sys.path.append('../../')
from localisation import schau_robinson_1987 as sr1987
# %% 
# Let's imagine a 4 channel array somewhere with a couple of sound sources around it.
# To check if my Schau & Robinson implementation is working I need a couple of things
# 1) range diffferences between the mics to sources
# 2) array geometry

# M matrix
raw_array_geom = np.array([[0,0,1],
                           [0,1,0],
                           [2,0,0],
                           [0.5,0.5,0]]) # each row is one mic

array_geom_ref_m4 = raw_array_geom-raw_array_geom[-1,:]
M = array_geom_ref_m4[:-1,:]

# %% Making the simulated data
source_position = np.random.rand(3)*np.random.choice([1,10,100],1)
source_position_ref_m4 = source_position-raw_array_geom[-1,:]
x = source_position_ref_m4.reshape(3,1)


#%% Range to source (unknown in experiment)
D_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom_ref_m4,
                               source_position_ref_m4)
# %% Range difference , ref. Range of mic 4 (obtained from TDOAs)
d_matrix = D_matrix - D_matrix[-1]
d = d_matrix[:-1].reshape(3,1)



#%% Mic to origin distances
R_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom_ref_m4,
                               np.zeros(3))

R = (R_matrix[:-1]).reshape(3,1)


# %% Delta matrix
Delta = R**2-d**2



# %% Solutions!!

# %% Equation 10
M_inv = np.linalg.inv(M)
Minv_transp_into_Minv = (M_inv.T).dot(M_inv)

a = 4 - 4*(d.T).dot(Minv_transp_into_Minv.dot(d))

b_leftterm = 2*(d.T).dot(Minv_transp_into_Minv.dot(Delta))
b_rightterm = 2*(Delta.T).dot(Minv_transp_into_Minv.dot(d))
b = b_leftterm + b_rightterm 

c = -( (Delta.T).dot(Minv_transp_into_Minv.dot(Delta)))

num1 = -b + np.sqrt(b**2-4*a*c)
num2 = -b - np.sqrt(b**2-4*a*c)
denom = 2*a

Rs1 = num1/denom
Rs2 = num2/denom
Rs = (Rs1, Rs2)
#print(Rs, D_matrix[-1])

x_solutions = []
for R_solution in Rs:
    x = 0.5*M_inv.dot(Delta-2*R_solution*d)
    x_solutions.append(x)

# final coordinates in real world frame of reference
real_world_locations = []
for each in x_solutions:
    real_world_locations.append(each.flatten()+raw_array_geom[-1,:].T)
#print(real_world_locations)

print(f'\n Expected \n: {real_world_locations}')


#%% Acctually testing the module all by itself

solutions = sr1987.schau_robinson_solution(raw_array_geom, d)
print(f'\n Obtained \n: {solutions}')
