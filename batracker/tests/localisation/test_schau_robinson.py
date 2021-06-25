# -*- coding: utf-8 -*-
"""
Schau robinson tests

"""


import unittest 
import numpy as np 
import scipy.spatial as spatial

from batracker.localisation import schau_robinson_1987 as sr1987

class SimpleTest(unittest.TestCase):
    
    def setUp(self):
        # a tristar config 
        R = 1.2
        theta = np.pi/3
        self.array_geom =  np.array([[0,0.01,0],
                           [-R*np.sin(theta),0.1,-R*np.cos(theta)],
                           [ R*np.sin(theta),0.1,-R*np.cos(theta)],
                           [0,0.1,R]])
        self.source_position = np.array([10,8,-2])
        self.calculate_d_matrix()

    def calculate_d_matrix(self):
        # create source positions
        
        source_position_ref_m4 = self.source_position-self.array_geom[-1,:]
        self.final_sourcepositions = source_position_ref_m4.reshape(3,1)
        array_geom_ref_m4 = self.array_geom - self.array_geom[-1,:]
        
        D_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom_ref_m4,
                                   source_position_ref_m4)
        # %% Range difference , ref. Range of mic 4 (obtained from TDOAs)
        d_matrix = D_matrix - D_matrix[-1]
        self.d = d_matrix[:-1].reshape(3,1)
        
        
    def test_simple(self):
        
        output_positions = sr1987.schau_robinson_solution(self.array_geom, self.d)
        print(self.source_position, output_positions)


if __name__ == '__main__':
    unittest.main()
        
        