# -*- coding: utf-8 -*-
"""
Tests for Spiesberger & Wahlberg 2002
"""


import unittest 
import numpy as np 
import scipy.spatial as spatial

from batracker.localisation import spiesberger_wahlberg_2002 as sw02

class SimpleTest(unittest.TestCase):
    
    def setUp(self):
        # a tristar config 
        R = 1.2
        theta = np.pi/3
        random_noise = np.random.normal(0,10**-6,4)
        self.array_geom =  np.array([[0,0.0,0],
                                     [-R*np.sin(theta),0.0,-R*np.cos(theta)],
                                     [ R*np.sin(theta),0.0,-R*np.cos(theta)],
                                     [0,0.0,R]])
        self.array_geom[:,1] = random_noise
        self.source_position = np.array([10,8,-2])
        self.calculate_d_matrix()

    def calculate_d_matrix(self):
        # create source positions
        D_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,self.array_geom, 
                                       self.source_position)
        # %% Range difference , ref. Range of mic 4 (obtained from TDOAs)
        d_matrix = D_matrix - D_matrix[0]
        self.d = d_matrix[1:].reshape(3,1)
        
    def test_simple(self):
        self.array_geom = self.array_geom - self.array_geom[0,:]
        output_positions = sw02.spiesberger_wahlberg_solution(self.array_geom, self.d)
        print(output_positions)
        num_matches = []
        for each in output_positions:
                num_matches.append(np.allclose(self.source_position, each, atol=10**-2))
        
        self.assertEqual(sum(num_matches), 1)

if __name__ == '__main__':
    unittest.main()
        
        