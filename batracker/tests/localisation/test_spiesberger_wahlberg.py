# -*- coding: utf-8 -*-
"""
Tests for Spiesberger & Wahlberg 2002
"""


import unittest 
import numpy as np 
np.random.seed(82319)
import scipy.spatial as spatial

from batracker.localisation import spiesberger_wahlberg_2002 as sw02

class SimpleTest(unittest.TestCase):
    '''With tristar and one position
    '''
    
    def setUp(self):
        # a tristar config 
       
        self.source_position = np.array([10,8,-2])

    def make_tristar_geom(self):
        R = 1.2
        theta = np.pi/3
        random_noise = np.random.normal(0,10**-6,4)
        array_geom =  np.array([[0,0.0,0],
                                     [-R*np.sin(theta),0.0,-R*np.cos(theta)],
                                     [ R*np.sin(theta),0.0,-R*np.cos(theta)],
                                     [0,0.0,R]])
        array_geom[:,1] = random_noise
        return array_geom

    def calculate_d_matrix(self, array_geom, source_position):
        # create source positions
        D_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom, 
                                       source_position)
        # %% Range difference , ref. Range of mic 4 (obtained from TDOAs)
        d_matrix = D_matrix - D_matrix[0]
        return d_matrix[1:].reshape(-1,1)
        

    def test_simple(self):
        #self.array_geom = self.array_geom - self.array_geom[0,:]
        self.array_geom = self.make_tristar_geom()
        self.d = self.calculate_d_matrix(self.array_geom, self.source_position)
        output_positions = sw02.spiesberger_wahlberg_solution(self.array_geom, self.d)
        num_matches = []
        for each in output_positions:
                num_matches.append(np.allclose(self.source_position, each, atol=10**-2))
        
        self.assertEqual(sum(num_matches), 1)

    def test_multiple_points(self):
        self.array_geom = self.make_tristar_geom()
        x_rand = np.random.choice(np.linspace(-10,10,50), 30)
        y_rand = np.random.choice(np.linspace(0.2,20,50), 30)
        z_rand = np.random.choice(np.linspace(-2,10,50), 30)
        
        multiple_points = np.column_stack((x_rand, y_rand, z_rand))
        
        all_calc_positions = []
        for i, position in enumerate(multiple_points):
            d = self.calculate_d_matrix(self.array_geom, position)
            output_positions = sw02.spiesberger_wahlberg_solution(self.array_geom, d)

            valid_position =  list(filter(lambda X: X[1]>0, output_positions))[0]
            all_calc_positions.append(valid_position)
        
        all_true = np.allclose(np.array(all_calc_positions).reshape(-1,3),
                               multiple_points, atol=1e-2)
        self.assertTrue(all_true)
        
    def test_withgreaterthan4mics(self):
        '''
        A 15 microphone array placed everywhere
        '''
        n_mics = 15
        num_sources = 20
        self.array_geom = np.random.normal(0,1,n_mics*3).reshape(n_mics,-1)
        self.array_geom *= 5
        
        x_rand = np.linspace(-10,13,num_sources)
        y_rand = np.linspace(0.2,20,num_sources)
        z_rand = np.linspace(-2,10,num_sources)
        for each in [x_rand, y_rand, z_rand]:
            np.random.shuffle(each)
    
        multiple_points = np.column_stack((x_rand, y_rand, z_rand))

        all_calc_positions = []
        for i in range(multiple_points.shape[0]):
            position = multiple_points[i,:]
            d = self.calculate_d_matrix(self.array_geom, position)
            output_positions = sw02.spiesberger_wahlberg_solution(self.array_geom, d)
            all_calc_positions.append(output_positions)
        
        match = []
        for i,(candidate1, candidate2) in enumerate(all_calc_positions):
            match_result = np.logical_or(np.allclose(candidate1, multiple_points[i,:], atol=1e-2),
                                         np.allclose(candidate2, multiple_points[i,:], atol=1e-2))
            match.append(match_result)

        self.assertTrue(np.all(match))
        

if __name__ == '__main__':
    unittest.main()      
        