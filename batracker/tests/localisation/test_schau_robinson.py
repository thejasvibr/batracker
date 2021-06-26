# -*- coding: utf-8 -*-
"""
Schau robinson tests

"""


import unittest 
import numpy as np 
np.random.seed(82319)
import scipy.spatial as spatial

from batracker.localisation import schau_robinson_1987 as sr87

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
        array_geom =  np.array([[0,0.0,R],
                                     [-R*np.sin(theta),0.0,-R*np.cos(theta)],
                                     [ R*np.sin(theta),0.0,-R*np.cos(theta)],
                                     [0,0.0,0]])
        array_geom[:,1] = random_noise
        return array_geom

    def calculate_d_matrix(self, array_geom, source_position):
        # create source positions
        D_matrix = np.apply_along_axis(spatial.distance.euclidean, 1,array_geom, 
                                       source_position)
        # %% Range difference , ref. Range of mic 4 (obtained from TDOAs)
        d_matrix = D_matrix - D_matrix[-1]
        return d_matrix[:-1].reshape(-1,1)
        

    def test_simple(self):
        #self.array_geom = self.array_geom - self.array_geom[0,:]
        self.array_geom = self.make_tristar_geom()
        self.d = self.calculate_d_matrix(self.array_geom, self.source_position)
        output_positions = sr87.schau_robinson_solution(self.array_geom, self.d)
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
            output_positions = sr87.schau_robinson_solution(self.array_geom, d)

            valid_position =  list(filter(lambda X: X[1]>0, output_positions))[0]
            all_calc_positions.append(valid_position)
        
        all_true = np.allclose(np.array(all_calc_positions).reshape(-1,3),
                               multiple_points, atol=1e-2)
        
    def test_withgreaterthan4mics(self):
        '''
        A 15 microphone array placed everywhere
        '''
        n_mics = 15
        num_sources = 1
        array_geom = np.random.normal(0,1,n_mics*3).reshape(n_mics,-1)
        position = np.random.normal(0,1,3)            
        d = self.calculate_d_matrix(array_geom, position)
        try:
            output_positions = sr87.schau_robinson_solution(array_geom, d)
            error_raised = False
        except ValueError:
            error_raised = True
        
        self.assertTrue(error_raised)

if __name__ == '__main__':
    unittest.main()
        
        