#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Test trajectory #1
==================
Create some example audio files with a trajectory in it. 


"""
import os
import numpy as np 
import pandas as pd
import yaml


## % 
# Make the required source files for the simulated audio 
 
mic_geometry  = np.array([[0, 0, 0],
                          [1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]])

source_x = np.linspace(-10,10,5)

source_positions = np.array(np.meshgrid(source_x, source_x, source_x,)).T.reshape(-1,3)

# save the array geometry file 
array_geom_file = 'array_geom.csv'
source_pos_file = 'source_pos.csv'
pd.DataFrame(mic_geometry, columns=['x','y','z']).to_csv(array_geom_file)
pd.DataFrame(source_positions, columns=['x','y','z']).to_csv(source_pos_file)


yaml_entries = {}
yaml_entries['array_geometry'] = array_geom_file
yaml_entries['source_position'] = source_pos_file
yaml_entries['sample_rate'] = 250000
yaml_entries['sim_name'] = 'batracker_simple'

with open('sim_audio_params.yaml', 'w') as f:
    yaml.dump(yaml_entries, f)
    
    
## %
# Call tacost and get the audio files
import os
stream = os.popen('bash -i make_sim_audio.sh')
output = stream.read()
output
