#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing the STL arbit 3D volume simulation option in pyroomacoustics


Code is a slightly modified version of that on the Github repo 
(https://github.com/LCAV/pyroomacoustics/blob/pypi-release/examples/room_from_stl.py)

"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
import scipy.signal as signal 
import pyroomacoustics as pra
from stl import mesh

file_path = './INRIA_MUSIS.stl'


# with numpy-stl
the_mesh = mesh.Mesh.from_file(file_path)
ntriang, nvec, npts = the_mesh.vectors.shape
size_reduc_factor = 500.0  # to get a realistic room size (not 3km)



#%% Plot the mesh 
# (code form https://numpy-stl.readthedocs.io/en/latest/usage.html)

figure = plt.figure()
axes = mplot3d.Axes3D(figure)

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(the_mesh.vectors))

# Auto scale to the mesh size
scale = the_mesh.points.flatten()
axes.auto_scale_xyz(scale, scale, scale)
scale = 0.01

# Show the plot to the screen
plt.show()

#%% Now state the wall specs

material = pra.Material(energy_absorption=0.2, scattering=0.1)

# create one wall per triangle
walls = []
for w in range(ntriang):
    walls.append(
        pra.wall_factory(
            the_mesh.vectors[w].T / size_reduc_factor,
            material.energy_absorption["coeffs"],
            material.scattering["coeffs"],
        )
    )
#%% Initialise the room - this seems to take some time.
durn = 0.001
fs = 192000
t = np.linspace(0,durn,int(fs*durn))
source_signal = signal.chirp(t, 20000, t[-1], 85000)

room = (
    pra.Room(walls, fs=fs, max_order=2, ray_tracing=True, air_absorption=True,)
    .add_source([-2.0, 2.0, 1.8], signal=source_signal, delay=0.005)
    .add_microphone_array(np.c_[[-6.5, 8.5, 3 + 0.1], [-6.5, 8.1, 3 + 0.1]])
)
#%% 
# compute the rir
room.image_source_model()
room.ray_tracing()
room.compute_rir()
room.plot_rir()

# show the room
room.plot(img_order=1)
plt.show()
#%% 
# Run sound simulations 

# simulate sound propagati"on 
room.simulate()
#%% Also show the spectrograms of the channels to understand how the sound propagated
plt.figure()
a0=plt.subplot(211)
plt.specgram(room.mic_array.signals[0, :], Fs=fs)
plt.ylim(15000,fs*0.5)

plt.subplot(212,sharex=a0,sharey=a0)
plt.specgram(room.mic_array.signals[1, :], Fs=fs)
plt.ylim(15000,fs*0.5)

