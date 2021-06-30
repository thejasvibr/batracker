'''
Enclosed volume simulation to test tracking system performance using pyroomacoustics. 
A large (1.2m) tristar array is simulated
'''
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import pyroomacoustics as pra
import numpy as np 
import scipy.signal as signal 
import soundfile as sf

# an chirp  20 ms long
fs = 192000
durn = 0.02
t = np.linspace(0, durn, int(fs*durn))
chirp = signal.chirp(t, 90000, t[-1], 25000)
chirp *= signal.hanning(chirp.size)
chirp *= 0.25

# crate room to match Orlova Chuka end chamber 
# code based on https://pyroomacoustics.readthedocs.io/en/pypi-release/pyroomacoustics.room.html
room_dims = [7,5,5] #m
room = pra.ShoeBox(room_dims, fs, max_order=2)

# mic positions - on the room wall
R = 1.2
theta = np.pi/3
tristar_geom = np.column_stack(([0,0,0],
                        [-R*np.sin(theta),0, -R*np.cos(theta)],
                        [R*np.sin(theta), 0, -R*np.cos(theta)],
                        [0,0,R]))
tristar_geom [-1,:] += 2.0
tristar_geom[1,:] += 0.3 # move the array a bit front of the wall
room.add_microphone_array(tristar_geom)


# sources
x = [2, 2.5,  2.2, 1.5]
y = [4.2, 4.1, 4.0, 3.95]
z = [3, 3.5, 3.7, 3 ]

source_positions = np.column_stack((x,y,z))

room.add_source(source_positions[0,:], signal=chirp)
room.add_source(source_positions[1,:], signal=chirp, delay=0.2)
room.add_source(source_positions[2,:], signal=chirp, delay=0.3)

# simulate sound propagation 
room.simulate()

channel_audio = []
for i in range(4):
    channel_audio.append(room.mic_array.signals[i,:])

plt.figure()
a1 = plt.subplot(411)
plt.specgram(channel_audio[0], Fs=fs)
for i in range(2,5):
    plt.subplot(410+i, sharex=a1, sharey=a1)
    plt.specgram(channel_audio[i-1], Fs=fs)
plt.savefig('multichannel_specgram.png')

f = plt.figure()
a0 = f.add_subplot(111, projection='3d')
a0.plot(source_positions[:,0], source_positions[:,1], source_positions[:,2],'g*')
m0 = tristar_geom[:,0]
for each in range(1,4):    
    line = np.column_stack((m0, tristar_geom[:,each]))
    x,y,z = tristar_geom[:,each]
    a0.text3D(x,y,z,'ch '+str(each+1))
    a0.plot(line[0,:], line[1,:], line[2,:], '-k')
plt.savefig('simulation_plot.png') 

# save the simulated audio 
sf.write('tristar120_roomsimulation.wav',np.column_stack(channel_audio), fs)



