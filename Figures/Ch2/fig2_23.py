# Beam pattern for 10-element uniform array (d=lambda/2) scanned to 30 degrees (60 degrees from broadside)

import numpy as np
import matplotlib.pyplot as plt

N = 10
n = np.array([np.arange((N-1)/-2, (N)/2)]).T # 10x1
theta = np.array([np.pi * np.arange(-1, 1, 0.001)]) # 1x2000
u = np.cos(theta)
d = 0.5
vv = np.exp(2j * np.pi * d * n @ u) # 10x2000
theta_T = 30 / 180 * np.pi
w = 1 / N * np.ones((N,1)) * np.exp(-1j * np.pi * n * np.cos(theta_T)) # 10x1  NOTE I HAD TO ADD THE -1 TO GET RESULTS TO MATCH MATLAB
B = w.T @ vv # 1x2000
B = B.squeeze() # make it a 1D array
B = 10 * np.log10(np.abs(B)**2)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(np.squeeze(theta),B) # MAKE SURE TO USE RADIAN FOR POLAR
ax.set_theta_zero_location('N') # make 0 degrees point up
ax.set_theta_direction(-1) # increase clockwise
ax.set_rgrids([0,-10,-20,-30,-40])
ax.set_ylim([-40,0])
ax.set_xticks(np.arange(0, 2*np.pi, np.pi/6)) # ticks every 30 degrees
ax.set_rlabel_position(75)  # Move grid labels away from other labels, in degrees
plt.show()

