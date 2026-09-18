import numpy as np
import matplotlib.pyplot as plt

import surf2stl as stl

# load data for layer height
M = np.loadtxt('postProcessing/layerHeight/0/traceSurface.dat', skiprows=2)
avg, std = np.mean(M[:, 2]), np.std(M[:, 2])
print(f'Mean height {avg:.5g}, std {std:.5g}')

# export surface to a stl format file
sz = (30, 20)
X = M[:, 0].reshape(*sz)
Y = M[:, 1].reshape(*sz)
Z = M[:, 2].reshape(*sz)
stl.write('surf.stl', X, Y, Z)

# plot surface
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_trisurf(M[:, 0], M[:, 1], M[:, 2], linewidth=0.2, antialiased=True)
ax.set_aspect('equal', 'box')
ax.set_axis_off()
fig.tight_layout()
plt.show()

# load data for fill fraction
M = np.loadtxt('postProcessing/fillFraction/0/traceSurface.dat', skiprows=2)
avg, std = np.mean(M[:, 2]), np.std(M[:, 2])
print(f'Mean fraction {avg:.5g}, std {std:.5g}')

# plot surface
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_trisurf(M[:, 0], M[:, 1], M[:, 2], linewidth=0.2, antialiased=True)
ax.set_axis_off()
fig.tight_layout()
plt.show()