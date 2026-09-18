#!./.venv/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# Base mesh
L, Nx = 1e-3, 500
x = np.linspace(0, L, Nx)

# Time directories
base = os.path.join('postProcessing', 'lineUniform')
timeDirs = sorted(os.listdir(base), key=float)

# Precomputed eigenfunctions
D = 1e-5
M = 20
k = np.arange(M)
w = (2*k+1)*np.pi/L
C = np.cos(np.outer(w, x))
coef = 2./(L*w) * (-1)**k

# Initial plot
fig, ax = plt.subplots()
c0 = np.ones_like(x)
c0[x<0.5*L] = 0
ax.plot(x, c0, 'k')
lines = ax.plot(x, c0, 'bo', x, c0, 'r', lw=2)

# Animation
def func(dir):
    x, y = np.loadtxt(os.path.join(base, dir, 'center.xy'), unpack=True)
    lines[0].set_data(x, y)

    t = float(dir)
    lines[1].set_ydata(0.5 - (coef*np.exp(-w**2*D*t))@C)

    return lines
ani = FuncAnimation(fig, func, timeDirs, interval=200)
ani.save("animation.mp4")