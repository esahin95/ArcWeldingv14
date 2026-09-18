import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

# Time directories
base = os.path.join('postProcessing', 'interface')
timeDirs = sorted(os.listdir(base), key=float)

# Initial surface
x, y, *_ = np.loadtxt(os.path.join('postProcessing', 'interface', timeDirs[0], 'isoSurface.xy'), unpack=True)

# Analytical solution
A, C, m = 0.2, 5, 0
xa = np.linspace(-0.5, 0.5, 100)
ya = 1 - 1/16*A*C*(xa*(4*xa**2 - 3) + 4/3*m*(12*xa**2 - 1))

fig, ax = plt.subplots()
lines = ax.plot(x*1e3, y*1e3, 'b', xa+0.5, ya*A, 'k', lw=2)
ax.set_xlim([0, 1])
ax.set_ylim([0, 0.4])

# Animation
def func(dir):
    x, y, *_ = np.loadtxt(os.path.join(base, dir, 'isoSurface.xy'), unpack=True)
    lines[0].set_data(x*1e3, y*1e3)
    return lines
ani = FuncAnimation(fig, func, timeDirs, interval=100, repeat=False)
ani.save("animation.mp4")