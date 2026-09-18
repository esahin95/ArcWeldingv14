#!./.venv/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import erf, erfc
from scipy.optimize import root
import sys

# Numerical solution
t, x = np.loadtxt(os.path.join('postProcessing', 'data_0.csv'), unpack=True)

# Analytical solution
cp = 123
L = 44600
Ts = 544
Tl = 545
Tm = 0.5*(Ts+Tl)

TL = 552.5
TS = 540.5

k = 10.35
rho = 9780.0
a = k / rho / cp

class F:
    def __init__(self, nu, StL, StS):
        self.nu = nu
        self.StL = StL
        self.StS = StS
    def __call__(self, x):
        return (
            self.StL/(np.exp(x**2)*erf(x))
          - self.StS/(self.nu*np.exp(self.nu**2*x**2)*erfc(self.nu*x))
          - x*np.sqrt(np.pi)
        )
lam = root(
    F(1, cp*(TL-Tm)/L, cp*(Tm-TS)/L),
    1.0,
).x

ta = np.linspace(0, np.max(t), 100)
xa = 2*lam*np.sqrt(a*ta)

# Plot data
fig, ax = plt.subplots()
ax.plot(ta, xa, 'k', t, x, 'b', lw=2)
plt.savefig(
    sys.argv[1],
    format='pdf',
    bbox_inches='tight',
    pad_inches=0.0,
    transparent=True
)