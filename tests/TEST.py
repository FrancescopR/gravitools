
# %%
# 3) Wrap the modules in a single package 
#  - pip install in anaconda env 
#  -
# 4) Add tests

from units import choose_units, MSTAR, Astronomical, Standard
choose_units(Astronomical)

# %%
import pytest


# %%
import twobody_tools 
from star_tools import Star


S1 = Star.MainSequence(mass=1.1, ID=0, z=0.001)
S2 = Star.MainSequence(mass=1.1, ID=0, z=0.02)

ecc  = 0.5
peri = 1.0

print(S1.pos, S1.vel)

print()
twobody_tools.put_2stars_in_orbit(S1,S2, ecc, peri, theta=0.0)
print()

print(S1.pos, S1.vel)


# %%
import numpy as np
print(np.array(4))
print(np.array([4]))
print(np.array(np.array([4])))

# %%
import setting

setting.PARALLEL    = True
setting.NUM_THREADS = 4

import manybody_tools as mt
import numpy as np
import time


# %%

np.random.seed(23)
N = 70000

M   = np.random.rand(N)
pos = np.random.rand(N, 3)
vel = np.random.rand(N, 3)

start = time.time()
K, U  = mt.kinetic_gravitational_potential(M, pos, vel)
end   = time.time()

print(end-start)

# %%

# %% TESTINGS
# 3) Test ecc semi and peri in twobody_tools
import twobody_tools as tt

Mtot = 1

t, x, y, Vx, Vy = tt.keplerian_trajectory(Mtot, peri=1.0, ecc=0.9, theta=np.array([1.1]))

dR = np.array([x[0], y[0], 0])
dV = np.array([Vx[0], Vy[0], 0])

ecc, a, peri, d = tt.ecc_semi_peri_distance(Mtot, dR, dV)

print(ecc, a, peri, d)

