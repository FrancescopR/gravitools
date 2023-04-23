


# %%

import matplotlib.pyplot as plt

from gravitools.density_profile_generator import uniform_sphere
from gravitools.binding_energy_calculator import binding_energy_of_uniform_sphere
from gravitools.manybody_tools import kinetic_gravitational_potential
import numpy as np



# %% Sample coordinated from uniform sphere

R = 5
N = int(1e+4)


X,Y,Z = uniform_sphere(N=N, R=R)
pos   = np.stack([X, Y, Z], axis=1)
mass  = np.full(N, 1)


# %%
_, U = kinetic_gravitational_potential(M=mass, pos=pos, vel=pos)
U0   = binding_energy_of_uniform_sphere(R=R, M=N)



# %%
print( abs((U-U0)/U0) )

#


# %%



# %%