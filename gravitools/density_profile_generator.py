
import numpy as np
from gravitools.units import Consts
import numba as nb
import math
import random
from gravitools import setting


nb.set_num_threads(setting.NUM_THREADS)



@nb.jit(nopython=setting.NOPYTHON, parallel=setting.PARALLEL)
def uniform_sphere(N=100, R=1):
    """Sample X, Y, Z coordinates form uniform sphere

    Args:
        N (int, optional): Number of samples. Defaults to 100.
        R (int, optional): Size of the sphere. Defaults to 1.

    Returns:
        X (array of int): x-coordinates of the N sampled particles.
        Y (array of int): y-coordinates of the N sampled particles.
        Z (array of int): z-coordinates of the N sampled particles.
    """

    X = np.empty(N)
    Y = np.empty(N)
    Z = np.empty(N)

    for i in range(N):

        r = R * random.random()**(1/3)       
  
        Z[i] = 2 * r * ( random.random() - 0.5 )

        theta = 2.0 * math.pi * random.random()

        X[i] = math.sqrt( r**2 - Z[i]**2 ) * np.cos(theta)  
        Y[i] = math.sqrt( r**2 - Z[i]**2 ) * np.sin(theta)

    return X, Y, Z

