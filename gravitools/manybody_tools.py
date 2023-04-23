import numba as nb
import numpy as np
from gravitools.units import Consts
from gravitools import twobody_tools
from gravitools import setting

#nb.set_num_threads(nb.config.NUMBA_NUM_THREADS)
nb.set_num_threads(setting.NUM_THREADS)


@nb.jit(nopython=setting.NOPYTHON, parallel=setting.PARALLEL)
def _kinetic_gravitational_potential(M, pos, vel):
	'''
	
	Parameters
	----------
	M : 1D float arrays
		masses.
	pos : 3D float arrays
		positions of the N particles (3xN matrix).
	vel : 3D float arrays
		velocities of the N particles (3xN matrix)..

	Returns
	-------
	K : float
		Total Kinetic Energy.
	U : float
		Total Potential Energy.

	'''
	N = M.shape[0]
	U = 0.0
	K = 0.0; 
	for i in nb.prange(N):
		K  += 0.5 * M[i] * np.dot(vel[i], vel[i])
		for j in nb.prange(N):
			if(j>i):
				r  = np.sqrt( ( pos[i, 0] - pos[j, 0] )**2 +  
							  ( pos[i, 1] - pos[j, 1] )**2 +  
							  ( pos[i, 2] - pos[j, 2] )**2 
							)  
				U -= M[j] * M[i] / r
	return K, U

def kinetic_gravitational_potential(M, pos, vel):
    G = Consts.G
    K, U = _kinetic_gravitational_potential(M, pos, vel)
    return K, G*U


@nb.jit(nopython=setting.NOPYTHON)
def ecc_semi_peri_distance_many_(Mb, dR, dV, ecc, semi, peri):
    
    for i in nb.prange(Mb.size):
        ecc[i], semi[i], peri[i], _ = twobody_tools.ecc_semi_peri_distance(Mb[i], dR[i], dV[i])        


@nb.jit(nopython=setting.NOPYTHON)
def ecc_semi_peri_distance_many(Mb, dR, dV):
    
    ecc  = np.empty(Mb.size, dtype=np.float64)
    peri = np.empty(Mb.size, dtype=np.float64)
    semi = np.empty(Mb.size, dtype=np.float64)
    
    ecc_semi_peri_distance_many_(Mb, dR, dV, ecc, semi, peri)
    
    return ecc, semi, peri
    


