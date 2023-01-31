"""
Created on Wed Aug 20 11:17:28 2020

@author: Francesco Paolo Rizzuto
"""


import numpy as np
import numba as nb
import math
from units import Consts
from scipy.spatial.ckdtree	import cKDTree	


# %% ########## ANALYTICAL DENSITY PROFILES ##########

class Plummer(object):
    
    def __init__(self, M0, Rh):
        self.M0   = M0
        self.Rh   = Rh

        a = Rh * ( 0.5**(-2/3) - 1 )**0.5

        G = Consts.G

        self.Emax = G * M0 / a
        self.a    = a    

    def density_profile(self, r):

        M0 = self.M0
        a  = self.a
        A  = 3 * M0 / 4 / math.pi / a**3

        return A * (1 + (r/a)**2 )**(-5/2)  

    def potential_energy(self, r):
        M0 = self.M0
        a  = self.a
        G  = Consts.G

        return -G * M0 / np.sqrt(r*r + a*a) 

    def sigma(self, r):
        M0 = self.M0
        a  = self.a
        G  = Consts.G

        return math.sqrt(0.5 * G * M0 / np.sqrt(r*r + a*a)) 
        
    def energy(self, r):
        M0 = self.M0
        a  = self.a
        G  = Consts.G

        return 0.75 * G * M0 / np.sqrt(r*r + a*a) 




class StoneOstriker(object):
    
    def __init__(self, M0, Rh, rc):

        self.M0 = M0  
        self.Rh = Rh
        self.rc = rc

        # derived properties
        self.rho_c = self._core_density()
        self.rhalo = self._outer_halo_radius()

    def _outer_halo_radius(self, ):
        
        Rh = self.Rh
        rc = self.rc
        C  = rc * (math.pi -  math.atan(Rh/rc)) / Rh
     
        x = Rh / Rh

        return (math.pi - math.atan(x)) / x - C

    def _core_density(self):

        M0 = self.M0
        Rh = self.Rh
        rc = self.rc

        Den = 2*(math.pi*rc*Rh)**2 

        return M0 * (Rh+rc) / Den   


    def density_profile(self, r):

        rho_c = self.rho_c
        Rh    = self.Rh
        rc    = self.rc

        Den = (1.0 + (r/rc)**2 ) * (1.0 + (r/Rh)**2 ) 
        
        return rho_c / Den

    def mass_profile(self, r):
        M0 = self.M0
        rh = self.rhalo
        rc = self.rc

        return 2 * M0 * (rh * np.arctan(r/rh) - rc * np.arctan(r/rc)) / math.pi / (rh-rc)


    def potential_energy(self, r):

        rho_c = self.rho_c
        Rh    = self.Rh
        rc    = self.rc

        G = Consts.G

        A = 4.0 * math.pi * G * rc*rc * Rh*Rh * rho_c / ( Rh*Rh - rc*rc  )
        
        B  = Rh * math.atan(r/Rh) / r  
        B -= rc * math.atan(r/rc) / r
        B += 0.5 * math.log( (r*r + Rh*Rh) / (r*r + rc*rc) )    

        return A * B 

    def central_escape_velocity(self):

        M0 = self.M0
        Rh = self.Rh
        rc = self.rc

        G = Consts.G

        V = 2*math.sqrt( G * M0 * math.log( Rh / rc ) / (Rh - rc) / math.pi )

        return V 

    def core_dispersion_velocity(self):

        rho_c = self.rho_c
        Rh    = self.Rh
        rc    = self.rc

        G = Consts.G

        sig2 = 12 * math.pi * G * rho_c * rc*rc * (math.pi*math.pi/8 - 1)

        return np.sqrt(sig2)


# %% ########## EXTRACT STAR CLUSTERS PROPERTIES ##########


def lagrangian_radius(R, M, frac=0.5):
	
	R = np.array(R)
	M = np.array(M)
	
	N	  = M.shape[0]
	M_tot = M.sum()
	
	# sort masses according to their distance 
	indx_sort = np.argsort(R) 
	R = R[indx_sort]
	M = M[indx_sort]

	
	for i in range(0, N, 100):
		Mshell = M[:i].sum()
	
		if Mshell/M_tot > frac:
			if((Mshell/M_tot) - frac > 0.01):
				print("Warning Mshell/M_tot=", Mshell/M_tot)
			
			return R[i] 
 
 
 
def dispersion_velocity(Vx, Vy, Vz):
   
	sigx = np.std(Vx)
	sigy = np.std(Vy)
	sigz = np.std(Vz)
 
	sigma = (sigx**2 + sigy**2 + sigz**2)**0.5
 
	return sigma 
 
 
 
 
def dispersion_velocity_profile(R, Vx, Vy, Vz, step=1000):
	
    i_sort = np.argsort(R)
    R = R[i_sort]
    Vx = Vx[i_sort]
    Vy = Vy[i_sort]
    Vz = Vz[i_sort]

    Sig = np.array([]) # dispersion velocity
    d_c = np.array([]) # central distance
		
    for i in range(0, R.size(), step):
        Sig = np.append(Sig, dispersion_velocity(Vx[i:i+step],
                                                 Vy[i:i+step], Vz[i:i+step] ))
        d_c = np.append(d_c, R[i:i+step].values.mean())		  
		
    return d_c, Sig

 
 
 
def mass_and_particle_density(R, M, step = 1000):

    i_sort = np.argsort()

    R = R[i_sort]    
    M = M[i_sort]

    r, Rho, Rho_n = [], [], []	   
    for i in range(0, R.size() - step, step):
        DM = M[i:i+step].sum()
        DN = M[i:i+step].shape[0]
        # Compute the density
        volume = 4.0 * math.pi * (R[i+step-1]**3 - R[i]**3 ) / 3.0
        Rho    += [DM / volume]
        Rho_n  += [DN / volume]
        r	   += [0.5 * (R[i] + R[i+step-1])]
		#print(i, i+step, DM, r)
    return np.array(r), np.array(Rho), np.array(Rho_n) 



def relaxation_time(N=None, m_mean=None, Rh=None):

    G = Consts.G
    A = 0.138 * (N)**(0.5)
    B = (Rh**3 / G / m_mean)**0.5
    Trh = A * B / np.log(0.1*N)
    
    print("Half Mass Relax. Time = {0:.3e} Myr".format(Trh))
    
    return Trh


def sample_array(x, fraction=0.5, use_values_as_prob=False):
	size	= int(fraction * len(x))

	if(use_values_as_prob):
		indices = np.random.choice(len(x), size=size, p=x/x.sum(), replace=False)
	else:
		indices = np.random.choice(len(x), size=size, p=x/x.sum(), replace=False)
	
	return indices




def density_center(pos, vel, mass, sampling=1.0):
	# take a sample of the original arrays
	i_sample = sample_array(mass, fraction=sampling, use_values_as_prob=True)	 
	pos_	 = pos[i_sample]
	vel_	 = vel[i_sample]
	mass_	 = mass[i_sample]

	# find 5 closest neighborhood particles
	tree = cKDTree(pos_)
	r, indices = tree.query((pos_), k=6)
	
	# find distance of the last neighborhood
	r6	   = r[:,-1]
	
	# find most massive
	m6_tot = mass_[indices].sum() 
	
	# compute local density for each particle
	Rho = m6_tot / r6**3
	
	# find density center
	Rho_tot = Rho.sum()
	
	xc = np.sum(pos_[:,0] * Rho / Rho_tot)
	yc = np.sum(pos_[:,1] * Rho / Rho_tot)
	zc = np.sum(pos_[:,2] * Rho / Rho_tot) 
	
	# the center of velocity must hcange as well
	Vxc = np.sum(vel_[:,0] * Rho / Rho_tot)
	Vyc = np.sum(vel_[:,1] * Rho / Rho_tot)
	Vzc = np.sum(vel_[:,2] * Rho / Rho_tot) 

	# padding Rho to the original size
	Rho_		   = np.zeros(len(mass))
	Rho_[i_sample] = Rho 
	
	return Rho_, [xc, yc, zc], [Vxc, Vyc, Vzc] 



def core_radius_mass_density(R, M, Rho): 

	rc = np.sqrt( (R**2 * Rho**2).sum() / (Rho**2).sum() )

	Mc = M[R < rc].sum()
	Nc = sum((R < rc))
	Rho_c = 3.0 * Mc / 4.0 / math.pi / rc**3		
	n_c   = 3.0 * Nc / 4.0 / math.pi / rc**3		

	return rc, Mc, Rho_c, n_c



