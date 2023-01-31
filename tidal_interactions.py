
"""
Created on Wed Aug 26 13:47:58 2020

@author: Francesco Paolo Rizzuto
"""



import numpy as np
from math import pi
import numba as nb
import twobody_tools
import setting
from units import Consts


@nb.jit(nopython=setting.NOPYTHON)
def tidal_capture_efficiency_factor2(polytropic_index=None, zeta=None):
	"""
	Compute the tidal energy factor2 using a fifth-order polynomial fit

	Parameters:
	polytropic_index (float): polytropic index of the star
	zeta (float): zeta value for the star

	Returns:
	float: tidal energy factor2
	"""

	c = np.zeros(6)
	eff_2 = 0.0

	if(polytropic_index==1.5):
		c[0] =-0.397
		c[1] = 1.678
		c[2] = 1.277
		c[3] =-12.42
		c[4] = 9.446
		c[5] =-5.550
	elif(polytropic_index==2.0):
		c[0] =-0.517
		c[1] =-0.906
		c[2] = 23.88
		c[3] =-93.49
		c[4] = 112.3
		c[5] =-44.15
	elif(polytropic_index==3.0):
		c[0] =-1.124
		c[1] = 0.877
		c[2] =-13.37
		c[3] = 21.55
		c[4] =-16.48
		c[5] = 4.124
	else:
		print("Incorrect polytropic index:", polytropic_index)
		#raise Exception("Incorrect polytropic index:", polytropic_index)


	# the tidal energy factor computed using fifth-order polynomial fit.
	x = np.log10(zeta)

	eff_2 = ( ( ( (c[5]*x + c[4])*x + c[3])*x + c[2])*x + c[1])*x + c[0]


	try:
		if eff_2 > 20:
			print("WARNING eff_3 IS TOO LARGE", eff_2)
	except:
		pass

	return 10.0**eff_2


@nb.jit(nopython=setting.NOPYTHON)
def tidal_capture_efficiency_factor3(polytropic_index=None, zeta=None):
	"""
	Compute the tidal energy factor3 using a fifth-order polynomial fit

	Parameters:
	polytropic_index (float): polytropic index of the star
	zeta (float): zeta value for the star

	Returns:
	float: tidal energy factor3
	"""	
	c = np.zeros(6)
	eff_3 = 0.0

	if(polytropic_index==1.5):
		c[0] =-0.909
		c[1] = 1.574
		c[2] = 12.37
		c[3] =-57.40
		c[4] = 80.10
		c[5] =-46.43
	elif(polytropic_index==2.0):
		c[0] =-1.040
		c[1] =-1.354
		c[2] = 37.64
		c[3] =-139.9
		c[4] = 168.2
		c[5] =-66.53
	elif(polytropic_index==3.0):
		c[0] =-1.703
		c[1] = 2.653
		c[2] =-14.34
		c[3] = 12.85
		c[4] =-0.492
		c[5] =-3.600
	else:
		print("Incorrect polytropic index:", polytropic_index)
		#raise Exception("Incorrect polytropic index:", polytropic_index)

	#		the tidal energy factor computed using fifth-order polynomial fit.
	x = np.log10(zeta)

	eff_3 = ( ( ( (c[5]*x + c[4])*x + c[3])*x + c[2])*x + c[1])*x + c[0]
		
	if eff_3 > 20:
		print("WARNING eff_3 IS TOO LARGE", eff_3)

	return 10.0**eff_3


@nb.jit(nopython=setting.NOPYTHON)
def tidal_radius(radius, m1, m2):
	# Kochanek (1992): tidal distruption criteria with damping	
	Rt = radius * 1.3 * ( (m1+m2) / (2*m1) )**(1/3)
	return Rt


@nb.jit(nopython=setting.NOPYTHON)
def Zeta(eta=None, ecc=None):
	alpha = 0.5 + 0.25 * np.abs(0.5 * eta - 1.0)**1.5
	y = eta * (2.0/(ecc + 1.0))**alpha

	return max(1.0, y)



@nb.jit(nopython=setting.NOPYTHON) 
def tidal_energy_loss(peri=None, ecc=None, m=None, radius=None, kw=None, 
					  m_perturber=None):
  

	Mb	   = m + m_perturber
	
	# No tidal oscilations for BHs
	if kw == 14:
		dE = 0.0
		return dE

	#Determine polytropic index
	if(kw == 0):
		polytropic_index = 1.5
	elif(kw > 9):
		polytropic_index = 1.5
	else:
		polytropic_index = 3.0
	
	eta = (peri/radius)**1.5 * (m/Mb)**0.5	 

	# Mardling 2001
	zeta = Zeta(eta=eta, ecc=ecc)
	
	if zeta < 1.0:
		# note zeta cannot be lt 1.0: the fitting fnction are valid for zeta>1.0
		print("WARNING IN TIDAL INTERACTION zeta < 1.0, zeta", zeta)
		zeta=1.0

	psi2 = tidal_capture_efficiency_factor2(polytropic_index=polytropic_index, zeta=zeta)
	psi3 = tidal_capture_efficiency_factor3(polytropic_index=polytropic_index, zeta=zeta)

	R2	= (radius/peri)**2
	R6	= R2*R2*R2 * (Mb - m)**2 / radius

	G = Consts.G	
	dE2 = G * R6 * psi2			   # l = 2
	dE3 = G * R6 * psi3 * R2	   # l = 3
	dE	= dE2 + dE3
	
	return dE



@nb.jit(nopython=setting.NOPYTHON)
def drag_force_coefficient(m1, m2, R1, kw1, dr, dv, G=0.00430125637339471887):

    Mtot = m1 + m2    

    ecc, a, peri, d = twobody_tools.ecc_semi_peri_distance(Mtot, dr, dv)    

    f = 0.5 * pi * (2.0 + 7.0 * ecc**2 + ecc**4 )

    dE = tidal_energy_loss(peri, ecc, m=m1, radius=R1, kw=kw1, m_perturber=m2)
    C  = 0.5 * dE * ( peri * (1+ecc) )**3.5 / (Mtot*G)**0.5 / f
	
    return C



