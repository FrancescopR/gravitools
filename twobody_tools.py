"""
Created on Wed Aug 20 12:37:11 2021

@author: Francesco Paolo Rizzuto
"""


import numpy as np
import pandas as pd
import numba as nb
from units import Consts
from math import pi, fabs
import setting

####################### ORBITAL ELEMENTS ########################


@nb.jit(nopython=setting.NOPYTHON)	
def _ecc_semi_peri_distance(Mb, dR, dV, G=None):
	"""
	Compute eccentricity, semi-major axis, periapsis distance and separation distance
	from relative position and velocity vectors.

	Parameters
	----------
	Mb : float
		Total mass.
	dR : numpy.ndarray
		Relative position vector.
	dV : numpy.ndarray
		Relative velocity vector.
	G : float
		Gravitational constant. Default is None, which means it will use the value of Consts.G.

	Returns
	-------
	eccentricity : float
		Eccentricity of the orbit.
	semi_major_axis : float
		Semi-major axis of the orbit.
	periapsis_distance : float
		Periapsis distance of the orbit.
	separation_distance : float
		Separation distance of the two bodies.

	"""

	d	 = np.dot(dR, dR)**0.5

	if(fabs(d)<1e-18):
		return -1.0, 0.0, 0.0, 0.0

	VxR  = np.cross(dR, dV)  
	VxR2 = np.dot(VxR, VxR)

	# FIND SEMI-MAJOR AXES
	a_inv  = ( 2.0/d - np.dot(dV, dV) / Mb / G )	   
	if (a_inv<0.0) and (fabs(a_inv)<1e-16): # to prevent semi to be infinity
		a = -1e+16
	elif (a_inv>=0.0) and (fabs(a_inv)<1e-16): 
		a =  1e+16
	else:
		a = 1.0/a_inv	 

	# FIND ECCENTRICITY
	r0	   = VxR2 / Mb / G	# semi-latus rectum 
	ecc2   = 1.0  - a_inv * r0	
	if(ecc2<1e-16):  # to avoid negative values for ecc
		ecc = 0.0
	else:
		ecc = ecc2**0.5
		
	# FIND PERICENTER (note that peri=a * (1-ecc) would lead to 
	#				   erros when a is not well defined, around e=1.0)	  
	peri = r0 / (1.0 + ecc)

	return ecc, a, peri, d 


def ecc_semi_peri_distance(Mb, dR, dV):
    return _ecc_semi_peri_distance(Mb, dR, dV, G=Consts.G)


def true_anomaly(r, peri, ecc):
    """Calculate the true anomaly from radial distance and orbital parameters.
    
    Parameters:
    ----------
    r : float
        radial distance
    peri : float
        periapsis
    ecc : float
        eccentricity
        
    Returns:
    -------
    float
        True anomaly
    """

    r0	= peri * (1 + ecc)
    Arg = (r0/r -1) / ecc
    
    if(Arg>1.0):
        return pi
    elif(Arg<-1.0):
        return -pi 
    
    return np.arccos(Arg)


def eccentric_anomaly(theta=None, ecc=None):
        """Calculate the eccentric anomaly given the true anomaly and eccentricity.
        
        Parameters:
        ----------
        theta : float, optional
            True anomaly (default is None)
        ecc : float, optional
            Eccentricity (default is None)
            
        Returns:
        -------
        float
            Eccentric anomaly
        """
        real = np.sqrt(1+ecc) * np.cos(theta/2.0)
        imag = np.sqrt(1-ecc) * np.sin(theta/2.0)

        #z = complex(real, imag)

        E = 2 * np.arccos( real/np.sqrt( real**2 + imag**2 ) )

        return E


def hyperbolic_anomaly(theta=None, ecc=None):
        """Calculate the hyperbolic anomaly given the true anomaly and eccentricity.
        
        Parameters:
        ----------
        theta : float, optional
            True anomaly (default is None)
        ecc : float, optional
            Eccentricity (default is None)
            
        Returns:
        -------
        float
            Hyperbolic anomaly
        """

        H = 2 * np.arctanh( np.sqrt( (ecc - 1) / (ecc + 1) ) * np.tan(theta/2) )

        return H

###################################################################################


def potential_energy(m1, m2, r):
    """
    Calculates the potential energy of a two-body system.
    
    Parameters:
    m1 (float): mass of the first body
    m2 (float): mass of the second body
    r (float): distance between the two bodies
    
    Returns:
    float: potential energy
    """
    G = Consts.G
    return -G * m1 * m2 / r

 
def energy_2body(m1, m2, semi):
    """
    Calculates the energy of a two-body system.
    
    Parameters:
    m1 (float): mass of the first body
    m2 (float): mass of the second body
    semi (float): semi-major axis of the binary system
    
    Returns:
    float: energy of the two-body system
    """	
    G = Consts.G
    return -0.5*G * m1 * m2 / semi


####################### TIME RELATED VARIABLES ##################

def binary_period(Mtot=None, semi=None):
	"""
	Calculates the period of a binary system.

	Parameters:
	Mtot (float): total mass of the binary system
	semi (float): semi-major axis of the binary system

	Returns:
	float: period of the binary system
	"""
	G     = Consts.G
	try:
		P = 2.0 * pi * semi**1.5 / ( G * Mtot )**0.5		
	except:
		#raise Error("ecc > 1.0 ==> the period cannot be calculated")
		P = -10
	return P


def time_at_true_anomaly(theta=None, Mtot=None, ecc=None, peri=None):
	"""
	Calculates the time elapsed from periapsis passage given the true anomaly.

	Parameters:
	theta (float): True anomaly in radians.
	Mtot (float): Total mass of the system in kilograms.
	ecc (float): Eccentricity of the orbit.
	peri (float): Periapsis distance in meters.

	Returns:
	float: Time elapsed from periapsis passage in seconds.
	"""
	G     = Consts.G
	alpha = G * Mtot

	semi = peri / (1.0-ecc)

	if(ecc < 1.0):
		E = eccentric_anomaly(theta=theta, ecc=ecc)
		t     = semi * np.sqrt(semi/alpha) * ( E - ecc * np.sin(E))
	elif(ecc>1.0):
		semi  = np.abs(semi)
		H     = hyperbolic_anomaly(theta=theta, ecc=ecc)
		t     = semi * np.sqrt(semi/alpha) * ( ecc * np.sinh(H) - H )

	return t



def find_time(theta=None, Mtot=None, ecc=None, peri=None):
	"""
	Calculates the time elapsed from periapsis passage given the mean anomaly.

	Parameters:
	theta (float): Mean anomaly in radians.
	Mtot (float): Total mass of the system in kilograms.
	ecc (float): Eccentricity of the orbit.
	peri (float): Periapsis distance in meters.

	Returns:
	float: Time elapsed from periapsis passage in seconds.
	"""
	theta1 = theta % (2.0*pi)
	t  = time_at_true_anomaly(theta=theta1, Mtot=Mtot, ecc=ecc, peri=peri)
	t0 = time_at_true_anomaly(theta=2*pi, Mtot=Mtot, ecc=ecc, peri=peri)

	I  = (theta / (2.0*pi))
	I  = I.astype(int)

	t += t0*I.astype(int)

	return t

####################################################################

def keplerian_trajectory(Mtot, peri, ecc, theta):
    """
   
    Parameters
    ----------
    Mtot : float
        total mass of the two-body system.
    peri : float
        pericenter.
    ecc : float
        eccentricity.
    theta : numpy array/ float
        True anomaly (trajectory sampled from -th0 to +th0).
    N : integer
        Number of points sampled
   
    Returns
    -------
    x : numpy array / flaot
        x coordinate.
    y : numpy array / flaot
        y coordinate.
    Vx : numpy array / flaot
        x-velocity component.
    Vy : numpy array / flaot
        y-velocity component.
   
    """
    G = Consts.G

    theta = np.array(theta)

    if(ecc>=1.0):
        theta_max = np.abs( np.arccos(-1/ecc) )
        theta     = theta[ np.abs(theta) < theta_max ]
        # assert( np.max(np.abs(theta)) < np.abs( np.arccos(-1/ecc) ) )

    l     = peri * (1+ecc)   # semi latus-rectum
    alpha = G * Mtot
    Vx    = -(alpha / l )**0.5 * np.sin(theta)
    Vy    =  (alpha / l )**0.5 * ( ecc + np.cos(theta) )


    x = np.cos(theta) * peri * (1+ecc) / (1.0 + ecc * np.cos(theta))
    y = np.sin(theta) * peri * (1+ecc) / (1.0 + ecc * np.cos(theta))


    t = find_time(theta=theta, Mtot=Mtot, ecc=ecc, peri=peri)

    return t, x, y, Vx, Vy



def keplerian_trajectories(m1, m2, peri, ecc, theta):
	"""
	This function computes the position and velocity of two stars based on their masses, periapsis, eccentricity, and true anomaly angle.

	Parameters:
	m1 (float): Mass of star 1.
	m2 (float): Mass of star 2.
	peri (float): Periapsis of the orbit.
	ecc (float): Eccentricity of the orbit.
	theta (float): True anomaly angle.

	Returns:
	Tuple of two pandas DataFrames, representing the position and velocity of the two stars at each time.
	"""

	t, x, y, Vx, Vy = keplerian_trajectory(Mtot=m1+m2, peri=0.01, ecc=0.4, theta=theta)

	ID   = np.full(t.size, 0)
	mass = np.full(t.size, m1)
	mu   = m1*m2 / (m1+m2)

	f1  = mu/m1
	d1  = {"ID":ID, "t":t, "m":mass, 
			"x" : f1*x,  "y" : f1*y,  "z" : 0*x,
			"vx": f1*Vx, "vy": f1*Vy, "vz": 0*Vx}
	df1 = pd.DataFrame(d1)

	ID   = np.full(t.size, 1)
	mass = np.full(t.size, m2)

	f2   = mu/m2
	d2 = {"ID":ID, "t":t, "m":mass, 
			"x" : -f2*x,  "y" : -f2*y,  "z" : -0*x,
			"vx": -f2*Vx, "vy": -f2*Vy, "vz": -0*Vx}
	df2 = pd.DataFrame(d2)

	return  df1, df2
    
################## INITIALIZATION ################


def put_2stars_in_orbit(s1,s2, ecc, peri, theta=0.0):
    """
    This function updates the position and velocity of two stars based on their initial conditions (eccentricity, periapsis, and true anomaly angle).

    Parameters:
    s1 (Star object): Star 1.
    s2 (Star object): Star 2.
    ecc (float): Eccentricity of the orbit.
    peri (float): Periapsis of the orbit.
    theta (float): True anomaly angle (default: 0.0).

    Returns:
    None
    """	

    d1, d2 = keplerian_trajectories(s1.m, s2.m, ecc, peri, theta=[theta])

    x,y,z,vx,vy,vz = d1["x y z vx vy vz".split()].values[0]
    s1.update_pos_vel(x,y,z,vx,vy,vz)

    x,y,z,vx,vy,vz = d2["x y z vx vy vz".split()].values[0]
    s2.update_pos_vel(x,y,z,vx,vy,vz)


    
    