
from units import Consts
import double3
import numpy as np

class Star(object):
    def __init__(self, mass, kw, radius, ID=None):
        self.id   = ID
        self.m    = mass
        self.r    = radius
        self.type = kw
        self.pos  = double3.Double3()
        self.vel  = double3.Double3()
        self.acc  = double3.Double3()

    @classmethod
    def MainSequence(cls, *, mass=None, ID=None, z=None):
        radius = rzamsf_Tout1996(mass, z) # Warning this work only if the mass is in Msun 
        if(mass<0.7):
            kw = 0
        elif(mass>=0.7):
            kw=1
        return cls(mass, kw, radius, ID)

    @classmethod
    def BlackHole(cls, *, mass=None, ID=None):
        G       = Consts.G
        c_light = Consts.c_light
        R_sun   = Consts.R_sun

        kw = 14
        radius = 2 * G * mass / c_light**2 
        return cls(mass, kw, radius, ID)

    def update_pos_vel(self, x, y, z, vx, vy, vz):
        self.pos.set_values(x=x, y=y, z=z)
        self.vel.set_values(x=vx, y=vy, z=vz)

    def __repr__(self):
        string1 = 'id={0} type={1} mass={2:.3e} radius={3:.3e}\n'.format(self.id, self.type, self.m, self.r) 
        string2 = 'pos = {}\n'.format(self.pos)
        string3 = 'vel = {}\n'.format(self.vel)
        return string1 + string2 + string3

    def as_string(self):
        ID = self.id
        m  = self.m
        r  = self.r
        kw = self.type
        pos_str = self.pos.as_string()
        vel_str = self.vel.as_string()

        return f"{ID} {m:16e} {r:16e} " + pos_str + " " + pos_str + f" {kw}"

        



def Star_binding_energy(M, R, polytropic_index):
    """Estimate gravitational binding energy of a polytropic object 

    Args:
        M (double): Stellar Mass
        R (double): Stellar Radius
        polytropic_index (double)

    Returns:
        double: binding energy
    """
    G = Consts.G
    return (3.0  * G * M * M ) / R / (5.0 - polytropic_index)



def angular_momentum(star):
	L =  star.m * double3.cross_product(star.pos,  star.vel)
	return L 


def kinetic_energy(star):
	T = 0.5 * star.m * double3.scalar_product(star.vel,  star.vel)
	return T 






# fitting parameters
xz = [666, 0.39704170, -0.32913574,  0.34776688,  0.37470851, 0.09011915, 8.52762600,-24.41225973, 56.43597107, 37.06152575, 5.45624060,	0.00025546, -0.00123461, -0.00023246,  0.00045519, 
	  0.00016176,	 5.43288900, -8.62157806, 13.44202049, 14.51584135, 3.39793084,    5.56357900,-10.32345224, 19.44322980, 18.97361347, 4.16903097,	 0.78866060, -2.90870942,  6.54713531, 
	  4.05606657, 0.53287322,	 0.00586685, -0.01704237,  0.03872348,	0.02570041, 0.00383376,    1.71535900,	0.62246212, -0.92557761, -1.16996966,-0.30631491,	 6.59778800, -0.42450044,
	  -12.13339427,-10.73509484,-2.51487077,   10.08855000, -7.11727086,-31.67119479,-24.24848322,-5.33608972,	  1.01249500,  0.32699690, -0.00923418, -0.03876858,-0.00412750,	0.07490166, 
	  0.02410413,  0.07233664,	0.03040467, 0.00197741,    0.01077422,	  3.08223400,  0.94472050, -2.15200882, -2.49219496,-0.63848738,   17.84778000, -7.45345690,-48.96066856,
	  -40.05386135,-9.09331816,   0.00022582, -0.00186899,	0.00388783,  0.00142402,-0.00007671]


def fitting_parameters_Tout1996(z):
	lzs = np.log10(z) - np.log10(0.02)
	return {8 : xz[36]+lzs*(xz[37]+lzs*(xz[38]+lzs*(xz[39]+lzs*xz[40]))),
			9 : xz[41]+lzs*(xz[42]+lzs*(xz[43]+lzs*(xz[44]+lzs*xz[45]))),
			10:  xz[46]+lzs*(xz[47]+lzs*(xz[48]+lzs*(xz[49]+lzs*xz[50]))),
			11:  xz[51]+lzs*(xz[52]+lzs*(xz[53]+lzs*(xz[54]+lzs*xz[55]))),
			12:  xz[56]+lzs*(xz[57]+lzs*(xz[58]+lzs*(xz[59]+lzs*xz[60]))),
			13:  xz[61],
			14:  xz[62]+lzs*(xz[63]+lzs*(xz[64]+lzs*(xz[65]+lzs*xz[66]))),
			15:  xz[67]+lzs*(xz[68]+lzs*(xz[69]+lzs*(xz[70]+lzs*xz[71]))),
		    16:  xz[72]+lzs*(xz[73]+lzs*(xz[74]+lzs*(xz[75]+lzs*xz[76])))}
    
    

def rzamsf_Tout1996(m, z):
	"""A function to estimate the radius of a zero-age main sequence star
		( from Tout et al., 1996, MNRAS, 281, 257 ).

	Args:
		m (double): Mass of the star in [MSUN]
		z (double): stellar metallicity.

	Returns:
		double: zero-age main sequence star in [RSUN] 
	"""
	main_sequence_parameters = fitting_parameters_Tout1996(z)
	m2	 = m*m
	m05  = m**0.5
	m25  = m05*m2
	m65  = m25*m2*m2
	m90  = m25*m65
	m190 = m90*m90*m
	m195 = m190*m05
	A = (main_sequence_parameters[8]*m25 +	main_sequence_parameters[9]*m65 +  main_sequence_parameters[10]*m2*m90 +
		     main_sequence_parameters[11]*m190 + main_sequence_parameters[12]*m195)
	B = (main_sequence_parameters[13] + m2*(main_sequence_parameters[14] 
		     +	m65*(main_sequence_parameters[15] + m90*m)) + main_sequence_parameters[16]*m195)

	rzams = A/ B

	return rzams

