import math
from gravitools.units import Consts

# %%


def binding_energy_of_uniform_sphere(R=1, M=1):

    G = Consts.G

    return -G * 3*M*M/5/R



def density_of_uniform_sphere(r, R=1, M=1):

    if(r>R):
        return 0
    else:
        Vol = 4 * math.pi * R**3 / 3.0
        return M / Vol

    return None



def binding_energy(rho=density_of_uniform_sphere, R=1.0, M=1.0):
    """
    Calculate the binding energy of a spherical body with a given
    3D spherically symmetric density profile.

    Parameters:
    rho (function): A function representing the spherically symmetric density profile of the body as a function of
                    radial distance (r), total mass (M), and radius (R).
    R (float, optional): Radius of the spherical body (default is 1.0).
    M (float, optional): Total mass of the spherical body (default is 1.0).

    Returns:
    float: The binding energy of the spherical body.

    """  
    G = Consts.G

    dr = 0.0001
    r  = 0.0001
    U  = 0.0
    while(r < R):

        dm_ext = 4 * math.pi * rho(r=r, M=M, R=R) * r*r * dr
        m_int  = 4 * math.pi * r**3 * rho(r=r, M=M, R=R) / 3.0 

        dU = - G * dm_ext * m_int / r

        U += dU 
        r += dr

    return U

