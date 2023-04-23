
# %%
import pytest
import numpy as np
import gravitools.double3 as d3
import gravitools.twobody_tools as tt 
from gravitools.density_profile_generator import uniform_sphere
from gravitools.binding_energy_calculator import binding_energy_of_uniform_sphere
from gravitools.manybody_tools import kinetic_gravitational_potential

#####################################################

def test_double3_scalar_product():
    a = d3.Double3(1,2,3)
    b = d3.Double3(3,2,1)
    res = d3.scalar_product(a, b)

    assert(res==10)


def test_double3_sum():
    a = d3.Double3(1,2,3)
    b = d3.Double3(3,2,1)
    c = a + b 

    assert(c.x==4)
    assert(c.y==4)
    assert(c.z==4)


def test_double3_sub():
    a = d3.Double3(1,2,3)
    b = d3.Double3(3,2,1)
    c = a - b 

    assert(c.x==-2)
    assert(c.y==0)
    assert(c.z==2)

#####################################################

def test_star_tools_orbit():
    Mtot  = 1
    peri0 = 1.0
    ecc0  = 0.9

    t, x, y, Vx, Vy = tt.keplerian_trajectory(Mtot, peri=peri0, ecc=ecc0, theta=np.array([1.1]))

    dR = np.array([x[0], y[0], 0])
    dV = np.array([Vx[0], Vy[0], 0])

    ecc, a, peri, d = tt.ecc_semi_peri_distance(Mtot, dR, dV)

    Dperi = abs(peri0-peri)
    Decc = abs(ecc0-ecc)

    assert(Dperi<1e-5)
    assert(Decc<1e-5)

#####################################################

def test_potential_energy():

    R = 2
    N = int(3e+4)

    X, Y, Z = uniform_sphere(N=N, R=R)
    pos     = np.stack([X, Y, Z], axis=1)
    mass    = np.full(N, 1)

    _, U = kinetic_gravitational_potential(M=mass, pos=pos, vel=pos)
    U0   = binding_energy_of_uniform_sphere(R=R, M=N)

    assert( abs((U-U0)/U0)< 1e-2 )


# %%
