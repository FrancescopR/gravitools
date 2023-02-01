
# %%
import pytest
import numpy as np
import gravitools.double3 as d3
import gravitools.twobody_tools as tt 


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


# %%
