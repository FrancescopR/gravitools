#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 15:08:52 2020

@author: francesco
"""

from math import pi


class Standard(object):

    R_NB  = 1    # m
    V_NB  = 1    # m/s 
    T_NB  = 1    # sec
    M_NB  = 1    # kg


    M_sun   = 1.98847e+30 # kg
    R_sun   = 6.95510e8 # met
    AU      = 1.495978707e+11   # met
    PC      = 3.08567758128e+16 # met 
    YEAR    = 3.154e+7  # sec
    Myr     = YEAR * 1e+6

    # constants
    c_light = 2.997e+8 # m/s
    G       = 6.67430e-11 # N m^2 / kg^2


class Astronomical(object):
    """
        [L] = pc, [M] = Msun, [V] = km/s
    """
    # to convert in standard units
    R_NB  = Standard.PC       # m
    V_NB  = 1000.0            # m/s 
    T_NB  = R_NB/V_NB         # sec
    M_NB  = Standard.M_sun    # kg
    ###
    G       = 0.0043009         # pc M_sun^-1 (km/s)^2
    R_sun   = Standard.R_sun / R_NB
    c_light = Standard.c_light / V_NB


class MSTAR(object):

    # to convert in standard units
    R_NB  = Standard.AU       # m
    V_NB  = 4740.571712255485 # m/s 
    T_NB  = Standard.YEAR     # sec
    M_NB  = Standard.M_sun    # kg
    ###
    G       = 4.0 * pi * pi
    R_sun   = Standard.R_sun / R_NB
    c_light = Standard.c_light / V_NB



class SolarSystem(object):
    """
        [L] = AU, [M] = Msun, [T] = 4.74 yr, [V] = km/s
    """    
    
    # to convert in standard units
    R_NB = Standard.AU                        # m
    V_NB = 1000.0                             # m/s 
    T_NB = 4.743115748256183 * Standard.YEAR  # s
    M_NB = Standard.M_sun    
    ###
    G       = 887.3785227272725
    R_sun   = Standard.R_sun / R_NB
    c_light = Standard.c_light / V_NB


class Consts(object):
    R_NB = Astronomical.R_NB                     
    V_NB = Astronomical.V_NB                             
    T_NB = Astronomical.T_NB 
    M_NB = Astronomical.M_NB    
    ###
    G       = Astronomical.G
    R_sun   = Astronomical.R_sun 
    c_light = Astronomical.c_light


def choose_units(UnitInput):
    Consts.R_NB = UnitInput.R_NB                     
    Consts.V_NB = UnitInput.V_NB                             
    Consts.T_NB = UnitInput.T_NB 
    Consts.M_NB = UnitInput.M_NB    
    ###
    Consts.G       = UnitInput.G
    Consts.R_sun   = UnitInput.R_sun 
    Consts.c_light = UnitInput.c_light




#PC_TO_AU   = 206264.806
#PC_TO_RSUN = 44365682.5
#AU_TO_RSUN = 215.090898





