# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:37:39 2014

@author: sm1fg

    Construct the magnetic network and generate the adjustments to the
    non-magnetic atmosphere for mhs equilibrium.

"""

import numpy as np
import astropy.units as u

#============================================================================
# Derive the hydrostatic profiles and include the magneto adjustments
#============================================================================
def mhs_3D_profile(z,
                   pressure_z,
                   rho_z,
                   pressure_m,
                   rho_m
                  ):
    """Return the vertical profiles for thermal pressure and density in 3D"""
    #Make things 3D
    rho_0 = np.empty(z.shape) * u.kg /u.m**3
    rho_0[:] = rho_z
    #hydrostatic vertical profile
    pressure_0 = np.empty(z.shape) * u.Pa
    pressure_0[:] = pressure_z
    #magnetohydrostatic adjusted full 3D profiles
    pressure =  pressure_m + pressure_0
    rho = rho_m+rho_0

    return pressure, rho

#============================================================================
# Calculate internal energy
#============================================================================
def get_internal_energy(pressure, magp, physical_constants):
    """ Convert pressures to internal energy -- this may need revision if an
    alternative equation of state is adopted.
    """
    energy = pressure/(physical_constants['gamma']-1) + magp
    return energy

