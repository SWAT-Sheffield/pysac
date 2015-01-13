# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 14:02:22 2014

@author: sm1fg
"""
import numpy as np
import astropy.constants as asc
import astropy.units as u
#import pysac.mhs_atmosphere.parameters.set_up as set_up

def get_parameters():
#============================================================================
    # Dimensional units in terms of SI
#============================================================================
    scales   = {
            'length':         1*u.Mm,
            'density':        1e-6*u.kg/u.m**3,
            'velocity':       1e3*u.m/u.s,
            'temperature':    1.0*u.K, 
            'magnetic':       1e-3*u.mT
           }
    scales['energy density'] = scales['density'] * scales['velocity']**2
    scales['time'] = scales['length'] / scales['velocity'] 
    scales['mass'] = scales['density'] * scales['length']**3 
    scales['force density'] = scales['density'] * scales['velocity'] / \
                       scales['time'] #D momentum/Dt force density balance 
#============================================================================
# physical constants
#============================================================================
    physical_constants = {'gamma':      5.0/3.0         , 
                          'mu':         0.602           , 
                          'mu0':        asc.mu0         , 
                          'boltzmann':  asc.k_B         ,
                          'proton_mass':asc.m_p         ,
                          'gravity':    -274.0*u.km/u.s/u.s
                         }

    return scales, physical_constants
