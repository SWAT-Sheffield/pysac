# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 14:02:22 2014

@author: sm1fg
"""
import numpy as np
import astropy.constants as asc
import astropy.units as u
#import pysac.mhs_atmosphere.parameters.set_up as set_up

def get_parameters(model, l_mpi, logical_pars, size):
#============================================================================
    # Dimensional units in terms of SI
#============================================================================
    if logical_pars['l_SI']:
        lscale   = 1e6                     #m
        rhoscale = 1e-6                    #kg/m3
        uscale   = 1e3                     #m/s
        escale   = rhoscale*uscale**2      #
        tmscale  = lscale/uscale           #s
        mscale   = rhoscale * lscale **3   #kg
        Tscale   = 1.0                     #K 
        Bscale   = 1e-3#np.sqrt(4*np.pi*1e-7)   #mTesla mT
        Fscale   = rhoscale*uscale/tmscale #D momentum/Dt force density balance 
    
        scales   = {
                'length':lscale, 
                'density':rhoscale,
                'velocity':uscale,
                'energy density':escale, 
                'time':tmscale, 
                'mass':mscale,
                'temperature':Tscale,
                'magnetic':Bscale,
                'force density':Fscale
               }
    if logical_pars['l_CGS']:
        lscale   = 1e8                     #cm
        rhoscale = 1e-9                    #g/cm3
        uscale   = 1e5                     #cm/s
        escale   = rhoscale*uscale**2      #
        tmscale  = lscale/uscale           #s
        mscale   = rhoscale * lscale **3   #kg
        Tscale   = 1.0                     #K 
        Bscale   = 1.0#np.sqrt(4*np.pi*1e-7)   #Gauss G
        Fscale   = rhoscale*uscale/tmscale #D momentum/Dt force density balance 
    
        scales   = {
                'length':lscale, 
                'density':rhoscale,
                'velocity':uscale,
                'energy density':escale, 
                'time':tmscale, 
                'mass':mscale,
                'temperature':Tscale,
                'magnetic':Bscale,
                'force density':Fscale
               }
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

    print"parameters for model "+model                            
    return scales, physical_constants
