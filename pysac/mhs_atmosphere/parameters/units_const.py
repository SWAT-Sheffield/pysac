# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 14:02:22 2014

@author: sm1fg
"""
import numpy as np
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
                'energy':escale, 
                'time':tmscale, 
                'mass':mscale,
                'temperature':Tscale,
                'magnetic':Bscale,
                'force':Fscale
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
                'energy':escale, 
                'time':tmscale, 
                'mass':mscale,
                'temperature':Tscale,
                'magnetic':Bscale,
                'force':Fscale
               }
#============================================================================
# physical constants
#============================================================================
    mu0_SI = 4 * np.pi * 1e-7 # magnetic vacuum permeability 
    kB_SI  = 1.3806503e-23 # Boltsmann constant m^2kg s^-2 K^-1
    mp_SI  = 1.6728e-27 # proton mass kg
    g0_SI  = -274.0 #m/s/s gravitational acceleration at solar radius
    gamma0 = 5.0/3.0 # ratio of specific heats
    mu     = 0.602   #mean molecullar weight fully ionized 0.908 H 0.092 He
    mu0    = mu0_SI / Bscale**2*escale
    kB     = kB_SI/scales['mass']\
                  /scales['velocity']**2*scales['temperature'] 
    mp     = mp_SI/scales['mass'] 
    g0     = g0_SI/scales['velocity']*scales['time'] 
    if not logical_pars['l_SI']:
        if logical_pars['l_CGS']:
            mu0    = mu0_SI / (Bscale*1e-4)**2*escale
            kB     = kB_SI/(scales['mass']*1e-3)\
                          /(scales['velocity']*1e-2)**2*scales['temperature'] 
            mp     = mp_SI/(scales['mass']*1e-3) 
            g0     = g0_SI/(scales['velocity']*1e-2)*scales['time'] 
        
    physical_constants = {'gamma':     gamma0 , 
                          'mu':        mu     , 
                          'mu0':       mu0    , 
                          'boltzmann': kB     ,
                          'proton_mass':mp    ,
                          'gravity':   g0
                         }

    print"parameters for model "+model                            
    return scales, physical_constants
