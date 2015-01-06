# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 11:37:39 2014

@author: sm1fg

"""
##============================================================================
##logical parameters
##============================================================================
def get_logical(model, l_mpi, l_SI=True, l_gdf=True):

    """This module assigns the logical options for the model. If adding 
    new models with additional logical arguments add it to the default 
    list as false, include an if statement for True update the
    dictionary logical_pars 
    
    """    
    #default arguments
    l_CGS         = False
    if not l_SI:
        l_CGS     = True
    if l_gdf:
        suffix = '.gdf'
    else:
        suffix = ''
    l_hdonly      = False  # set magnetic field to zero to check background
    l_ambB        = False  # include some ambient magnetic field b_z
    l_spruit      = False  # thin flux tube model to check Spruit
    l_const       = False  # axial Alfven speed const  dependence on Z (Spruit)
    l_sqrt        = False  # axial Alfven speed sqrt   dependence on Z (Spruit)
    l_linear      = False  # axial Alfven speed linear dependence on Z (Spruit)
    l_square      = False  # axial Alfven speed square dependence on Z (Spruit)
    l_B0_expz     = False  # vertical strength of Bz(r=0) follows exponentials
    l_B0_quadz    = False  # vertical strength of Bz(r=0) follows polynomials 
                           # + coronal exponential    
    l_single      = False  # only one flux tube
    l_hmi         = False  # construct photopheric map of Bz from HMI/SDI
    l_tube_pair   = False # pair of flux tubes
    l_multi_bps   = False # multiple flux tubes as described in GFE (2014)
    l_multi_lanes = False # multiple flux tubes as described in GFE (2014)
    l_multi_twist = False # multiple flux tubes as described in GFE (2014)
    l_2D_loop     = False # make a 2D loop with sinusoidal Bz(x,0,0)
    l_mfe         = False # model Viktors model from MFE (2014)

    l_atmos_val3c_mtw = False # interpolate composite VAL3c+MTW atmosphere

    if model == 'mfe_setup':
        l_single      = True # only one flux tube
        l_mfe         = True # model Viktors model from MFE (2014)
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'spruit':    
        l_single      = True # only one flux tube
        l_spruit      = True  # thin flux tube model to check Spruit
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
    if model == 'paper1':
        l_ambB        = True # include some ambient magnetic field b_z
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
        l_single      = True # only one flux tube
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'paper2a':
        l_ambB        = True # include some ambient magnetic field b_z
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
        l_tube_pair   = True # pair of flux tubes
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'paper2b':
        l_ambB        = True # include some ambient magnetic field b_z
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
        l_multi_twist = True # multiple flux tubes as described in GFE (2014)
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'paper2c':
        l_ambB        = True # include some ambient magnetic field b_z
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
        l_multi_bps   = True # multiple flux tubes as described in GFE (2014)
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'paper2d':
        l_ambB        = True # include some ambient magnetic field b_z
        l_B0_expz     = True # vertical strength of Bz(r=0) follows exponentials
        l_multi_lanes = True # multiple flux tubes as described in GFE (2014)
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'hmi_model':
        l_B0_quadz    = True # vertical strength of Bz(r=0) follows polynomials 
                           # + coronal exponential    
        l_single      = True  # only one flux tube
        l_hmi         = True  # construct photopheric map of Bz from HMI/SDI
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere
    if model == 'loop_model':
        l_B0_quadz    = True # vertical strength of Bz(r=0) follows polynomials 
                           # + coronal exponential    
        l_single      = True  # only one flux tube
        l_2D_loop     = True # make a 2D loop with sinusoidal Bz(x,0,0)
        l_atmos_val3c_mtw = True # interpolate composite VAL3c+MTW atmosphere

    logical_pars = {
                    'l_mpi':            l_mpi               ,
                    'l_SI':             l_SI                , 
                    'l_CGS':            l_CGS               , 
                    'l_hdonly':         l_hdonly            , 
                    'l_ambB':           l_ambB              , 
                    'l_spruit':         l_spruit            ,
                    'l_const':          l_const             ,         
                    'l_sqrt':           l_sqrt              ,         
                    'l_linear':         l_linear            ,         
                    'l_square':         l_square            ,          
                    'l_B0_expz':        l_B0_expz           ,   
                    'l_B0_quadz':       l_B0_quadz          ,   
                    'l_single':         l_single            ,  
                    'l_hmi':            l_hmi               ,   
                    'l_tube_pair':      l_tube_pair         ,     
                    'l_multi_bps':      l_multi_bps         ,
                    'l_multi_lanes':    l_multi_lanes       ,
                    'l_multi_twist':    l_multi_twist       ,
                    'l_2D_loop':        l_2D_loop           ,
                    'l_mfe':            l_mfe               , 
                    'l_atmos_val3c_mtw':l_atmos_val3c_mtw   ,
                    'suffix':suffix
                    }
    return logical_pars
