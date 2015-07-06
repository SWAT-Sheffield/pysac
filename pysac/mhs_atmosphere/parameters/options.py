# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 11:37:39 2014

@author: sm1fg

"""
##============================================================================
##option parameters
##============================================================================
def set_options(model, l_mpi, l_gdf=True):

    """This module assigns the logical options for the model. If adding 
    new models with additional logical arguments add it to the default 
    list as false, include an if statement for True update the
    dictionary option_pars 
    
    """    
    #default arguments
    option_pars = {
        'l_hdonly':      False,# set mag field zero to check background
        'l_ambB':        False,# include some ambient magnetic field b_z
        'l_spruit':      False,# thin flux tube model to check Spruit 
        'l_const':       False,# axial Alfven speed const  Z-depend (Spruit)
        'l_sqrt':        False,# axial Alfven speed sqrt   Z-depend (Spruit)
        'l_linear':      False,# axial Alfven speed linear Z-depend (Spruit)
        'l_square':      False,# axial Alfven speed square Z-depend (Spruit)
        'l_B0_expz':     False,# Z-depend of Bz(r=0) exponentials
        'l_B0_quadz':    False,# Z-depend of Bz(r=0) polynomials + exponential 
        'l_B0_rootz':    False,# Z-depend of Bz(r=0) sqrt polynomials 
        'l_single':      False,# only one flux tube
        'l_hmi':         False,# construct photopheric map of Bz from HMI/SDI
        'l_tube_pair':   False,# pair of flux tubes
        'l_multi_netwk': False,# multiple flux tubes as described in GFE (2014)
        'l_multi_lanes': False,# multiple flux tubes as described in GFE (2014)
        'l_multi_twist': False,# multiple flux tubes as described in GFE (2014)
        'l_2D_loop':     False,# make a 2D loop with sinusoidal Bz(x,0,0)
        'l_mfe':         False,# model Viktor's model from MFE (2014)      
        'l_atmos_val3c_mtw':False,# interpolate composite VAL3c+MTW atmosphere
        'suffix':        '.gdf'
    }
    #revise optional parameters depending on configuration required
    if model['model'] == 'hmi_model':
        option_pars['l_hmi']             = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'mfe_setup':
        option_pars['l_single']          = True 
        option_pars['l_mfe']             = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'spruit':    
        option_pars['l_single']          = True 
        option_pars['l_spruit']          = True 
    if model['model'] == 'paper1':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True
        option_pars['l_single']          = True
        option_pars['l_atmos_val3c_mtw'] = True
    if model['model'] == 'paper2a':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_tube_pair']       = True 
        option_pars['l_atmos_val3c_mtw'] = True
    if model['model'] == 'paper2b':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_multi_twist'    ] = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'paper2c':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_multi_netwk']     = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'paper2d':
        option_pars['l_ambB']            = True 
        option_pars['l_B0_expz']         = True 
        option_pars['l_multi_lanes'    ] = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'hmi_model':
        option_pars['l_B0_quadz']        = True 
        option_pars['l_single']          = True 
        option_pars['l_hmi']             = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if model['model'] == 'loop_model':
        option_pars['l_B0_quadz']        = True 
        option_pars['l_single']          = True 
        option_pars['l_2D_loop']         = True 
        option_pars['l_atmos_val3c_mtw'] = True 
    if l_mpi:
        option_pars['l_mpi'] = True
    else:
        option_pars['l_mpi'] = False
    return option_pars
