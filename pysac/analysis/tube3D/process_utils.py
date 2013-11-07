# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 10:55:54 2012

@author: Stuart Mumford

3D Visualisation and Analysis Utils
"""
import os
import copy
import numpy as np

import tvtk_tube_functions as ttf
from mayavi.tools.sources import vector_field, scalar_field

def read_argv(argv):
    """
    Processes script arguments
    
    arguments are in the form:
        driver post_amp period raw_amp tube_r exp_fac [data_prefix] [h5_dir]
    
    
    Arguments
    ---------
    argv list:
            Script command args
    
    Returns
    -------
    either:
    data_dir, h5_dir, [driver, post_amp, period, raw_amp]
    
    or:
    [driver, post_amp, period, raw_amp]
    
    dependant upon the length of argv    
    """
    
    driver = argv[1]
    post_amp = argv[2]
    period = argv[3]
    raw_amp = argv[4]
    tube_r = argv[5]
    exp_fac = argv[6]
    
    if len(argv) >8:
        data_prefix = argv[7]
        h5_dir = argv[8]
        
        if driver in ['Slog', 'Sarch', 'Suni']:
            data_dir = os.path.join(data_prefix,"%s/%s_%s_%s_%s/"%(driver,period,post_amp,tube_r,exp_fac))
        else:
            data_dir = os.path.join(data_prefix,"%s/%s_%s_%s/"%(driver,period,post_amp,tube_r))
            
        return data_dir, h5_dir, [driver, post_amp, period,
                                  raw_amp, tube_r, exp_fac]
    else:
        return [driver, post_amp, period, raw_amp, tube_r, exp_fac]

def get_hdf5_tvtk(f, cube_slice, flux=True):
    """
    Reads in useful variables from a hdf5 file to tvtk data structures
    
    Arguments
    ---------
    f   hdf5 file handle:
        SAC HDF5 file
    flux    boolean:
        Read variables for fluc calculation?
    cube_slice  np.slice:
        Slice to apply to the arrays
    
    Returns
    -------
    if flux:
        bfield, vfield, density, valf, cs, beta
    else:
        bfield, vfield        
    """
    
    #Do this before convert_B
    if flux:
        va_f = f.get_va()
        cs_f = f.get_cs()
        thermal_p,mag_p = f.get_thermalp(beta=True)
        beta_f = mag_p / thermal_p

    #Convert B to Tesla
    f.convert_B()
    
    # Create TVTK datasets
    bfield = ttf.vector_field(f.w_sac['b3'][cube_slice] * 1e3,
                              f.w_sac['b2'][cube_slice] * 1e3, 
                              f.w_sac['b1'][cube_slice] * 1e3)
    
    vfield = ttf.vector_field(f.w_sac['v3'][cube_slice] / 1e3,
                              f.w_sac['v2'][cube_slice] / 1e3,
                              f.w_sac['v1'][cube_slice] / 1e3)
    
    if flux:
        density = ttf.scalar_field(f.w_sac['rho'][cube_slice])
                                             
        valf = ttf.scalar_field(va_f)
        cs = ttf.scalar_field(cs_f)
        beta = ttf.scalar_field(beta_f)
        
        return bfield, vfield, density, valf, cs, beta
    
    else:
        return bfield, vfield

def get_hdf5_mlab(f, cube_slice, flux=True):
    """
    Reads in useful variables from a hdf5 file to vtk data structures
    
    Arguments
    ---------
    f   hdf5 file handle:
        SAC HDF5 file
    flux    boolean:
        Read variables for flux calculation?
    cube_slice  np.slice:
        Slice to apply to the arrays
    
    Returns
    -------
    if flux:
        bfield, vfield, density, valf, cs, beta
    else:
        bfield, vfield        
    """
    
    #Do this before convert_B
    if flux:
        va_f = f.get_va()
        cs_f = f.get_cs()
        thermal_p,mag_p = f.get_thermalp(beta=True)
        beta_f = mag_p / thermal_p

    #Convert B to Tesla
    f.convert_B()
    
    # Create TVTK datasets
    bfield = vector_field(f.w_sac['b3'][cube_slice] * 1e3,
                                        f.w_sac['b2'][cube_slice] * 1e3, 
                                        f.w_sac['b1'][cube_slice] * 1e3,
                                        name="Magnetic Field",figure=None)
    
    vfield = vector_field(f.w_sac['v3'][cube_slice] / 1e3,
                                        f.w_sac['v2'][cube_slice] / 1e3,
                                        f.w_sac['v1'][cube_slice] / 1e3,
                                        name="Velocity Field",figure=None)
    
    if flux:
        density = scalar_field(f.w_sac['rho'][cube_slice],
                                             name="Density", figure=None)
                                             
        valf = scalar_field(va_f, name="Alven Speed", figure=None)
        cs = scalar_field(cs_f, name="Sound Speed", figure=None)
        beta = scalar_field(beta_f, name="Beta", figure=None)
        
        return bfield, vfield, density, valf, cs, beta
    
    else:
        return bfield, vfield

def get_yt_mlab(ds, cube_slice, flux=True):
    """
    Reads in useful variables from yt to vtk data structures
    
    Arguments
    ---------
    ds: yt dataset
        with derived fields
    flux    boolean:
        Read variables for flux calculation?
    cube_slice  np.slice:
        Slice to apply to the arrays
    
    Returns
    -------
    if flux:
        bfield, vfield, density, valf, cs, beta
    else:
        bfield, vfield        
    """
    cg = ds.h.grids[0]
    
    # Create TVTK datasets
    bfield = vector_field(cg['mag_field_x'][cube_slice] * 1e3,
                          cg['mag_field_y'][cube_slice] * 1e3, 
                          cg['mag_field_z'][cube_slice] * 1e3,
                          name="Magnetic Field",figure=None)
    
    vfield = vector_field(cg['velocity_x'][cube_slice] / 1e3,
                          cg['velocity_y'][cube_slice] / 1e3,
                          cg['velocity_z'][cube_slice] / 1e3,
                          name="Velocity Field",figure=None)
    
    if flux:
        density = scalar_field(cg['density'][cube_slice],
                                             name="Density", figure=None)
                                             
        valf = scalar_field(cg['alfven_speed'], name="Alven Speed", figure=None)
        cs = scalar_field(cg['sound_speed'], name="Sound Speed", figure=None)
        beta = scalar_field(cg['plasma_beta'], name="Beta", figure=None)
        
        return bfield, vfield, density, valf, cs, beta
    
    else:
        return bfield, vfield

#def get_hdf5(f, cube_slice, flux=True, method='mlab'):
#    """
#    Reads in useful variables from a hdf5 file to tvtk data structures
#    
#    Arguments
#    ---------
#    f   hdf5 file handle:
#        SAC HDF5 file
#    flux    boolean:
#        Read variables for fluc calculation?
#    cube_slice  np.slice:
#        Slice to apply to the arrays
#    method: 'mlab', 'tvtk'
#        method to use to read data
#    
#    Returns
#    -------
#    if flux:
#        bfield, vfield, density, valf, cs, beta
#    else:
#        bfield, vfield        
#    """
#    if method == 'mlab':
#        return get_hdf5_mlab(f, cube_slice, flux=flux)
#    
#    if method == 'tvtk':
#        return get_hdf5_tvtk(f, cube_slice, flux=flux)
#    
#    raise ValueError("Invalid Method")
#    
#def process_next_step_tvtk(f, cube_slice, bfield, vfield, density, valf, cs, beta):
#    """ Update all vtk arrays from current file state including flux"""
#    va_f = f.get_va()
#    cs_f = f.get_cs()
#    thermal_p,mag_p = f.get_thermalp(beta=True)
#    beta_f = mag_p / thermal_p
#    density_f = f.w_sac['rho']
#    f.convert_B()
#    
#    # Update Datasets
#    bfield = ttf.vector_field(f.w_sac['b3'][cube_slice] * 1e3,
#                              f.w_sac['b2'][cube_slice] * 1e3,
#                              f.w_sac['b1'][cube_slice] * 1e3)
#                              
#    bfield = ttf.vector_field(f.w_sac['v3'][cube_slice] * 1e3,
#                              f.w_sac['v2'][cube_slice] * 1e3,
#                              f.w_sac['v1'][cube_slice] * 1e3)
#      
#                        
#    density = ttf.scalar_field(density_f)                                   
#    valf = ttf.scalar_field(va_f)
#    cs = ttf.scalar_field(cs_f)
#    beta = ttf.scalar_field(beta_f)
#    
#    return bfield, vfield, density, valf, cs, beta
def process_next_step_yt(ds, cube_slice, bfield, vfield, density, valf, cs, beta):
    """ Update all vtk arrays from current file state including flux"""
    cg = ds.h.grids[0]
    
    # Update Datasets
    bfield.set(vector_data = np.rollaxis(np.array([cg['mag_field_x'][cube_slice] * 1e3,
                                                   cg['mag_field_y'][cube_slice] * 1e3,
                                                   cg['mag_field_z'][cube_slice] * 1e3]),
                                                   0, 4))
    vfield.set(vector_data = np.rollaxis(np.array([cg['velocity_x'][cube_slice] / 1e3,
                                                   cg['velocity_y'][cube_slice] / 1e3,
                                                   cg['velocity_z'][cube_slice] / 1e3]),
                                                   0, 4))
    valf.set(scalar_data = cg['alfven_speed'])
    cs.set(scalar_data = cg['sound_speed'])
    beta.set(scalar_data = cg['plasma_beta'])
    density.set(scalar_data = cg['density'])
    
    return bfield, vfield, density, valf, cs, beta

def process_next_step(f, cube_slice, bfield, vfield, density, valf, cs, beta):
    """ Update all vtk arrays from current file state including flux"""
    va_f = f.get_va()
    cs_f = f.get_cs()
    thermal_p,mag_p = f.get_thermalp(beta=True)
    beta_f = mag_p / thermal_p
    density_f = f.w_sac['rho']
    f.convert_B()
    
    # Update Datasets
    bfield.set(vector_data = np.rollaxis(np.array([f.w_sac['b3'][cube_slice] * 1e3,
                                                   f.w_sac['b2'][cube_slice] * 1e3,
                                                   f.w_sac['b1'][cube_slice] * 1e3]),
                                                   0, 4))
    vfield.set(vector_data = np.rollaxis(np.array([f.w_sac['v3'][cube_slice] / 1e3,
                                                   f.w_sac['v2'][cube_slice] / 1e3,
                                                   f.w_sac['v1'][cube_slice] / 1e3]),
                                                   0, 4))
    valf.set(scalar_data = va_f)
    cs.set(scalar_data = cs_f)
    beta.set(scalar_data = beta_f)
    density.set(scalar_data = density_f)
    
    return bfield, vfield, density, valf, cs, beta
    
def process_next_step_mlab(f, cube_slice):
    """ Update all vtk arrays from current file state including flux"""
    #Do this before convert_B
    va_f = f.get_va()
    cs_f = f.get_cs()
    thermal_p,mag_p = f.get_thermalp(beta=True)
    beta_f = mag_p / thermal_p

    #Convert B to Tesla
    f.convert_B()
    
    # Create TVTK datasets
    bfield = vector_field(f.w_sac['b3'][cube_slice] * 1e3,
                          name="Magnetic Field",figure=None)
    
    vfield = vector_field(f.w_sac['v3'][cube_slice] / 1e3,
                          f.w_sac['v2'][cube_slice] / 1e3,
                          f.w_sac['v1'][cube_slice] / 1e3,
                          name="Velocity Field",figure=None)
    
    density = scalar_field(f.w_sac['rho'][cube_slice],
                           name="Density", figure=None)
                                         
    valf = scalar_field(va_f, name="Alven Speed", figure=None)
    cs = scalar_field(cs_f, name="Sound Speed", figure=None)
    beta = scalar_field(beta_f, name="Beta", figure=None)
    
    return bfield, vfield, density, valf, cs, beta
#
#def process_next_step(f, cube_slice, bfield, vfield, density, valf, cs, beta, method='mlab'):
#    if method == 'mlab':
#        return process_next_step_mlab(f, cube_slice)
#    
#    if method == 'tvtk':
#        return process_next_step_tvtk(f, cube_slice, bfield, vfield, density, valf, cs, beta)
#    
#    raise ValueError("Invalid Method")
