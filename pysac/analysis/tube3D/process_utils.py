# -*- coding: utf-8 -*-
"""
This module contains routines to convert from SACData and yt DataSets into
tvtk fields via mayavi.
"""
import numpy as np

from mayavi.tools.sources import vector_field, scalar_field

__all__ = ['get_sacdata_mlab', 'get_yt_mlab', 'process_next_step_yt',
           'process_next_step_sacdata']

def get_sacdata_mlab(f, cube_slice, flux=True):
    """
    Reads in useful variables from a hdf5 file to vtk data structures
    
    Parameters
    ----------
    f : hdf5 file handle
        SAC HDF5 file

    flux : boolean
        Read variables for flux calculation?
    
    cube_slice : np.slice
        Slice to apply to the arrays
    
    Returns
    -------
    bfield, vfield[, density, valf, cs, beta]
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

def process_next_step_yt(ds, cube_slice, bfield, vfield, density, valf, cs, beta):
    """
    Update all vtk arrays from current file state including flux
    
    """
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

def process_next_step_sacdata(f, cube_slice, bfield, vfield, density, valf, cs, beta):
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