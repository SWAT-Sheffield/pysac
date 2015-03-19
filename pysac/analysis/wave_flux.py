# -*- coding: utf-8 -*-
r"""
Routines for the calculation of "Wave Energy Flux" as described in (Bogdan 2003).

This module used the following equations (from Bogdan 2003):

.. math::

    \vec{F}_{wave} \equiv \widetilde{p}_k \vec{v} + \frac{1}{\mu_0} \left(\vec{B}_b \cdot \vec{\widetilde{B}}\right) \vec{v} - \frac{1}{\mu_0}\left(\vec{v} \cdot \vec{\widetilde{B}} \right) \vec{B}_b
    
    F = pk*v + 1/mu[(Bb . B)v] * 1/mu[(v.B) Bb]
"""
import yt
import numpy as np

import pysac.io.yt_fields

__all__ = ['get_wave_flux', 'get_wave_flux_yt']

def get_wave_flux(f, pk):
    """
    Calculate the wave energy flux from a SACData instance and the kinetic pressure
        
    Parameters
    ----------
    f: VACdata instance
        RAW data file
    
    pk: np.ndarray
        Thermal pressure
    
    Returns
    -------
    Fwave: np.ndarray
        The wave energy flux
    """
    Bb = np.array([f.w_dict['bg3'], f.w_dict['bg2'], f.w_dict['bg1']])
    Bp = np.array([f.w_dict['b3'], f.w_dict['b2'], f.w_dict['b1']])
    V = np.array([f.w_sac['v3'], f.w_sac['v2'], f.w_sac['v1']])
    
    #Calculate wave flux
    Fp = 0.25*np.pi * (np.sum(Bb*Bp, axis=0)[None] * V) - (np.sum(V*Bp, axis=0)[None] * Bb)
    Fa = pk[None]*V
    
    Fwave = Fa + Fp
    return Fwave

def get_wave_flux_yt(ds):  # , B_to_SI=1, V_to_SI=1, Pk_to_SI=1):
    """
    Calculate the wave energy flux from a yt dataset.
    
    Parameters
    ----------
    ds: yt dataset
        with derived fields
    
    Returns
    -------
    Fwave: np.ndarray
        The wave energy flux
    """
    cg = ds.index.grids[0]
    
    Bb = np.array([cg['mag_field_x_bg'], cg['mag_field_y_bg'], cg['mag_field_z_bg']])
    Bp = np.array([cg['mag_field_x_pert'], cg['mag_field_y_pert'], cg['mag_field_z_pert']])
    V = np.array([cg['velocity_x'], cg['velocity_y'], cg['velocity_z']])
    
    #Bb *= B_to_SI
    #Bp *= B_to_SI
    Pk = cg['thermal_pressure']# * Pk_to_SI
    
    #Calculate wave flux
    Fp = 0.25*np.pi * (np.sum(Bb*Bp, axis=0)[None] * V) - (np.sum(V*Bp, axis=0)[None] * Bb)
    Fa = Pk[None]*V
    
    Fwave = Fa + yt.YTArray(Fp, 'Pa')  # Cast to Pa
    return Fwave
