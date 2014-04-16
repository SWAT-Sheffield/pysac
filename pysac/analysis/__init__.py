"""
General routines for the analysis of HD and MHD simulations run with SAC
"""

import numpy as np

from . import tube3D

__all__ = ['get_wave_flux', 'get_wave_flux_yt', 'tube3D']

def get_wave_flux(f, pk):
    r"""
    This routine calculates the wave energy flux as defined in (Bogdan 2003).
    
    .. math::
    
        \vec{F}_{wave} \equiv \widetilde{p}_k \vec{v} + \frac{1}{\mu_0} \left(\vec{B}_b \cdot \vec{\widetilde{B}}\right) \vec{v} - \frac{1}{\mu_0}\left(\vec{v} \cdot \vec{\widetilde{B}} \right) \vec{B}_b
        
        F = pk*v + 1/mu[(Bb . B)v] * 1/mu[(v.B) Bb]
        
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

def get_wave_flux_yt(ds):
    r"""
    This routine calculates the wave energy flux as defined in (Bogdan 2003).
    
    .. math::
    
        \vec{F}_{wave} \equiv \widetilde{p}_k \vec{v} + \frac{1}{\mu_0} \left(\vec{B}_b \cdot \vec{\widetilde{B}}\right) \vec{v} - \frac{1}{\mu_0}\left(\vec{v} \cdot \vec{\widetilde{B}} \right) \vec{B}_b
        
        F = pk*v + 1/mu[(Bb . B)v] * 1/mu[(v.B) Bb]
    
    Parameters
    ----------
    ds: yt dataset
        with derived fields
    
    Returns
    -------
    Fwave: np.ndarray
        The wave energy flux
    """
    cg = ds.h.grids[0]
    
    Bb = np.array([cg['mag_field_x_bg'], cg['mag_field_y_bg'], cg['mag_field_z_bg']])
    Bp = np.array([cg['mag_field_x_pert'], cg['mag_field_y_pert'], cg['mag_field_z_pert']])
    V = np.array([cg['velocity_x'], cg['velocity_y'], cg['velocity_z']])
    
    #Calculate wave flux
    Fp = 0.25*np.pi * (np.sum(Bb*Bp, axis=0)[None] * V) - (np.sum(V*Bp, axis=0)[None] * Bb)
    Fa = cg['thermal_pressure'][None]*V
    
    Fwave = Fa + Fp
    return Fwave