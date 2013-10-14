import numpy as np

def get_wave_flux(f, pk):
    """
    This routine calculates the wave energy flux as defined in (Bogdan 2003)
    
    $\vec{F}_{wave} \equiv \widetilde{p}_k \vec{v} + \frac{1}{\mu_0} \left(\vec{B}_b \cdot \vec{\widetilde{B}}\right) \vec{v} - \frac{1}{\mu_0}\left(\vec{v} \cdot \vec{\widetilde{B}} \right) \vec{B}_b$
    
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