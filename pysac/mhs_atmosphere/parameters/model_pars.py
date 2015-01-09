# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 17:45:48 2014

@author: sm1fg
"""

import numpy as np
import astropy.units as u

hmi_model = {'photo_scale': 0.6*u.Mm,
             'chrom_scale': 0.31*u.Mm,
             'corona_scale': 100*u.Mm,                    #scale height for the corona
             'coratio': 0.06*u.one,
             'phratio': 0.15*u.one,
             'pixel': 0.36562475*u.Mm,              #(HMI pixel)
             'radial_scale': 0.044*u.Mm,
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 4.250e-4*u.T}
hmi_model['chratio'] = 1*u.one - hmi_model['coratio'] - hmi_model['phratio']

mfe_setup = {'photo_scale': 0.6*u.Mm,
             'chrom_scale': 0.31*u.Mm,
             'corona_scale': 100*u.Mm,                    #scale height for the corona
             'coratio': 0.075*u.one,
             'phratio': 0.0*u.one,
             'pixel': 0.36562475*u.Mm,              #(HMI pixel)
             'radial_scale': 0.044*u.Mm,
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 4.250e-4*u.T}
mfe_setup['chratio'] = 1*u.one - mfe_setup['coratio'] - mfe_setup['phratio']

spruit = {'photo_scale': 1.5*u.Mm,
          'chrom_scale': 1.5*u.Mm,
          'corona_scale': 100*u.Mm,                    #scale height for the corona
          'coratio': 0.0*u.one,
          'phratio': 0.0*u.one,
          'pixel': 0.1*u.Mm,              #(HMI pixel)
          'radial_scale': 0.025*u.Mm,
          'nftubes': 1,
          'B_corona': 0.*u.T,
          'pBplus': 4.250e-4*u.T}
spruit['chratio'] = 1*u.one - spruit['coratio'] - spruit['phratio']

paper1 = {'photo_scale': 0.6*u.Mm,
          'chrom_scale': 0.42*u.Mm,
          'corona_scale': 175*u.Mm,         #scale height for the corona
          'coratio': 0.0225*u.one,
          'phratio': 0.0*u.one,
          'pixel': 0.36562475*u.Mm,              #(HMI pixel)
          'radial_scale': 0.044*u.Mm,
          'nftubes': 1,
          'B_corona': 2.00875e-4*u.T,
          'pBplus': 4.250e-4*u.T}
paper1['chratio'] = 1*u.one - paper1['coratio'] - paper1['phratio']


def get_coords(Nxyz, xyz):
    """
    get_coords returns a non-dimensional dictionary describing the domain
    coordinates.
    """

    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = np.linspace(xyz[4],
                        xyz[5],Nxyz[2])
    Zext = np.linspace(Z.min()-4.*dz, Z.max()+4.*dz, Nxyz[2]+8)
    coords = {
              'dx':(xyz[1]-xyz[0])/(Nxyz[0]-1),
              'dy':(xyz[3]-xyz[2])/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xyz[0],
              'xmax':xyz[1],
              'ymin':xyz[2],
              'ymax':xyz[3],
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords
