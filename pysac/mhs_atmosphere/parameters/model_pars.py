# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 17:45:48 2014

@author: sm1fg
"""

import numpy as np
import astropy.units as u
hmi_model = {'photo_scale': 0.6*u.Mm,
             'chrom_scale': 0.31*u.Mm,
             'corona_scale': 100*u.Mm,      #scale height for the corona
             'coratio': 0.06*u.one,
             'model': 'hmi_model',
             'phratio': 0.15*u.one,
             'pixel': 0.36562475*u.Mm,      #(HMI pixel)
             'radial_scale': 0.044*u.Mm,
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 4.250e-4*u.T}
hmi_model['chratio'] = 1*u.one - hmi_model['coratio'] - hmi_model['phratio']

mfe_setup = {'photo_scale': 0.6*u.Mm,
             'chrom_scale': 0.31*u.Mm,
             'corona_scale': 100*u.Mm,  #scale height for the corona
             'coratio': 0.075*u.one,
             'model': 'mfe_setup',
             'phratio': 0.0*u.one,
             'pixel': 0.36562475*u.Mm,  #(HMI pixel)
             'radial_scale': 0.044*u.Mm,
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 4.250e-4*u.T}
mfe_setup['chratio'] = 1*u.one - mfe_setup['coratio'] - mfe_setup['phratio']
#if 1D or 2D set unused dimensions to 0, and unrequired xyz limits to 1.
mfe_setup['Nxyz'] = [128,128,128] # 3D grid
mfe_setup['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.5e-2*u.Mm,1.6*u.Mm] #grid size

spruit = {'photo_scale': 1.5*u.Mm,
          'chrom_scale': 0.5*u.Mm,
          'corona_scale': 100*u.Mm,      #scale height for the corona
          'coratio': 0.0*u.one,
          'model': 'spruit',
          'phratio': 0.0*u.one,
          'pixel': 0.1*u.Mm,              #(HMI pixel)
          'radial_scale': 0.075*u.Mm,
          'nftubes': 1,
          'p0': 117200.0 * u.dyne/u.cm**2,
          'B_corona': 0.*u.T,
          'pBplus': 4.250e-4*u.T}
spruit['chratio'] = 1*u.one - spruit['coratio'] - spruit['phratio']
spruit['Nxyz'] = [128,128,256] # 3D grid
spruit['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.0*u.Mm,25.5*u.Mm] #grid size


paper1 = {'photo_scale': 0.6*u.Mm,
          'chrom_scale': 0.42*u.Mm,
          'corona_scale': 175*u.Mm,         #scale height for the corona
          'coratio': 0.0225*u.one,
          'model': 'paper1',
          'phratio': 0.0*u.one,
          'pixel': 0.36562475*u.Mm,              #(HMI pixel)
          'radial_scale': 0.044*u.Mm,
          'nftubes': 1,
          'B_corona': 2.00875e-4*u.T,
          'pBplus': 4.250e-4*u.T}
paper1['chratio'] = 1*u.one - paper1['coratio'] - paper1['phratio']
paper1['Nxyz'] = [128,128,432] # 3D grid
paper1['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size

paper2a = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.42*u.Mm,
           'corona_scale': 175*u.Mm,         #scale height for the corona
           'coratio': 0.0225*u.one,
           'model': 'paper2a',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.044*u.Mm,
           'nftubes': 4,
           'B_corona': 2.00875e-4*u.T,
           'pBplus': 4.250e-4*u.T}
paper2a['chratio'] = 1*u.one - paper2a['coratio'] - paper2a['phratio']
paper2a['Nxyz'] = [160,80,432] # 3D grid
paper2a['xyz']  = [-1.59*u.Mm,1.59*u.Mm,-0.79*u.Mm,0.79*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size

paper2b = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.42*u.Mm,
           'corona_scale': 175*u.Mm,         #scale height for the corona
           'coratio': 0.0225*u.one,
           'model': 'paper2b',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.044*u.Mm,
           'nftubes': 4,
           'B_corona': 2.00875e-4*u.T,
           'pBplus': 4.250e-4*u.T}
paper2b['chratio'] = 1*u.one - paper2b['coratio'] - paper2b['phratio']
paper2b['Nxyz'] = [48,48,140] # 3D grid
paper2b['xyz']  = [-0.47*u.Mm,0.47*u.Mm,-0.47*u.Mm,0.47*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2c = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.42*u.Mm,
           'corona_scale': 175*u.Mm,         #scale height for the corona
           'coratio': 0.0225*u.one,
           'model': 'paper2c',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.044*u.Mm,
           'nftubes': 15,
           'B_corona': 2.00875e-4*u.T,
           'pBplus': 4.250e-4*u.T}
paper2c['chratio'] = 1*u.one - paper2c['coratio'] - paper2c['phratio']
paper2c['Nxyz'] = [224,224,140] # 3D grid
paper2c['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-2.23*u.Mm,2.23*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2d = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.42*u.Mm,
           'corona_scale': 175*u.Mm,         #scale height for the corona
           'coratio': 0.0225*u.one,
           'model': 'paper2d',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.044*u.Mm,
           'nftubes': 18,
           'B_corona': 2.00875e-4*u.T,
           'pBplus': 4.250e-4*u.T}
paper2d['chratio'] = 1*u.one - paper2d['coratio'] - paper2d['phratio']
paper2d['Nxyz'] = [224,224,140] # 3D grid
paper2d['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-0.79*u.Mm,0.79*u.Mm,0*u.Mm,2.78*u.Mm] #grid size


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
