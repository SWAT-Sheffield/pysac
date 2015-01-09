# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 17:45:48 2014

@author: sm1fg
"""

import numpy as np
import sys

def get_model(model, 
              Nxyz,
              xyz_SI,
              scales, 
              logical_pars
             ):

#==============================================================================
    dz=(xyz_SI[5]-xyz_SI[4])/(Nxyz[2]-1)/scales['length']
    Z    = np.linspace(xyz_SI[4]/scales['length'],
                        xyz_SI[5]/scales['length'],Nxyz[2])
    Zext = np.linspace(Z.min()-4.*dz, Z.max()+4.*dz, Nxyz[2]+8)
    coords = {
              'dx':(xyz_SI[1]-xyz_SI[0])/(Nxyz[0]-1)/scales['length'],
              'dy':(xyz_SI[3]-xyz_SI[2])/(Nxyz[1]-1)/scales['length'],
              'dz':(xyz_SI[5]-xyz_SI[4])/(Nxyz[2]-1)/scales['length'],
              'xmin':xyz_SI[0]/scales['length'],
              'xmax':xyz_SI[1]/scales['length'],
              'ymin':xyz_SI[2]/scales['length'], 
              'ymax':xyz_SI[3]/scales['length'],
              'zmin':xyz_SI[4]/scales['length'],
              'zmax':xyz_SI[5]/scales['length'],
              'Z':Z,
              'Zext':Zext
             }
    if model == 'hmi_model':
        save_gdf   = model+logical_pars['suffix']
        z_phot_SI  = 6e5
        z_chrom_SI = 3.1e5
        z_cor_SI   = 1e8                    #scale height for the corona
        coratio    = 0.06 
        phratio    = 0.15
        chratio    = 1 - coratio - phratio
        pixel_SI   = 365624.75              #(HMI pixel)
        f0_SI      = 4.4e4
        nftubes    = 1
        B_corona_SI= 0.
        pBplus_SI  = 4.250e-4
    elif model == 'mfe_setup':
        save_gdf   = model+logical_pars['suffix']
        z_phot_SI  = 6e5
        z_chrom_SI = 3.1e5
        z_cor_SI   = 1e8                    #scale height for the corona
        coratio    = 0.075 
        phratio    = 0.0
        chratio    = 1 - coratio - phratio
        pixel_SI   = 365624.75              #(HMI pixel)
        f0_SI      = 4.4e4
        nftubes    = 1
        B_corona_SI= 0.
        pBplus_SI  = 4.250e-4
    elif model == 'spruit':
        save_gdf   = model+logical_pars['suffix']
        z_phot_SI  = 1.5e6
        z_chrom_SI = 1.5e6
        z_cor_SI   = 1e8                    #scale height for the corona
        coratio    = 0.0 
        phratio    = 0.0
        chratio    = 1 - coratio - phratio
        pixel_SI   = 1e5              #(HMI pixel)
        f0_SI      = 2.5e4
        nftubes    = 1
        B_corona_SI= 0.
        pBplus_SI  = 4.250e-4
    elif model == 'paper1':
        save_gdf   = model+logical_pars['suffix']
        z_phot_SI  = 6e5
        z_chrom_SI = 4.2e5
        z_cor_SI   = 1.75e8                    #scale height for the corona
        coratio    = 0.0225 
        phratio    = 0.0
        chratio    = 1 - coratio - phratio
        pixel_SI   = 365624.75              #(HMI pixel)
        f0_SI      = 4.4e4
        nftubes    = 1
        B_corona_SI= 2.00875e-4
        pBplus_SI  = 4.250e-4
    else:
        sys.exit('in parameters.model_pars.get_model, \
                  set model parameters for model '+model)
    if logical_pars['l_CGS']:
        B_corona_SI, pBplus_SI = B_corona_SI*1e4, pBplus_SI*1e4
        z_phot_SI, z_chrom_SI = z_phot_SI*1e2, z_chrom_SI*1e2
        z_cor_SI, f0_SI, pixel_SI = z_cor_SI*1e2, f0_SI*1e2, pixel_SI*1e2
    model_pars = {
                  'save_gdf':save_gdf,
                  'pixel':pixel_SI/scales['length'],
                  'photo_scale':z_phot_SI/scales['length'], 
                  'chrom_scale':z_chrom_SI/scales['length'], 
                  'corona_scale':z_cor_SI/scales['length'], 
                  'radial_scale':f0_SI/scales['length'], 
                  'phratio':phratio,
                  'coratio':coratio,
                  'chratio':chratio,
                  'B_corona':B_corona_SI/scales['magnetic'],
                  'pBplus':pBplus_SI/scales['magnetic'],
                  'nftubes':nftubes
                 }

#==============================================================================             
    return coords, model_pars             
