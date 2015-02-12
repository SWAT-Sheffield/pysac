# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: stuart
"""
import os
import glob
import numpy as np
from yt.mods import *
import astropy.units as u
import yt
from pysac.mhs_atmosphere.parameters.model_pars import spruit as model_pars
from yt.visualization.api import Streamlines
import pysac.mhs_atmosphere as atm

l_mpi=False
scales, physical_constants = \
    atm.get_parameters()
coords = atm.get_coords(model_pars['Nxyz'], u.Quantity(model_pars['xyz']))
option_pars = atm.set_options(model_pars, l_mpi, l_gdf=True)
if option_pars['l_spruit']:
    option_pars['l_linear'] = True
    if option_pars['l_const']:
        model_pars['chrom_scale'] *= 1.
#        model_pars['xyz'][5] *= 1.
    if option_pars['l_sqrt']:
        model_pars['chrom_scale'] *= 1.
#        model_pars['xyz'][5] *= 2e4
    elif option_pars['l_linear']:
        model_pars['chrom_scale'] *= 1.
# #       model_pars['xyz'][5] *= 10.75
    elif option_pars['l_square']:
        model_pars['chrom_scale'] *= 1.
#        model_pars['xyz'][5] *= 1.
    else:
        model_pars['chrom_scale'] *= 1.
#        model_pars['xyz'][5] *= 1.
    pressure_Z, rho_Z, Rgas_Z = atm.get_spruit_hs(coords['Z'],
                                                  model_pars,
                                                  physical_constants,
                                                  option_pars
                                                  )
datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                                       model_pars['model']+'/')

files = glob.glob(datadir+'/*')
files.sort()

print(files)


ds = yt.load(files[0])
zcent = (ds.parameters['domain_right_edge'][2]+ds.parameters['domain_left_edge'][2])/2.
nx,ny,nz = (ds.domain_dimensions[0],
            ds.domain_dimensions[1],
            ds.domain_dimensions[2])
nxseeds, nzseeds = 8, 1

slc = yt.SlicePlot(ds, fields='mag_field_z_bg', center=[0., 0., zcent],
                   origin="native", normal='x')
slc.set_cmap('mag_field_z_bg', 'YlOrBr')
slc2 = yt.SlicePlot(ds, fields='density_bg', center=[0., 0., zcent],
                   origin="native", normal='x')
slc2.set_cmap('density_bg', 'YlOrBr')
#c = np.array([0.0e6,0.0e6,4.0e6])
#N = 10
#scale = 1.0e6
#pos_dx = np.random.random((N,3))*scale-scale/2.
#pos = c+pos_dx

#pos = np.array([nx/2,ny/2,nz/40])
##pos = zip([np.linspace(nz/40,nz/40,nzseeds)]*nxseeds,
##                  np.linspace(0,nx/10-1,nxseeds))
#streamlines = Streamlines(ds, pos,
#                  'mag_field_x_bg',
#                  'mag_field_y_bg',
#                  'mag_field_z_bg', get_magnitude=True)

figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+
                                                       model_pars['model']+'/')
if not os.path.exists(figsdir):
    os.makedirs(figsdir)
figname = model_pars['model']+'Bz.eps'
slc.save(figsdir+figname)
figname2 = model_pars['model']+'density.eps'
slc2.save(figsdir+figname2)
