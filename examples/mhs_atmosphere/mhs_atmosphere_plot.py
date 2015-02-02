# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: stuart
"""
import os
import glob
import numpy as np
from yt.mods import *
import yt
from pysac.mhs_atmosphere.parameters.model_pars import paper1 as model_pars
from yt.visualization.api import Streamlines

datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                                       model_pars['model']+'/')

files = glob.glob(datadir+'/*')
files.sort()

print(files)


ds = yt.load(files[0])
nx,ny,nz = (ds.domain_dimensions[0],
            ds.domain_dimensions[1],
            ds.domain_dimensions[2])
nxseeds, nzseeds = 8, 1

slc = yt.SlicePlot(ds, fields='density_bg', center=[0., 0., 4.31e6],
                   origin="native", normal='x')
slc.set_cmap('density_bg', 'YlOrBr')
c = np.array([0.0e6,0.0e6,4.0e6])
N = 10
scale = 1.0e6
pos_dx = np.random.random((N,3))*scale-scale/2.
pos = c+pos_dx

#pos = np.array([nx/2,ny/2,nz/40])
##pos = zip([np.linspace(nz/40,nz/40,nzseeds)]*nxseeds,
##                  np.linspace(0,nx/10-1,nxseeds))
streamlines = Streamlines(ds, pos,
                  'mag_field_x_bg',
                  'mag_field_y_bg',
                  'mag_field_z_bg', get_magnitude=True)

figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+
                                                       model_pars['model']+'/')
if not os.path.exists(figsdir):
    os.makedirs(figsdir)
figname = 'density.eps'
slc.save(figsdir+figname)
