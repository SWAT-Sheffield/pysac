# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 11:33:12 2015

@author: sm1fg
"""

import os
import yt
#import glob
#import h5py
from pysac.mhs_atmosphere.parameters.model_pars import paper1 as model_pars
import pysac.mhs_atmosphere as atm

datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                       model_pars['model']+'/')
filename = datadir + model_pars['model'] + '.gdf'
if not os.path.exists(datadir):
    os.makedirs(datadir)
sourcefile = datadir + model_pars['model'] + '_sources' + '.gdf'
aux3D = datadir + model_pars['model'] + '_3Daux' + '.gdf'
aux1D = datadir + model_pars['model'] + '_1Daux' + '.gdf'
#files = glob.glob(datadir+'/*')
#files.sort()

#==============================================================================
#1D slices
#==============================================================================
ds1D = yt.load(aux1D)
#ds1D = yt.load(aux1D)
slc1D = ds1D.index.grids[0]
ds3D = yt.load(aux3D)
#ds1D = yt.load(aux1D)
slc3D = ds3D.index.grids[0]
dsSAC = yt.load(filename)
#ds1D = yt.load(aux1D)
slcSAC = dsSAC.index.grids[0]
empirical_data = atm.read_VAL3c_MTW()




#tmps = yt.load(filename)
#sac = tmps.covering_grid(level=0, left_edge=[0,0.0,0.0],
#                                      dims=tmps.domain_dimensions)


figsdir = os.path.expanduser('~/mhs_atmosphere/figs/'+model_pars['model']+'/')
if not os.path.exists(figsdir):
    os.makedirs(figsdir)