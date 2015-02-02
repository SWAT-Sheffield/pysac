# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 11:33:12 2015

@author: sm1fg
"""

import os
import yt
import glob


from pysac.mhs_atmosphere.parameters.model_pars import paper1 as model_pars

datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+model_pars['model']+'/')
filename = datadir + model_pars['model'] + '.gdf'
if not os.path.exists(datadir):
    os.makedirs(datadir)
sourcefile = datadir + model_pars['model'] + '_sources' + '.gdf'
auxfile = datadir + model_pars['model'] + '_aux' + '.gdf'
files = glob.glob(datadir+'/*')
files.sort()

tmpa = yt.load(auxfile)
aux = tmpa.covering_grid(level=0, left_edge=[0,0.0,0.0],
                                      dims=tmpa.domain_dimensions)
tmps = yt.load(filename)
sac = tmps.covering_grid(level=0, left_edge=[0,0.0,0.0],
                                      dims=tmps.domain_dimensions)


figsdir = os.path.expanduser('~/mhs_atmosphere/figs/'+model_pars['model']+'/')
if not os.path.exists(figsdir):
    os.makedirs(figsdir)