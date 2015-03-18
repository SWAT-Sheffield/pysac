# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: stuart
"""
import os
import glob
import numpy as np
import astropy.units as u
import matplotlib.colorbar as cb
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, FormatStrFormatter
import matplotlib.colors as colors
import pysac.yt as sacyt
#from pysac.mhs_atmosphere.parameters.model_pars import spruit as model_pars
from yt.visualization.streamlines import Streamlines
from yt.visualization.api import get_multi_plot
import pysac.mhs_atmosphere as atm

l_mpi=False
scales, physical_constants = \
    atm.get_parameters()
#define the models required
#spruits = ['spruit_const','spruit_sqrt','spruit_linear','spruit_square']
spruits = ['spruit_const']
oneD_arrays = {}
oned_dataset = []
figxz = [5.5,5.60]
# loop over all four models
for spruit in spruits:
    datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                 spruit+'/')
    figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'
                                              +spruit+'/')
    if not os.path.exists(figsdir):
        os.makedirs(figsdir)
    # open all gdf files in the model directory
    files = glob.glob(datadir+'/*')
#    files = glob.glob(datadir+'/'+spruits[0]+'_3Daux.gdf')
    files.sort()

    print(files)

    for file_ in files:
#        ds = yt.load(file_)
        ds = sacyt.SACGDFDataset(file_)
        xcent = (ds.parameters['domain_right_edge'][0] +
                 ds.parameters['domain_left_edge'][0])/2.
        ycent = (ds.parameters['domain_right_edge'][1] +
                 ds.parameters['domain_left_edge'][1])/2.
        zcent = (ds.parameters['domain_right_edge'][2] +
                 ds.parameters['domain_left_edge'][2])/2.
        zplot = (ds.parameters['domain_right_edge'][2] +
                 ds.parameters['domain_left_edge'][2])/8.
        nx,ny,nz = (ds.domain_dimensions[0],
                    ds.domain_dimensions[1],
                    ds.domain_dimensions[2])
        nxseeds, nzseeds = 8, 1
        c = np.array([xcent,ycent,zcent])
#        c = np.array([0.5]*3)
        N = 10
        scale = np.max(np.abs(ds.parameters['domain_left_edge']))
        pos_dx = np.random.random((N,3))*scale-scale/2.
        pos = c+pos_dx
        vars_ = ds.index.field_list
        lines, contours = False, False
        for var_ in vars_:
            var_field = var_[1]
            max_var = np.max(np.abs(ds.index.grids[0][var_field]))/\
                             ds.index.grids[0][var_field].unit_quantity
            var = ds.index.grids[0][var_field]
            if max_var > 0.:
                oneD_arrays = atm.make_1d_slices(ds, var_field, oneD_arrays)
                figname  = figsdir+spruit+'_'+var_field+'.eps'
                nx_2 = ds.domain_dimensions[1]/2
                if 'mag' in var_field:
                    lines = True
                elif 'density' or 'pressure' in var_field:
                    lines = True
                    if 'D' in file_:
                        lines = False
                else:
                    lines = False
                if '_HS' in var_field:
                    contours = False
                else:
                    contours = True
                atm.make_2d_plot(ds, var_field, figname, 
                                                    normal=['y',nx_2],
                                                    aspect=0.2, lines=lines,
                                                    contours=contours,
                                                    model=spruit)
    plot_label = figsdir+spruit+'_axis.eps'
    keys = ['alfven_speed','sound_speed','mag_field_z_bg']
    subkeys = ['axis']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = True
                                           )
            #c = np.array([0.0e6,0.0e6,4.0e6])
            #N = 10
            #scale = 1.0e6
            #pos_dx = np.random.random((N,3))*scale-scale/2.
            #pos = c+pos_dx


