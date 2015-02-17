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
import matplotlib.colorbar as cb
import matplotlib.pyplot as plt
from matplotlib import cm
import yt
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
# loop over all four models
for spruit in spruits:
    datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                 spruit+'/')
    # open all gdf files in the model directory
    files = glob.glob(datadir+'/*')
    files = glob.glob(datadir+'/'+spruits[0]+'_3Daux.gdf')
    files.sort()

    print(files)

    for file_ in files:
        ds = yt.load(file_)
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
#        pos[:,0] = xcent
        #pos = zip([np.linspace(nz/40,nz/40,nzseeds)]*nxseeds,
        #                  np.linspace(0,nx/10-1,nxseeds))
#        streamlines = Streamlines(ds, pos,
#                                  'mag_field_x_bg',
#                                  'mag_field_y_bg',
#                                  'mag_field_z_bg')

        vars_ = ds.index.field_list
        for var_ in vars_:
            var_field = var_[1]
            max_var = np.max(np.abs(ds.index.grids[0][var_field]))
            var = ds.index.grids[0][var_field]
            if max_var > 0.:
                oneD_arrays = atm.make_1d_slices(ds, var_field, oneD_arrays)
#                plt.figure(figsize=[3.5,5.60])
#                fig, axes, cbars = get_multi_plot(1,1,colorbar='vertical',bw=1)
                plots = []
                slc = yt.SlicePlot(ds, fields=var_field,
                                   center=[0., 0., zcent],
                                   normal='x', origin='native'
                                   )#, width = (zplot, 'm'))
                slc.figure_size *= 1.
                slc.set_axes_unit('Mm')
                slc.aspect = 1.
#                slc.zoom = 2.
                slc.set_font_size(12)
                slc.set_minorticks('all', 'off')
                colour = cm.YlOrBr
                if 'mag' in var_field:
#                    colour = 'RdYlGn_r'
                    colour = 'Blue-Red'
                if 'temperature' in var_field:
                    colour = cm.RdBu_r
                if 'beta' in var_field:
                    colour = cm.PiYG
                slc.set_cmap(var_field, colour)
 #               slc.plots[var_field].hide_colorbar()
#                slc.set_cbar_minorticks('all', 'off')
#                plt.figure(figsize=figxy)
#                plt.imshow(var.T, colour,
#                extent=[ds.parameters['domain_left_edge[0],
#                       ds.parameters['domain_right_edge[0],
#                       ds.parameters['domain_left_edge[2],
#                       ds.parameters['domain_right_edge[2]
#                       ],
#                norm=colors.LogNorm(),
#                origin='lower',
#                vmin=colarr.min(),
#                vmax=colarr.max())
                figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+
                                                                       spruit+'/')
                if not os.path.exists(figsdir):
                    os.makedirs(figsdir)
                figname  = spruit+'_'+var_field+'.eps'
                slc.save(figsdir+figname)

            #c = np.array([0.0e6,0.0e6,4.0e6])
            #N = 10
            #scale = 1.0e6
            #pos_dx = np.random.random((N,3))*scale-scale/2.
            #pos = c+pos_dx


