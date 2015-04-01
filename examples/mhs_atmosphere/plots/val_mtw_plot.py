# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: fred
"""
import os
import glob
import numpy as np
import pysac.yt as sacyt
import pysac.mhs_atmosphere as atm

l_mpi=False
scales, physical_constants = \
    atm.get_parameters()
#define the models required
papers = ['paper1','paper2a','paper2b','paper2c','paper2d','mfe_setup']
#papers = ['paper1']
oneD_arrays = {}
oned_dataset = []
#loop over all four models
for paper in papers:
    datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                 paper+'/')
    figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+paper+'/')
    if not os.path.exists(figsdir):
        os.makedirs(figsdir)
    #open all gdf files in the model directory
    files = glob.glob(datadir+'/*')
    #files = glob.glob(datadir+'/'+papers[0]+'_3Daux.gdf')
    files.sort()

    print(files)

    for file_ in files:
        ds = sacyt.SACGDFDataset(file_)
#        fkeys = ['thermal_pressure','plasma_beta',
        vars_ = ds.index.field_list
        for var_ in vars_:
            var_field = var_[1]
            max_var = np.max(np.abs(ds.index.grids[0][var_field]))/\
                             ds.index.grids[0][var_field].unit_quantity
            var = ds.index.grids[0][var_field]
            if max_var > 0.:
                # save 1D slices from each variable for plotting
                oneD_arrays = atm.make_1d_slices(ds, var_field, oneD_arrays)
                # select the central slice to plot normal to the y-plane
                plane, N_2 = 'y', ds.domain_dimensions[1]/2
                lines, contours = True, True
                if '_HS' in var_field:
                    lines, contours = False, False
                elif '1D' in file_:
                    lines, contours = False, False
                else:
                    lines, contours = True, True
                if '2c' in file_ or '2d' in file_:
                    aspect = 2.5
                    line_density = 0.7
                elif '2a' in file_:
                    aspect = 0.7
                    line_density = 0.9
                elif 'mfe' in file_:
                    aspect = 1.75
                    line_density = 1.0
                else:
                    aspect = 0.5
                    line_density = 1.6
                # save 2D plot in model's figures directory
                figname  = figsdir+paper+'_'+var_field+'.eps'
                atm.make_2d_plot(ds, var_field, figname,
                                 normal=[plane,N_2],
                                 aspect=aspect, lines=lines,
                                 contours=contours,
                                 model=paper, figxz=[5.5,5.6],
                                 line_density=line_density
                                 )
        if ('gas','density') in ds.derived_field_list:
            var_field = 'density'
            oneD_arrays = atm.make_1d_slices(ds, var_field, oneD_arrays)

        if ('gas','thermal_pressure') in ds.derived_field_list:
            var_field = 'thermal_pressure'
            oneD_arrays = atm.make_1d_slices(ds, var_field, oneD_arrays)
            figname  = figsdir+paper+'_'+var_field+'.eps'
            lines, contours = True, True
            atm.make_2d_plot(ds, var_field, figname,
                             normal=[plane,N_2],
                             aspect=aspect, lines=lines,
                             contours=contours,
                             model=paper, figxz=[5.5,5.6],
                             line_density=line_density
                             )
        if ('gas','mag_pressure') in ds.derived_field_list:
            var_field = 'mag_pressure'
            oneD_arrays = atm.make_1d_slices(ds, var_field, oneD_arrays)
            figname  = figsdir+paper+'_'+var_field+'.eps'
            lines, contours = True, True
            atm.make_2d_plot(ds, var_field, figname,
                             normal=[plane,N_2],
                             aspect=aspect, lines=lines,
                             contours=contours,
                             model=paper, figxz=[5.5,5.6],
                             line_density=line_density
                             )
    plot_label = figsdir+paper+'_axis.eps'
#    keys = ['alfven_speed','sound_speed','mag_field_z_bg']
    keys = ['thermal_pressure','mag_pressure','density','temperature']
    subkeys = ['axis']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, empirical=True
                                           )
    plot_label = figsdir+paper+'_edge.eps'
    keys = ['thermal_pressure','mag_pressure','density','temperature']
    subkeys = ['edge']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, empirical=True
                                           )
    plot_label = figsdir+paper+'_speeds.eps'
    keys = ['alfven_speed','sound_speed']
    subkeys = ['mean','min','max']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, loc_legend='lower right'
                                           )
    plot_label = figsdir+paper+'_meanz.eps'
    keys = ['thermal_pressure','mag_pressure','density']
    subkeys = ['mean','min','max']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False
                                           )
    plot_label = figsdir+paper+'_beta.eps'
    keys = ['plasma_beta','mag_pressure','thermal_pressure']
    subkeys = ['mean','min','max']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False
                                           )
