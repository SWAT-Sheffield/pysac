# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: stuart
"""
import os
import glob
import numpy as np
import pysac.yt as sacyt
import pysac.mhs_atmosphere as atm

l_mpi=False
scales, physical_constants = \
    atm.units_const.get_parameters()
#define the models required
spruits = ['spruit_const','spruit_sqrt','spruit_linear','spruit_square']
#spruits = ['spruit_const']
oneD_arrays = {}
oned_dataset = []
#loop over all four models
for spruit in spruits:
    datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+
                                 spruit+'/')
    figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+spruit+'/')
    if not os.path.exists(figsdir):
        os.makedirs(figsdir)
    #open all gdf files in the model directory
    files = glob.glob(datadir+'/*')
    #files = glob.glob(datadir+'/'+spruits[0]+'_3Daux.gdf')
    files.sort()

    print(files)

    for file_ in files:
        #ds = yt.load(file_)
        ds = sacyt.SACGDFDataset(file_)

        vars_ = ds.index.field_list
        for var_ in vars_:
            var_field = var_[1]
            max_var = np.max(np.abs(ds.index.grids[0][var_field]))/\
                             ds.index.grids[0][var_field].unit_quantity
            var = ds.index.grids[0][var_field]
            if max_var > 0.:
                # save 1D slices from each variable for plotting
                oneD_arrays = atm.mhs_plot.make_1d_slices(ds, var_field, oneD_arrays)
                # select the central slice to plot normal to the y-plane
                nx_2 = ds.domain_dimensions[1]/2
                if 'mag' in var_field:
                    lines = True
                elif 'density' in var_field or 'pressure' in var_field:
                    lines = True
                    if 'D' in file_:
                        lines = False
                else:
                    lines = False
                if '_HS' in var_field:
                    contours = False
                else:
                    contours = True
                # save 2D plot in model's figures directory
                figname  = figsdir+spruit+'_'+var_field+'.eps'
                atm.mhs_plot.make_2d_plot(ds, var_field, figname,
                                                    normal=['y',nx_2],
                                                    aspect=0.2, lines=lines,
                                                    contours=contours,
                                                    model=spruit)
        if ('gas','thermal_pressure') in ds.derived_field_list:
            var_field = 'thermal_pressure'
            figname  = figsdir+spruit+'_'+var_field+'.eps'
            lines, contours = True, True
            atm.mhs_plot.make_2d_plot(ds, var_field, figname,
                                                normal=['y',nx_2],
                                                aspect=0.2, lines=lines,
                                                contours=contours,
                                                model=spruit)
        if ('gas','mag_pressure') in ds.derived_field_list:
            var_field = 'mag_pressure'
            figname  = figsdir+spruit+'_'+var_field+'.eps'
            lines, contours = True, True
            atm.mhs_plot.make_2d_plot(ds, var_field, figname,
                                                normal=['y',nx_2],
                                                aspect=0.2, lines=lines,
                                                contours=contours,
                                                model=spruit)

    plot_label = figsdir+spruit+'_axis.eps'
    keys = ['alfven_speed','sound_speed','mag_field_z_bg']
    subkeys = ['axis']
    atm.mhs_plot.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = True
                                           )
