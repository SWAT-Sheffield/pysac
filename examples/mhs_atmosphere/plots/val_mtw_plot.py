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
#papers = ['paper1','paper2a','paper2b','paper2c','paper2d','mfe_setup']
#papers = ['mfe_setup']
papers = ['paper2b']
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
#    files = glob.glob(datadir+'/'+papers[0]+'_3Daux.gdf')
#    files = glob.glob(datadir+'/'+papers[0]+'.gdf')
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
#                lines, contours = True, True
                if '1D' in file_:
                    lines, contours = False, False
                else:
                    if '_HS' in var_field:
                        lines, contours = False, False
                    elif 'tension' in var_field or 'balancing' in var_field:
                        lines, contours = True, False
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
                    line_density = 1.9
                else:
                    aspect = 0.5
                    line_density = 1.1
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
    ymin = min(oneD_arrays['density']['axis'].in_units('kg / km**3').value.min(),
               oneD_arrays['temperature']['edge'].value.min(),
               oneD_arrays['temperature']['axis'].value.min(),
               oneD_arrays['density']['edge'].in_units('kg / km**3').value.min(),
               oneD_arrays['thermal_pressure']['edge'].value.min(),
               oneD_arrays['thermal_pressure']['axis'].value.min())*0.5
    ymax = max(oneD_arrays['density']['axis'].in_units('kg / km**3').value.max(),
               oneD_arrays['temperature']['edge'].value.max(),
               oneD_arrays['temperature']['axis'].value.max(),
               oneD_arrays['density']['edge'].in_units('kg / km**3').value.max(),
               oneD_arrays['thermal_pressure']['edge'].value.max(),
               oneD_arrays['thermal_pressure']['axis'].value.max())*2.
    if 'mfe_setup' in paper:
        loc_legend='lower left'
    else:
        loc_legend='center right'
    keys = ['thermal_pressure','mag_pressure','density','temperature']
    subkeys = ['axis']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, empirical=True,
                      loc_legend=loc_legend, ylim = (0.9*ymin,5*ymax)
                                           )
    plot_label = figsdir+paper+'_edge.eps'
    keys = ['thermal_pressure','mag_pressure','density','temperature']
    subkeys = ['edge']
    atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, empirical=True,
                      loc_legend=loc_legend, ylim = (0.9*ymin,5*ymax)
                                           )
    plot_label = figsdir+paper+'_speeds.eps'
    keys = ['alfven_speed','sound_speed']
    subkeys = ['mean','min','max']
    if 'mfe_setup' in paper:
        atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                          ylog = True, xlog = False, loc_legend='lower left',
                ylim = (2e3,oneD_arrays['sound_speed']['max'].value.max())
                                           )
    else:
        atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, loc_legend='lower right'
                                           )
    plot_label = figsdir+paper+'_meanz.eps'
    keys = ['thermal_pressure','mag_pressure','density']
    subkeys = ['mean','min','max']
    if 'mfe_setup' in paper:
        ymax = oneD_arrays['density']['max'].in_units('kg / km**3').value.max()
        atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, loc_legend='upper right',
                ylim = (1e-2,1.1*ymax)
                                           )
    else:
        atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, loc_legend='upper right'
                                           )
    plot_label = figsdir+paper+'_beta.eps'
    keys = ['plasma_beta','mag_pressure','thermal_pressure']
    subkeys = ['mean','min','max']
    if 'mfe_setup' in paper:
        atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, loc_legend='upper right',
                      ylim = (1e-3,1.1*oneD_arrays['plasma_beta']['max'].max())
                                           )
    else:
        atm.make_1d_zplot(oneD_arrays, plot_label, keys=keys, subkeys=subkeys,
                      ylog = True, xlog = False, loc_legend='upper right'
                                           )
    plot_label = figsdir+paper+'_compare.eps'
    mfe_Bz = np.load('mfe_Bz.npy')
    mfe_Z  = np.load('mfe_zz.npy')
    import matplotlib.pyplot as plt
    plt.figure(figsize=[6.47,4.0])
    if 'mfe_setup' in paper:
        plt.plot(mfe_Z, mfe_Bz/1000.,
                 'b-', label=r"MFE '15 $B_z$(axis)"
                )
        plt.plot(oneD_arrays['mag_field_z_bg']['Z'].in_units('Mm'),
                 oneD_arrays['mag_field_z_bg']['axis'].in_units('T'), 'm-.',
                 label=r"GFE '15 $B_z$(axis)", lw=2.0
                )
        plt.gca().set_yscale('log',subsy=[5,10])
        plt.xlabel('Height [Mm]')
        plt.ylabel(r'$B_z$ [T]')
        plt.legend(loc='lower left')
        plt.ylim(4e-4,1.25e-1)
        plt.subplots_adjust(bottom=0.125)
        plt.savefig(plot_label)
    else:
        plt.plot(mfe_Z, mfe_Bz/1000.,
                 'b-', label=r"MFE '15 $B_z$(axis)"
                )
        plt.plot(oneD_arrays['mag_field_z_bg']['Z'].in_units('Mm'),
                 oneD_arrays['mag_field_z_bg']['axis'].in_units('T'), 'm-.',
                 label=r"GFE '15 $B_z$(axis)", lw=2.0
                )
        plt.gca().set_yscale('log',subsy=[5,10])
        plt.xlabel('Height [Mm]')
        plt.ylabel(r'$B_z$ [T]')
        plt.legend(loc='lower left')
        ymin = oneD_arrays['mag_field_z_bg']['axis'].in_units('T').value.min()
        plt.ylim(0.05*ymin,1.25e-1)
        plt.subplots_adjust(bottom=0.125)
    plt.savefig(plot_label)

    if 'paper1' in paper:
        plt.figure(figsize=[6.47,4.0])
        import astropy.units as u
        from pysac.mhs_atmosphere.parameters.model_pars import paper1 as model_pars
        coords = atm.get_coords(model_pars['Nxyz'], u.Quantity(model_pars['xyz']))

        empirical_data = atm.read_VAL3c_MTW(mu=physical_constants['mu'])

        table = atm.interpolate_atmosphere(empirical_data,
                                           coords['Zext']
                                          )
        rho = u.Quantity(table['rho'], copy=True).to('kg km-3').value                                          
        Z = u.Quantity(table['Z'], copy=True).to('Mm').value                                          
        p = u.Quantity(table['p'], copy=True).to('Pa').value                                          
        plt.plot(Z,rho,'m-.',lw=2.0,label=r"density")
        plt.plot(Z,table['T'],'r:',lw=2.0,label=r"temperature")
        plt.plot(Z,p,'g-',lw=2.0,label=r"pressure") 
        plt.gca().set_yscale('log',subsy=[5,10])
        plt.xlabel('Height [Mm]')
        plt.ylabel(r'$p$ [Pa], $T$ [K], $\rho$ [kg km$^{-3}$]')
        plt.legend(loc='center right')
        plt.subplots_adjust(bottom=0.125)
        plt.xlim(0,8.8)
        plt.savefig(figsdir+'data_match.eps')
