# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:13:23 2015

@author: sm1fg
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import matplotlib.colorbar as cb
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.ticker import LogFormatterMathtext, FormatStrFormatter
import matplotlib.colors as colors
from collections import OrderedDict as od
y_axis_labels = od((
                ('temp',r'$T$ [K], '),
                ('dens',r'$\rho$ [kg km$^{-3}$], '),
                ('press',r'$p$ [Pa], '),
                ('mag_field_x',r'$B_x$ [T], '),
                ('mag_field_y',r'$B_y$ [T], '),
                ('mag_field_z',r'$B_z$ [T], '),
                ('tension',r'$B\cdot\nabla B/\mu_0$ [Pa m$^{-1}$], '),
                ('beta',r'$\beta$, '),
                ('alfven',r'$v_\mathrm{A}$ [km s$^{-1}$], '),
                ('sound',r'$c_\mathrm{s}$ [km s$^{-1}$], ')))
line_colours = od((('temp','red'), ('dens','purple'),
                ('mag_field_x','blue'),
                ('mag_field_y','blue'),
                ('mag_field_z','blue'),
                ('press','green'),
                ('beta','purple'),('alfven','cyan'),('sound','red')))

line_styles = {
               'min':':',
               'max':':',
               }
line_standard = ['-','-.','--',':']
line_widths = {'mean':3.0,
               'min':1.0,
               'max':1.0,
               'axis':3.0,
               'edge':3.0
              }

def make_1d_slices(ds, var_label, oneD_arrays):
    """extract 1D arrays from gdf data for plotting
    """
    var_ = ds.index.grids[0][var_label]

#    oneD_arrays.create[var_label]
    oneD_arrays.update({var_label:{}})
    oneD_arrays[var_label]['axis'] = var_[ds.domain_dimensions[0]/2,
                                          ds.domain_dimensions[1]/2,:]
    oneD_arrays[var_label]['edge'] = var_[0,        0        ,:]
    oneD_arrays[var_label]['mean'] = np.mean(np.mean(var_, axis=0), axis=0)
    oneD_arrays[var_label]['min' ] = np.min(np.min(var_, axis=0), axis=0)
    oneD_arrays[var_label]['max' ] = np.max(np.max(var_, axis=0), axis=0)
    oneD_arrays[var_label]['Z'   ] = np.linspace(ds.domain_left_edge[2],
                                                 ds.domain_right_edge[2],
                                                 ds.domain_dimensions[2])

    return oneD_arrays

def make_1d_zplot(f, plot_label,
                 keys = ['pressure_HS','density_HS','temperature'],
                 subkeys = ['axis'],
                 figxy=[6.47,4.0],
                 ylog = True, xlog = False
                ):
    """select 1D arrays and plot appropriate slices
    """
    plt.figure(figsize = figxy)
    tags = set()
    count = 0
    for key in keys:
        linestyle = line_standard[count]
        count +=1
        if count > 3:
            count -= 4
        for subkey in subkeys:
            color = None
            for tag in line_colours.keys():
                if tag in key:
                    tags.add(tag)
                    color = line_colours[tag]
            label = str.replace(subkey+' '+key, '_',' ')
            print(tags)
            for tag in line_styles.keys():
                if tag in subkeys:
                    linestyle = line_styles[tag]
            for tag in line_widths.keys():
                if tag in subkeys:
                    linewidth = line_widths[tag]
            rescale = np.string_(f[key][subkey].units)
            if 'density' in key:
                rescale = 'kg / km**3'
            plt.plot(f[key]['Z'].in_units('Mm'),
                     f[key][subkey].in_units(rescale),color=color,
                     linestyle=linestyle, lw=linewidth, label=label)
    ylabel = ''
    for tag in tags:
        ylabel += y_axis_labels[tag]
    if ylog:
        plt.gca().set_yscale('log',subsy=[5,10])
    if xlog:
        plt.gca().set_xscale('log',subsy=[5,10])
#    plt.gca().set_yticks([1.0e-3,1.0e-1,1.0e1,1.0e3,1.0e5])
    plt.xlabel('Height [Mm]')
    plt.ylabel(ylabel[:-2])
    plt.subplots_adjust(bottom=0.125)
#    import pdb; pdb.set_trace()
#    handles, labels = ax.Axes.get_legend_handles_labels()
    legend = plt.legend(loc='center right')
    plt.savefig(plot_label)

#-----------------------------------------------------------------------------
def make_2d_plot(ds, var_field, figname, normal = ['y',64],
                   lines=False, figxz=[4.5,5.60], aspect=1.0
                  ):
    if 'dens' in var_field:
        var = ds.index.grids[0][var_field].in_units('kg / km**3')
    else:
        var = ds.index.grids[0][var_field]
    if normal[0] is 'x':
        slc = var[normal[1],:,:]/var.unit_quantity
    if normal[0] is 'y':
        slc = var[:,normal[1],:]/var.unit_quantity
    if normal[0] is 'z':
        slc = var[:,:,normal[1]]/var.unit_quantity

    colour = cm.YlOrBr
    if 'mag_field' in var_field:
        colour = cm.BuPu
    if 'temperature' in var_field:
        colour = cm.RdBu_r
    if 'beta' in var_field:
        colour = cm.PiYG
    if 'tension' in var_field:
        colour = 'RdYlGn_r'
#    slc.set_cmap(var_field, colour)
 #   slc.plots[var_field].hide_colorbar()
#    slc.set_cbar_minorticks('all', 'off')
    plt.figure(figsize=figxz)
#    plt.figure()
    #lpmin=np.log10(var.min())
    #lpmax=np.log10(var.max())
    #colarr=10**(np.linspace(lpmin,lpmax,nz*2))
    if slc.min() < 0:
        norm = colors.SymLogNorm(slc.max())
        l_norm = True
    elif 'tension' in var_field:
        norm = colors.Normalize()
        l_norm = False
    else:
        if slc.max()/slc.min() > 1e2:
            norm = colors.LogNorm()
            l_norm = True
        else:
            norm = colors.Normalize()
            l_norm = False
    plt.imshow(slc.T, colour,
    extent=[ds.domain_left_edge[0].in_units('Mm'),
           ds.domain_right_edge[0].in_units('Mm'),
           ds.domain_left_edge[2].in_units('Mm'),
           ds.domain_right_edge[2].in_units('Mm')
           ],
    norm=norm,
    origin='lower',
    vmin=slc.min(),
    vmax=slc.max())
    plt.axes().set_aspect(aspect=aspect)
    plt.ylabel('Height [Mm]')
    plt.xlabel('Width [Mm]')
    if l_norm:
        l_f = LogFormatterMathtext(10, labelOnlyBase=False)
    else:
        l_f = FormatStrFormatter('%.0f')
    ylabel = ''
    for key in y_axis_labels:
        if key in var_field:
            ylabel += y_axis_labels[key][:-2]
    cbar = plt.colorbar(format = l_f,
                       #ticks=[var.min(),1,10,100,1000,var.max()],
                        fraction=0.2)
    cbar.ax.set_ylabel(ylabel)
    cbar.solids.set_edgecolor("face")
#    import pdb; pdb.set_trace()
#    slc.save(figsdir+figname)
    plt.savefig(figname)
