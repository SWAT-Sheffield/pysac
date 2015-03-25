# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:13:23 2015

@author: sm1fg

Use to plot 1D and 2D slices of gdf formatted data. 
Dictionaries are used so appropriate labels and line colors can be extracted 
from the data name lists.
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
import os
import warnings
import astropy.units as u
from mayavi_cust_streamlines import sStreamline
from mayavi import mlab
#match the field name to the appropriate axis/colorbar labels. Ordered to 
#control option on, for example, 'mag_pressure'.
#add to the dictionary if variable not specified
y_axis_labels = od((
                ('temp',r'$T$ [K], '),
                ('dens',r'$\rho$ [kg km$^{-3}$], '),
                ('press',r'$p$ [Pa], '),
                ('energy',r'$e$ [N m$^{-2}$], '),
                ('mag_pressure',r''),
                ('mag_field_x',r'$B_x$ [T], '),
                ('mag_field_y',r'$B_y$ [T], '),
                ('mag_field_z',r'$B_z$ [T], '),
                ('mag_strength',r'$|B|$ [T], '),
                ('tension',r'$B\cdot\nabla B/\mu_0$ [Pa m$^{-1}$], '),
                ('beta',r'$\beta$, '),
                ('alfven',r'$v_\mathrm{A}$ [km s$^{-1}$], '),
                ('sound',r'$c_\mathrm{s}$ [km s$^{-1}$], ')))

#match the field name to the appropriate line color for consitency between 
#plots. Ordered to control option on, for example, 'mag_pressure'.
#add to the dictionary if variable not specified
line_colours = od((('temp','red'), ('dens','purple'),
                ('press','green'),
                ('mag_field_x','blue'),
                ('mag_field_y','blue'),
                ('mag_field_z','blue'),
                ('mag_pressure','blue'),
                ('beta','purple'),('alfven','cyan'),('sound','red')))
#determine line styles and widths between mean profiles and outliers
line_styles = {
               'min':':',
               'max':':'
               }
line_standard = ['-','--','-.',':']
line_widths = {'mean':2.0,
               'min':1.0,
               'max':1.0,
               'axis':2.0,
               'edge':2.0
              }

#data may be required from more than one data file for 1D slices, so collate
#1D slices for use together in plotting. 
#include min, max, mean, axial and edge profiles.
def make_1d_slices(ds, var_label, oneD_arrays):
    """extract 1D arrays from gdf data for plotting.
    ds: gdf data set, var_label: gdf field name, oneD_arrays: collection of
    1D slices to which this field's 1D slices shall be appended
    """
    var_ = ds.index.grids[0][var_label]

    #append new 1D slice to existing
    oneD_arrays.update({var_label:{}})
    #vertical profile along the axis
    oneD_arrays[var_label]['axis'] = var_[ds.domain_dimensions[0]/2,
                                          ds.domain_dimensions[1]/2,:]
    #vertical profile along the edge of the domain, useful with axial symmetry
    oneD_arrays[var_label]['edge'] = var_[0,        0        ,:]
    oneD_arrays[var_label]['mean'] = np.mean(np.mean(var_, axis=0), axis=0)
    oneD_arrays[var_label]['min' ] = np.min(np.min(var_, axis=0), axis=0)
    oneD_arrays[var_label]['max' ] = np.max(np.max(var_, axis=0), axis=0)
    #integer size and range of vertical domain may vary so match for each array
    oneD_arrays[var_label]['Z'   ] = np.linspace(ds.domain_left_edge[2],
                                                 ds.domain_right_edge[2],
                                                 ds.domain_dimensions[2])
    return oneD_arrays

def make_1d_zplot(f, plot_label,
                 keys = ['pressure_HS','density_HS','temperature'],
                 subkeys = ['axis'],
                 figxy=[6.47,4.0],
                 ylog = True, xlog = False, loc_legend='center right',
                 empirical = False
                ):
    """select 1D arrays and plot appropriate slices['y',nx_2]
    f: the collated set of labelled 1D arrays
    keys: variable's field name
    subkeys: label discriminating which slice of variable to be plotted with Z
    figxy: size of figure
    ylog/xlog: True to use logarithmic scale
    """
    plt.figure(figsize = figxy)
    if empirical:
        import pysac.mhs_atmosphere as atm
        table = atm.read_VAL3c_MTW()
        table_keys = ['p','rho','T']
        table_sym = ['g*','m+','rx']
        table_units = ['Pa','kg/km3','K']
        Z = u.Quantity(table['Z'], unit=table['Z'].unit).to('Mm')
        for tab in range(0,3):
            empiric = u.Quantity(table[table_keys[tab]], 
                                 unit = table[table_keys[tab]].unit
                                ).to(table_units[tab])
            plt.plot(Z, empiric, table_sym[tab])
    #initialise set of tags to be plotted and count for line styles
    tags = set()
    count = 0
    for key in keys:
        #select line style for field if not overwritten by subkey
        linestyle = line_standard[count]
        if count > 3:
            count -= 4
        for subkey in subkeys:
            count += 1
            color = None
            for tag in line_colours.keys():
                if tag in key:
                    tags.add(tag)
                    color = line_colours[tag]
            #label field only when required for legend
            if subkey in 'min' or subkey in'max':
                nolabel = True
            else:
                nolabel = False
                label = str.replace(subkey+' '+key, '_',' ')
            for tag in line_styles.keys():
                #replace line style and wind back count for next field
                #if line_style to be used instead
                if tag in subkey:
                    linestyle = line_styles[tag]
                    count -= 1
            for tag in line_widths.keys():
                if tag in subkey:
                    linewidth = line_widths[tag]
            #rescale density so camparable in plot with pressure/temperature
            rescale = np.string_(f[key][subkey].units)
            if 'density' in key:
                rescale = 'kg / km**3'
            if nolabel: 
                plt.plot(f[key]['Z'].in_units('Mm'),
                         f[key][subkey].in_units(rescale),color=color,
                         linestyle=linestyle, lw=linewidth)
            else:
                plt.plot(f[key]['Z'].in_units('Mm'),
                         f[key][subkey].in_units(rescale),color=color,
                         linestyle=linestyle, lw=linewidth, label=label)
    #collate labels for y-axis and apply log scale if True
    ylabel = ''
    for tag in tags:
        ylabel += y_axis_labels[tag]
    if ylog:
        plt.gca().set_yscale('log',subsy=[5,10])
    if xlog:
        plt.gca().set_xscale('log',subsy=[5,10])
#    plt.gca().set_yticks([1.0e-3,1.0e-1,1.0e1,1.0e3,1.0e5])
    plt.xlabel('Height [Mm]')
    if empirical:
        #limit x axis to simulation data rather than empirical range
        xmax = f[keys[0]]['Z'].in_units('Mm').max()
        xmax /= xmax.unit_quantity
        plt.xlim(0,xmax)
    #tidy label string by stripping final ', ' from the end  
    plt.ylabel(ylabel[:-2])
    #raise plot x-axis to fit x-label 
    plt.subplots_adjust(bottom=0.125)
    #consider moving legend for different plots, add a loc to function call?
    legend = plt.legend(loc=loc_legend)
    plt.savefig(plot_label)

#-----------------------------------------------------------------------------
def make_2d_plot(ds, var_field, figname, normal = ['y',64],
                   lines=False, contours=False,
                   figxz=[6.5,5.60], aspect=1.0, model='', 
                   line_density = 0.75
                  ):
    """ds: a gdf data structure, from which the variable labelled by the
    string var_field is selected. 
    figname: name derived from the model name, variable name and directory 
    into which the plot will be saved.
    normal: specifies orientation of plane and which slice from 0 to N-1) 
    options allow magnetic field lines to be overplotted and/or plasma beta
    contours.
    """
    #convert density to units consistent with 1d plots and label
    if 'dens' in var_field:
        var = ds.index.grids[0][var_field].in_units('kg / km**3')
    else:
        var = ds.index.grids[0][var_field]

    #select the 2D slice from the 3D array (normal[0]: direction of plane
    #and normal[1] the integer index of the slice selected)
    #extent: range of plot domain
    if normal[0] is 'x':
        slc = var[normal[1],:,:]/var.unit_quantity
        extent=[ds.domain_left_edge[1].in_units('Mm'),
                ds.domain_right_edge[1].in_units('Mm'),
                ds.domain_left_edge[2].in_units('Mm'),
                ds.domain_right_edge[2].in_units('Mm')]
    if normal[0] is 'y':
        slc = var[:,normal[1],:]/var.unit_quantity
        extent=[ds.domain_left_edge[0].in_units('Mm'),
                ds.domain_right_edge[0].in_units('Mm'),
                ds.domain_left_edge[2].in_units('Mm'),
                ds.domain_right_edge[2].in_units('Mm')]
    if normal[0] is 'z':
        slc = var[:,:,normal[1]]/var.unit_quantity
        extent=[ds.domain_left_edge[0].in_units('Mm'),
                ds.domain_right_edge[0].in_units('Mm'),
                ds.domain_left_edge[1].in_units('Mm'),
                ds.domain_right_edge[1].in_units('Mm')]

    #select colour tables according to the variable being plotted
    colour = cm.YlOrBr
    if 'mag_field' in var_field:
        colour = cm.BuPu
    if 'temperature' in var_field:
        colour = cm.RdBu_r
    if 'beta' in var_field:
        colour = cm.PiYG_r
    if 'tension' in var_field:
        colour = cm.RdYlGn_r
    if 'velocity' in var_field:
        colour = cm.RdYlGn
    if 'alfven' in var_field:
        colour = cm.Blues
    if 'sound' in var_field:
        colour = cm.Greens_r

    #increase the horizontal extent to include second colorbar for field lines
    if lines:
        xz2 = figxz[0]+2
        figxz = [xz2,figxz[1]]
    plt.figure(figsize=figxz)

    #adjust color scaling appropriate for field profiles
    if slc.min() < 0:
        norm = colors.SymLogNorm(slc.max())
        l_norm = True
    elif 'tension' in var_field:
        norm = colors.Normalize()
        l_norm = False
    else:
        if slc.max()/slc.min() > 5e2:
            norm = colors.LogNorm()
            l_norm = True
        else:
            norm = colors.Normalize()
            if slc.max() > 1e3:
                l_norm = True
            else:
                l_norm = False

    #plot a 2D slice
    plt.imshow(slc.T, colour,
               extent=extent,
               norm=norm,
               origin='lower',
               vmin=slc.min(),
               vmax=slc.max()
              )
    plt.axes().set_aspect(aspect=aspect)
    plt.ylabel('Height [Mm]')
    plt.xlabel('Width [Mm]')
    plt.autoscale(False) # fix plot to extent

    # add field's colorbar with appropriate label
    if l_norm:
        l_f = LogFormatterMathtext(10, labelOnlyBase=False)
    else:
        l_f = FormatStrFormatter('%.0f')
    ylabel = ''
    for key in y_axis_labels:
        if key in var_field:
            ylabel += y_axis_labels[key][:-2]
    cbar1 = plt.colorbar(format = l_f,
                       #ticks=[var.min(),1,10,100,1000,var.max()],
                        fraction=0.2)
    cbar1.ax.set_ylabel(ylabel)
    cbar1.solids.set_edgecolor("face")

    # obtain data arrays for beta contours and magnetic field lines and save
    # for use in plots for other variables - magnetic field must data must be
    # included in gdf file of first run.
    if lines or contours:
        warnings.warn("magnetic field data required for first run to obtain"+
                      "field lines or Beta contours", Warning)
        #set up files where magnetic and beta slices will be stored
        if len(model) == 0:
            raise ValueError("in plot.mhs_plot.make_2d_plot specify the model")
        flinesdir = os.path.expanduser('~/Documents/mhs_atmosphere/flines/')
        if not os.path.exists(flinesdir):
            os.makedirs(flinesdir)
        v1_file = flinesdir+'v1_norm_'+normal[0]+'_slice_'+\
                           str(normal[1])+'_'+model+'.npy'
        v2_file = flinesdir+'v2_norm_'+normal[0]+'_slice_'+\
                           str(normal[1])+'_'+model+'.npy'
        beta_file = flinesdir+'beta_norm_'+normal[0]+'_slice_'+\
                           str(normal[1])+'_'+model+'.npy'
        #acquire the B field components and beta for the lines and contours
        if ('gdf','mag_field_x_bg') not in ds.field_list:
            v1 = np.load(v1_file)
            v2 = np.load(v2_file)
            beta = np.load(beta_file)
        else:
            if normal[0] is 'x':
                v1 = ds.index.grids[0]['mag_field_y_bg'][normal[1],:,:]
                v2 = ds.index.grids[0]['mag_field_z_bg'][normal[1],:,:]
                beta =  ds.index.grids[0]['plasma_beta'][normal[1],:,:]
            if normal[0] is 'y':                      
                v1 = ds.index.grids[0]['mag_field_x_bg'][:,normal[1],:]
                v2 = ds.index.grids[0]['mag_field_z_bg'][:,normal[1],:]
                beta =  ds.index.grids[0]['plasma_beta'][:,normal[1],:]
            if normal[0] is 'z':                      
                v1 = ds.index.grids[0]['mag_field_x_bg'][:,:,normal[1]]
                v2 = ds.index.grids[0]['mag_field_y_bg'][:,:,normal[1]]
                beta =  ds.index.grids[0]['plasma_beta'][:,:,normal[1]]
            v1 = np.array(v1)
            v2 = np.array(v2)
            beta = np.array(beta)
            #save B and beta for use with other gdf data sets
            np.save(v1_file, v1)
            np.save(v2_file, v2)
            np.save(beta_file, beta)
        #plot field lines as streamlines with color and thickness functions of
        #field magnitude
        """Consider replacing strength from 2D arrays with actual derived field
        of 'mag_magnitude' saved as fourth data set speed?
        Also include option to apply velocity and alternative to beta: FAG
        """
        if lines:
            clabel = y_axis_labels['mag_strength'][:-2]
            X = np.linspace(extent[0],extent[1],slc.shape[0])       
            Y = np.linspace(extent[2],extent[3],slc.shape[1])
            U = np.sqrt(v1**2 + v2**2)
            plt.streamplot(X.value, Y.value, v1.T, v2.T, color=U.T,
                           cmap=cm.BuPu, density=line_density, 
                           linewidth = 3*U.T/U.max(), minlength = 0.99
                          ) 
            #add second colorbar for field strength 
            cbar2 = plt.colorbar(format = l_f,
                               #ticks=[var.min(),1,10,100,1000,var.max()],
                                fraction=0.2)
            cbar2.ax.set_ylabel(clabel)
            cbar2.solids.set_edgecolor("face")
        # plot contours with inline labels 
        if contours:
            X = np.linspace(extent[0],extent[1],slc.shape[0])       
            Y = np.linspace(extent[2],extent[3],slc.shape[1])
            CS=plt.contour(X,Y,beta.T,4,
                           levels=[beta.min(),
                                   1e-2,1e-1,1.,1e1,1e2,
                                   beta.max()],
                           linewidths=1,cmap=cm.PiYG_r,norm=LogNorm())
            if beta.min() < 1e-2:
                fomt = '%.3f'
            elif beta.min() < 1e-1:
                fomt = '%.2f'
            elif beta.min() < 1:
                fomt = '%.1f'
            else:
                fomt = '%.0f'
            plt.clabel(CS, fmt=fomt)
    plt.savefig(figname)

##============================================================================
## Fieldline Generation
##============================================================================
#def mad_3d_plot(ds, fields, figname,
#                figxy=[900,950]
#               ):
#    
#    cube_slice = np.s_[:,:,:]
#    x_slice = np.s_[:,:,:,:]
#    scene = mlab.figure(1, bgcolor=(1, 1, 1),
#                    fgcolor=(0.5, 0.5, 0.5),size=figxy)
#
#    mlab.savefig(figname)
#    import pdb; pdb.set_trace()

##============================================================================
## Fieldline Generation
##============================================================================
##def get_seeds(model):
"""Depricated version to calculate field lines, subsequent to inbuilt 
streamlines with color and thickness analogue
"""
#    
#
#
#def make_field_lines(v1,v2,seeds,extent,flines_file,dS=1):
#    """ Simple Euler fieldline integrator
#    v1,v2 - y and x vectors*lscale
#    seeds - array of coordinate (array indexes) pairs
#    dS [optional] - step size
#    """
#    from scipy import ndimage
#    
#    field = []
#    #Get extent of domain
#    max1 = v1.shape[0]
#    max2 = v2.shape[1]
#    min2 = 0
#    min1 = 0
#    vmax = np.sqrt(np.max(v1**2+v2**2))/v1.unit_quantity
#    #integrate for each fieldline
#    for seed1 in seeds[1]:
#        for seed0 in seeds[0]:
#            c1,c2 = seed0,seed1
#            out1 = [c1] #first point is seed
#            out2 = [c2]
#            coords = np.array([[c1],[c2]])
#            v1i = ndimage.map_coordinates(v1,coords)[0]
#            v2i = ndimage.map_coordinates(v2,coords)[0]
#            out3 = [np.sqrt(v1i**2 + v2i**2)/vmax]
#            cnt = 0
#            while (c1 <= max1-0.5 and c1 >= min1+0.5) and (c2 <= max2-0.5
#                                                    and c2 >= min2+0.25):
#                #Interpolate the vector field to the coord downwards
#                coords = np.array([[c1],[c2]])
#                v1i = ndimage.map_coordinates(v1,coords)[0]
#                v2i = ndimage.map_coordinates(v2,coords)[0]
#                vi = np.sqrt(v1i**2 + v2i**2)
#                d1 = ( v1i * dS )/ vi
#                d2 = ( v2i * dS )/ vi
#                c1 -= d1
#                c2 -= d2
#                c3 = vi/vmax
#                out1.append(c1)
#                out2.append(c2)
#                out3.append(c3)
#                cnt += 1
#                if cnt > max2*5:
#                    print "limit seeds = ",seed0,seed1
#                    break
##            out = np.zeros([len(out1),len(out2),len(out3)])
##            out[0] = out1 * (extent[1] - extent[0])/(max1-1) + extent[0]
##            out[1] = out2 * (extent[3] - extent[2])/(max2-1) + extent[2]
##            out[2] = out3
##            field.append(out)
#            c1,c2 = seed0,seed1
#            cnt = 0
#            while (c1 <= max1-0.5 and c1 >= min1+0.5) and (c2 <= max2-0.5
#                                                    and c2 >= min2+0.25):
#                #Interpolate the vector field to the coord upwards
#                coords = np.array([[c1],[c2]])
#                v1i = ndimage.map_coordinates(v1,coords)[0]
#                v2i = ndimage.map_coordinates(v2,coords)[0]
#                vi = np.sqrt(v1i**2 + v2i**2)
#                d1 = ( v1i * dS )/ vi
#                d2 = ( v2i * dS )/ vi
#                c1 += d1
#                c2 += d2
#                c3 = vi/vmax
#                out1.append(c1)
#                out2.append(c2)
#                out3.append(c3)
#                cnt += 1
#                if cnt > max2*5:
#                    print "limit seeds = ",seed0,seed1
#                    break
#            out = np.zeros([3,len(out1)])
#            out[0] = out1 * (extent[1] - extent[0])/(max1-1) + extent[0]
#            out[1] = out2 * (extent[3] - extent[2])/(max2-1) + extent[2]
#            out[2] = out3
#            field.append(out)
##    import pdb; pdb.set_trace()
#    np.save(flines_file, field)        
#
#    return np.array(field)

#============================================================================
# Fieldline Generation
#============================================================================
