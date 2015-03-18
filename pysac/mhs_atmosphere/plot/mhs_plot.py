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
import os
import warnings
y_axis_labels = od((
                ('temp',r'$T$ [K], '),
                ('dens',r'$\rho$ [kg km$^{-3}$], '),
                ('press',r'$p$ [Pa], '),
                ('mag_field_x',r'$B_x$ [T], '),
                ('mag_field_y',r'$B_y$ [T], '),
                ('mag_field_z',r'$B_z$ [T], '),
                ('mag_strength',r'$|B|$ [T], '),
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
    """select 1D arrays and plot appropriate slices['y',nx_2]
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
                   lines=False, contours=False,
                   figxz=[6.5,5.60], aspect=1.0, model=''
                  ):
    """ds: a gdf data structure, from which the variable labelled by the
    string var_field is selected. 
    figname: name derived from the model name, variable name and directory 
    into which the plot will be saved.
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

    #increase the horizontal extent to include second colorbar for field lines
    if lines:
        xz2 = figxz[0]+2
        figxz = [xz2,figxz[1]]
    plt.figure(figsize=figxz)

    #adjust color scaling for variable profiles
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

    #plot a 2D 
    plt.imshow(slc.T, colour,
    extent=extent,
    norm=norm,
    origin='lower',
    vmin=slc.min(),
    vmax=slc.max())
    plt.axes().set_aspect(aspect=aspect)
    plt.ylabel('Height [Mm]')
    plt.xlabel('Width [Mm]')
    plt.autoscale(False)
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
    # included in first run.
    if lines or contours:
        warnings.warn("magnetic field data required for first run to obtain"+
                      "field lines or Beta contours", Warning)
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
        if 'mag_field_x_bg' not in np.array(ds.field_list)[:,1]:
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
            np.save(v1_file, v1)
            np.save(v2_file, v2)
            np.save(beta_file, beta)
        if lines:
            clabel = y_axis_labels['mag_strength'][:-2]
            X = np.linspace(extent[0],extent[1],slc.shape[0])       
            Y = np.linspace(extent[2],extent[3],slc.shape[1])
            U = np.sqrt(v1**2 + v2**2)
            plt.streamplot(X.value, Y.value, v1.T, v2.T, color=U.T,
                           cmap=cm.BuPu, density=0.5, 
                           linewidth = 2*np.exp(np.log(U.T)-np.log(U.max()))
                          ) 
            cbar2 = plt.colorbar(format = l_f,
                               #ticks=[var.min(),1,10,100,1000,var.max()],
                                fraction=0.2)
            cbar2.ax.set_ylabel(clabel)
            cbar2.solids.set_edgecolor("face")
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
#            import pdb; pdb.set_trace()
            plt.clabel(CS, fmt=fomt)
#    import pdb; pdb.set_trace()
    plt.savefig(figname)

#============================================================================
# Fieldline Generation
#============================================================================
#def get_seeds(model):
    


def make_field_lines(v1,v2,seeds,extent,flines_file,dS=1):
    """ Simple Euler fieldline integrator
    v1,v2 - y and x vectors*lscale
    seeds - array of coordinate (array indexes) pairs
    dS [optional] - step size
    """
    from scipy import ndimage
    
    field = []
    #Get extent of domain
    max1 = v1.shape[0]
    max2 = v2.shape[1]
    min2 = 0
    min1 = 0
    vmax = np.sqrt(np.max(v1**2+v2**2))/v1.unit_quantity
    #integrate for each fieldline
    for seed1 in seeds[1]:
        for seed0 in seeds[0]:
            c1,c2 = seed0,seed1
            out1 = [c1] #first point is seed
            out2 = [c2]
            coords = np.array([[c1],[c2]])
            v1i = ndimage.map_coordinates(v1,coords)[0]
            v2i = ndimage.map_coordinates(v2,coords)[0]
            out3 = [np.sqrt(v1i**2 + v2i**2)/vmax]
            cnt = 0
            while (c1 <= max1-0.5 and c1 >= min1+0.5) and (c2 <= max2-0.5
                                                    and c2 >= min2+0.25):
                #Interpolate the vector field to the coord downwards
                coords = np.array([[c1],[c2]])
                v1i = ndimage.map_coordinates(v1,coords)[0]
                v2i = ndimage.map_coordinates(v2,coords)[0]
                vi = np.sqrt(v1i**2 + v2i**2)
                d1 = ( v1i * dS )/ vi
                d2 = ( v2i * dS )/ vi
                c1 -= d1
                c2 -= d2
                c3 = vi/vmax
                out1.append(c1)
                out2.append(c2)
                out3.append(c3)
                cnt += 1
                if cnt > max2*5:
                    print "limit seeds = ",seed0,seed1
                    break
#            out = np.zeros([len(out1),len(out2),len(out3)])
#            out[0] = out1 * (extent[1] - extent[0])/(max1-1) + extent[0]
#            out[1] = out2 * (extent[3] - extent[2])/(max2-1) + extent[2]
#            out[2] = out3
#            field.append(out)
            c1,c2 = seed0,seed1
            cnt = 0
            while (c1 <= max1-0.5 and c1 >= min1+0.5) and (c2 <= max2-0.5
                                                    and c2 >= min2+0.25):
                #Interpolate the vector field to the coord upwards
                coords = np.array([[c1],[c2]])
                v1i = ndimage.map_coordinates(v1,coords)[0]
                v2i = ndimage.map_coordinates(v2,coords)[0]
                vi = np.sqrt(v1i**2 + v2i**2)
                d1 = ( v1i * dS )/ vi
                d2 = ( v2i * dS )/ vi
                c1 += d1
                c2 += d2
                c3 = vi/vmax
                out1.append(c1)
                out2.append(c2)
                out3.append(c3)
                cnt += 1
                if cnt > max2*5:
                    print "limit seeds = ",seed0,seed1
                    break
            out = np.zeros([3,len(out1)])
            out[0] = out1 * (extent[1] - extent[0])/(max1-1) + extent[0]
            out[1] = out2 * (extent[3] - extent[2])/(max2-1) + extent[2]
            out[2] = out3
            field.append(out)
#    import pdb; pdb.set_trace()
    np.save(flines_file, field)        

    return np.array(field)

#============================================================================
# Fieldline Generation
#============================================================================
