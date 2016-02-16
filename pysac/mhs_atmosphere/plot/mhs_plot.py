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
#match the field name to the appropriate axis/colorbar labels. Ordered to 
#control option on, for example, 'mag_pressure'.
#add to the dictionary if variable not specified
plt.ioff()

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
                ('balancing',r'$B\cdot\nabla B/\mu_0$ [Pa m$^{-1}$], '),
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
    if 'density' in var_label:
        var_ = ds.index.grids[0][var_label].in_units('kg / km**3')
    if 'sound' in var_label or 'alfven' in var_label:
        var_ = ds.index.grids[0][var_label].in_units('km / s')

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
    #oneD_arrays[var_label]['Z'   ] = u.Quantity(
    #                                   np.linspace(ds.domain_left_edge[2].value,
    #                                   ds.domain_right_edge[2].value,
    #                                   ds.domain_dimensions[2]),
    #                                   unit=ds.domain_left_edge[2].unit)
    oneD_arrays[var_label]['Z'] = np.linspace(
                                ds.domain_left_edge[2].in_units('Mm'),
                                ds.domain_right_edge[2].in_units('Mm'),
                                ds.domain_dimensions[2])
    return oneD_arrays

def make_1d_zplot(f, plot_label,
                 keys = ['pressure_HS','density_HS','temperature'],
                 subkeys = ['axis'],
                 figxy=[6.47,4.0],
                 ylog = True, xlog = False, loc_legend='center right',
                 empirical = False, ylim = (0,0), xlim = (0,0)
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
        table = atm.hs_atmosphere.read_VAL3c_MTW()
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
            #rescale Bz so camparable in plot with speeds
            if key == 'mag_field_z':
                rescale = 1e3
            else:
                rescale = 1
            if nolabel: 
                plt.plot(f[key]['Z'],
                         f[key][subkey]*rescale,color=color,
                         linestyle=linestyle, lw=linewidth)
            else:
                plt.plot(f[key]['Z'],
                         f[key][subkey],color=color,
                         linestyle=linestyle, lw=linewidth, label=label)
    #collate labels for y-axis and apply log scale if True
    if 'plasma_beta' in keys:
        plt.axhline(1., color='k', ls=':')
    ylabel = ''
    for tag in tags:
        if tag == 'mag_field_z':
            ylabel += r'$B_z$ [mT], '
        else:
            ylabel += y_axis_labels[tag]
    if ylog:
        plt.gca().set_yscale('log',subsy=[5,10])
    if xlog:
        plt.gca().set_xscale('log',subsy=[5,10])
#    plt.gca().set_yticks([1.0e-3,1.0e-1,1.0e1,1.0e3,1.0e5])
    plt.xlabel('Height [Mm]')
    if empirical:
        #limit x axis to simulation data rather than empirical range
        xmax = f[keys[0]]['Z'].max()
        plt.xlim(0,xmax)
    #tidy label string by stripping final ', ' from the end  
    plt.ylabel(ylabel[:-2])
    # allow the axes range to be set manually if required
    if not ylim == (0,0):
        plt.ylim(ylim)
    if not xlim == (0,0):
        plt.xlim(xlim)
    #raise plot x-axis to fit x-label 
    plt.subplots_adjust(bottom=0.125)
    #consider moving legend for different plots, add a loc to function call?
    plt.legend(loc=loc_legend)
    plt.savefig(plot_label)
    plt.close()

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
        colour = cm.bwr
    if 'temperature' in var_field:
        colour = cm.RdBu_r
    if 'beta' in var_field:
        colour = cm.PiYG_r
    if 'tension' in var_field:
        colour = cm.RdYlGn_r
    if 'balancing' in var_field:
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
        norm = colors.SymLogNorm(max(slc.max(),-slc.min()))
        l_norm = True
        cmax = max(slc.max(),-slc.min())
        cmin = -cmax
    elif slc.min() == 0:
        norm = colors.Normalize()
        l_norm = True
        cmin = slc.min()
        cmax = slc.max()
    else:
        norm = colors.LogNorm()
        l_norm = True
        cmin = slc.min()
        cmax = slc.max()

    #plot a 2D slice
    plt.imshow(slc.T, colour,
               extent=extent,
               norm=norm,
               origin='lower',
               vmin=cmin,
               vmax=cmax
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
#    if var_field == "mag_field_y_bg":
#        import pdb; pdb.set_trace()
    cbar1 = plt.colorbar(format = l_f,
                       #ticks=[var.min(),1,10,100,1000,var.max()],
                        fraction=0.2)
    cbar1.ax.set_ylabel(ylabel)
    cbar1.solids.set_edgecolor("face")

    # obtain data arrays for beta contours and magnetic field lines and save
    # for use in plots for other variables - magnetic field must data must be
    # included in gdf file of first run.
    if lines or contours:
        warnings.warn("magnetic field data required for first run to \n"
                     +"obtain field lines or Beta contours", Warning)
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
            X = np.linspace(extent[0], extent[1], slc.shape[0])
            Y = np.linspace(extent[2], extent[3], slc.shape[1])
            #X = u.Quantity(np.linspace(extent[0].value,extent[1].value,
            #               slc.shape[0]), unit=extent[0].unit)   
            #Y = u.Quantity(np.linspace(extent[2].value,extent[3].value,
            #               slc.shape[1]), unit=extent[2].unit)
            U = np.sqrt(v1**2 + v2**2)
            plt.streamplot(X, Y, v1.T, v2.T, color=U.T,
                           cmap=cm.winter, density=line_density, 
                           linewidth = 4*U.T/U.max(), minlength = 0.99
                          ) 
            #add second colorbar for field strength 
            cbar2 = plt.colorbar(format = l_f,
                               #ticks=[var.min(),1,10,100,1000,var.max()],
                                fraction=0.2)
            cbar2.ax.set_ylabel(clabel)
            cbar2.solids.set_edgecolor("face")
        # plot contours with inline labels 
        if contours:
            X = np.linspace(extent[0], extent[1], slc.shape[0])
            Y = np.linspace(extent[2], extent[3], slc.shape[1])
            #X = u.Quantity(np.linspace(extent[0].value,extent[1].value,
            #               slc.shape[0]), unit=extent[0].unit)   
            #Y = u.Quantity(np.linspace(extent[2].value,extent[3].value,
            #               slc.shape[1]), unit=extent[2].unit)
            CS=plt.contour(X,Y,beta.T,4,
                           levels=[
                                   1e-2,1e-1,1.,1e1,1e2,
                                   ],
                           linewidths=1,cmap=cm.PiYG_r,norm=LogNorm())
            if beta.min() < 1e-2:
                fomt = '%.2f'
            elif beta.min() < 1e-1:
                fomt = '%.1f'
            else:
                fomt = '%.0f'
            plt.clabel(CS, fmt=fomt)
    plt.savefig(figname)
    plt.close()

##============================================================================
## Fieldline Generation
##============================================================================
def make_3d_plot(ds, figname, 
                 fields= ['thermal_pressure','plasma_beta',
                           'mag_field_x','mag_field_y','mag_field_z'], 
                 figxy=[500,550],
                 view=(-45., 90., 20., np.array([0,0,3.75])),
                 seeds=np.array([[0,0,1]])
                ):
    """Make a 3D rendition of the atmosphere including volume filling, iso
    surfaces and field lines. 
    ds: gdf data set
    figname: string with path of file to save image
    fields: list of strings indicating the fields to be plotted
            first volume filling, second iso surfaces, third:fifth vector field
    """
    from mayavi import mlab
    from pysac.plot.mayavi_seed_streamlines import SeedStreamline

#    # extract field labels for plot variables
#    vector_field = 'magnetic field'
#    if 'velocity' in fields[-1]:
#        vector_field = 'velocity field'
#    if 'current' in fields[-1]:
#        vector_field = 'electric current'
#    if 'vort' in fields[-1]:
#        vector_field = 'vorticity'
            
#    mlab.options.offscreen = True
    
    scene = mlab.figure(1, bgcolor=(1, 1, 1),
                    fgcolor=(0.5, 0.5, 0.5),size=figxy)
    x,y,z = np.mgrid[ ds.domain_left_edge[0].in_units('Mm'):
                     ds.domain_right_edge[0].in_units('Mm'):
                     1j*ds.domain_dimensions[0],
                      ds.domain_left_edge[1].in_units('Mm'):
                     ds.domain_right_edge[1].in_units('Mm'):
                     1j*ds.domain_dimensions[1],
                      ds.domain_left_edge[2].in_units('Mm'):
                     ds.domain_right_edge[2].in_units('Mm'):
                     1j*ds.domain_dimensions[2]
                    ]

    fill = np.array(ds.index.grids[0][fields[0]])
    surf = np.array(ds.index.grids[0][fields[1]])
    vec1 = np.array(ds.index.grids[0][fields[2]])
    vec2 = np.array(ds.index.grids[0][fields[3]])
    vec3 = np.array(ds.index.grids[0][fields[4]])

    volfill = mlab.pipeline.scalar_field(x,y,z, fill,
                                     name=fields[0])

    isosurf = mlab.pipeline.scalar_field(x,y,z, surf,
                                     name=fields[1])
#    # Set up projection outline to 2D analogue plot
#    o1 = mlab.outline()
#    o1.manual_bounds = True
#    o1.bounds = [x.min(),x.max(),y.min(),y.max(),z.min(),z.max()]
#    o1.actor.property.line_width = 3
#    o1.actor.property.color = (0.5,0.5,0.5)
#    o1.actor.property.line_stipple_pattern = 0x0FF00

    # Generate isosurfaces
    iso = mlab.pipeline.iso_surface(isosurf)
    iso.contour.auto_contours = False
    iso.module_manager.scalar_lut_manager.lut_mode = 'PiYG'  
    iso.module_manager.scalar_lut_manager.lut.scale = 'log10'
    iso.module_manager.scalar_lut_manager.reverse_lut = True 
    iso.contour.contours = [1e-2,1e-1,1,100,10000]   
    iso.module_manager.scalar_lut_manager.data_range = [1e-3,1e3]          
    iso.actor.property.opacity = 0.4

    # Plot the flow of trajectories with suitable parameters.
#    vfield = mlab.flow(x, y, z, vec1, vec2, vec3, 
#                       line_width=1, colormap='winter', seed_scale=2.,
#                       #seedtype='plane',
#                       seed_visible=False, opacity=0.9)
    vecfield = mlab.pipeline.vector_field(x, y, z, vec1, vec2, vec3, 
                                        name=fields[-1][:-3])
    vecmag = mlab.pipeline.extract_vector_norm(vecfield, name="Field line Normals")
    field_lines = SeedStreamline(seed_points = np.array(seeds))
    vecmag.add_child(field_lines)
#    vfield.seed.widget.enabled = False
    field_lines.module_manager.scalar_lut_manager.lut_mode = 'winter'
    field_lines.module_manager.scalar_lut_manager.reverse_lut = False
    field_lines.module_manager.scalar_lut_manager.lut.scale = 'log10'
    field_lines.stream_tracer.integration_direction = 'both'
#    field_lines.stream_tracer.maximum_propagation = 150
#    field_lines.streamline_type = 'tube'
##    field_lines.ribbon_filter.vary_width = True
##    field_lines.ribbon_filter.width_factor = 0.01
#    field_lines.tube_filter.radius_factor = 0.1
#    field_lines.tube_filter.vary_radius = True
#    field_lines.tube_filter.number_of_sides = 3
##    field_lines.streamline_type = 'tube'
    
    # generate rear and lower background surfaces    
    cut = mlab.pipeline.scalar_cut_plane(volfill)
    cut.implicit_plane.normal = [1,0,0]
    cut.implicit_plane.origin = [x.min()*0.999,0,0]
    cut.actor.property.lighting = False
    cut.implicit_plane.widget.enabled = False
    cut.actor.property.opacity = 1.0
    
    cut2 = mlab.pipeline.scalar_cut_plane(volfill)
    cut2.implicit_plane.origin = [0,y.max()*0.999,0]
    cut2.implicit_plane.normal = [0,1,0]
    cut2.actor.property.lighting = False
    cut2.implicit_plane.widget.enabled = False
    cut2.actor.property.opacity = 1.0

    cut3 = mlab.pipeline.scalar_cut_plane(volfill)
    cut3.implicit_plane.origin = [0,0,z.min()+.0001]
    cut3.implicit_plane.normal = [0,0,1]
    cut3.parent.scalar_lut_manager.lut_mode = 'YlOrBr'
    cut3.actor.property.lighting = False
    cut3.parent.scalar_lut_manager.lut.scale = 'log10'
    cut3.implicit_plane.widget.enabled = False
    cut3.parent.scalar_lut_manager.reverse_lut=False
    cut3.actor.property.opacity = 1.0
    
    #Draw the axes and axis labels
    o = mlab.outline()
    o.manual_bounds = True
    o.bounds = [x.min(),x.max(),y.min(),y.max(),z.min(),z.max()]
    o.actor.property.line_width = 3
    o.actor.property.color = (0.5,0.5,0.5)
    o.actor.property.line_stipple_pattern = 0x0FF00
    mlab.axes(x_axis_visibility=True,
              y_axis_visibility=True,
              z_axis_visibility=False,
              zlabel=" Height\n [Mm]  ",
              xlabel="Width\n [Mm]",
              ylabel="Width\n [Mm]",
              nb_labels=3,
              extent=[x.min(),x.max(),y.min(),y.max(),z.min(),z.max()])
    #Define the camera projection for the 3D image
    mlab.view(view[0],view[1],view[2],view[3])
#    mlab.view(-45.0, 90.0, 20.0,
#              np.array([0,  0,   3.7500000e+00]))
 
    #Add colorbars to left (volume fill) and right (isosurfaces)
    unit_str0 = str(ds.index.grids[0][fields[0]].units)
    if str(ds.index.grids[0][fields[0]].units) == 'dimensionless':
        unit_str0 = ''
    label0 = str.replace(fields[0]+'_'+unit_str0, '_','\n')
 #   label0 = str.replace(fields[0], '_','\n')
    cbar0=mlab.scalarbar(object=cut, title=label0, orientation='vertical', 
              nb_labels=None, nb_colors=None, label_fmt=None)        
    cbar0.scalar_bar_representation.position  = [0.80, 0.15]
    cbar0.scalar_bar_representation.position2 = [0.10, 0.80]
    unit_str1 = str(ds.index.grids[0][fields[1]].units)
    if str(ds.index.grids[0][fields[1]].units) == 'dimensionless':
        unit_str1 = ''
    label1 = str.replace(fields[1]+'_'+unit_str1, '_','\n')
 #   label1 = str.replace(fields[1], '_','\n')
    cbar1 = mlab.scalarbar(object=iso, title=label1, orientation='vertical', 
              nb_labels=None, nb_colors=None, label_fmt=None) 
    cbar1.scalar_bar_representation.position  = [0.05, 0.15]
    cbar1.scalar_bar_representation.position2 = [0.10, 0.80]

    mlab.savefig(figname)
    mlab.close() 

