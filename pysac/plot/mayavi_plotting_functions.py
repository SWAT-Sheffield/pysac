# -*- coding: utf-8 -*-
"""
This module provides a set of helper functions to make mayavi scenes into 
nice plots. It also provides routines for the visualisation of flux surfaces.
"""
import numpy as np

from mayavi import mlab

__all__ = ['set_text', 'set_cbar_text', 'add_axes', 'add_cbar_label',
           'add_colourbar', 'draw_surface', 'change_surface_scalars']

def set_text(text_property):
    """
    Set the text to sane defaults
    
    Parameters
    ----------
    text_property: mayavi TextProperty object
    """
    text_property.bold = False
    text_property.italic = False
    text_property.font_family = 'times'
    text_property.color = (0,0,0)
    
def set_cbar_text(cbar,lut_manager='scalar'):
    if lut_manager == 'scalar':
        manager = cbar.parent.scalar_lut_manager
    elif lut_manager == 'vector':
        manager = cbar.parent.vector_lut_manager
    else:
        raise Exception()
    set_text(manager.title_text_property)
    set_text(manager.label_text_property)

def add_axes(ranges, obj=None):
    axes1 = mlab.axes(obj)
    axes1.axes.ranges = ranges
    axes1.axes.use_ranges = True
    set_text(axes1.axes.axis_label_text_property)
    axes1.axes.label_format = '%3.1f'
    axes1.axes.number_of_labels = 5
    axes1.axes.property.color = (0,0,0)
    axes1.axes.property.line_width = 2.0
    set_text(axes1.axes.axis_title_text_property)
    axes1.axes.x_label = "X [Mm]"
    axes1.axes.y_label = "Y [Mm]"
    axes1.axes.z_label = "Z [Mm]"
    o = mlab.outline(color=(0,0,0))
    o.actor.actor.property.line_width = 1.5
    return axes1, o

def add_cbar_label(cbar,title):
    position = cbar.scalar_bar_representation.position
    position2 = cbar.scalar_bar_representation.position2
    
    x = position[0] + position2[0] #+ 0.002
    y = position[1]
    
    text = mlab.text(x,y,title)
    text.property.font_size = 24
    set_text(text.property)
    
    text.property.orientation = 90.
    text.width = 0.205
    
    return text


def add_colourbar(module, position ,position2, title,label_fstring='%#4.2f',
                  number_labels=5, orientation=1,lut_manager='scalar'):
    set_cbar_text(module,lut_manager=lut_manager)
    if lut_manager == 'scalar':
        bar = module.parent.scalar_lut_manager.scalar_bar_widget
        bar_parent = module.parent.scalar_lut_manager
    elif lut_manager == 'vector':
        bar = module.parent.vector_lut_manager.scalar_bar_widget
        bar_parent = module.parent.vector_lut_manager
    else:
        raise ValueError("Please provide 'scalar' or 'vector' for LUT")
    bar.enabled = True
    bar.scalar_bar_representation.orientation = orientation
    bar.repositionable = False
    bar.resizable = False
    bar.scalar_bar_representation.position = position
    bar.scalar_bar_representation.position2 = position2
    bar_parent.number_of_labels = number_labels
    bar_parent.data_name = title
    bar.scalar_bar_actor.label_format = label_fstring
    return bar

def _parse_limits(lim):
    #Parse lim
    if isinstance(lim, float) or isinstance(lim, int):
        lim = [-lim,lim]
    elif np.shape(lim) == (2,):
        pass #lim = lim
    else:
        raise ValueError("lim will accept [None, float, int or [min,max]")
    return lim

def draw_surface(surf_poly, cmap, scalar='vperp', lines=False, lim=None, log10=False,
                 colorbar_label='Velocity Perpendicular\n    to Surface [km/s]',
                 **colourbar_args):
    """
    Draw a mayavi surface from a PolyData object.
    
    Parameters
    ----------
    surf_poly: tvtk.PolyData Object
        The polydata object containing the surface information.
    
    cmap: ndarray or string
        The colour map to scale the scalar data with.
    
    scalar: string
        The name of the tvtk scalar to plot.
    
    lines: (optional) bool
        If True then lines are not removed from the PolyData object and will be plotted.
    
    lim: [None, float, int or [min,max]]
        If lim is None symmetric limits will be created covering the whole dynamic range,
        if the input is a number symmetric limts will e made with that as a maximum,
        if a [min,max] pair is passed that will be the limits.
    
    log10: bool
        Scale the scalar with log scaling.
    
    colorbar_label: string
        Label the colorbar.
    
    **colorar_args: dict
        Other kwargs will be passed to mayavi_plotting_functions.add_colourbar()
    
    Returns
    -------
    new_tube: mayavi.modules.surface.Surface
        The surface object
    
    surf_bar: ScalarBarWidget
        The colorbar object
    
    surf_bar_label: mayavi.modules.text.Text
        The object created to label the colorbar
    """
    #Remove lines
    if not lines:
        surf_poly.lines = None
    
    #Create the surface
    new_tube = mlab.pipeline.surface(surf_poly)
    
    #Set the surface scalar
    new_tube.parent.parent.point_scalars_name = scalar
    
    #Parse the colour map
    if isinstance(cmap, np.ndarray):
        new_tube.module_manager.scalar_lut_manager.lut.table = cmap
    elif isinstance(cmap, basestring):
        new_tube.module_manager.scalar_lut_manager.lut_mode = cmap
    else:
        raise TypeError("cmap should be an array or a string")

    #Set the arguments to the colour bar function
    cbar_args = {'position':[0.84, 0.35],
                 'position2':[0.11,0.31],
                 'title':''}
    cbar_args.update(colourbar_args)
    surf_bar = add_colourbar(new_tube, **cbar_args)
    surf_bar_label = add_cbar_label(surf_bar,colorbar_label)

    if log10:
        new_tube.module_manager.scalar_lut_manager.lut.scale = 'log10'
    if not log10:
        new_tube.module_manager.scalar_lut_manager.lut.scale = 'linear'

    #Parse lim
    if lim is None:
        lim = np.max([np.nanmax(surf_poly.point_data.scalars),
                      np.abs(np.nanmin(surf_poly.point_data.scalars))])
        lim = [-lim,lim]
    lim = _parse_limits(lim)
    new_tube.module_manager.scalar_lut_manager.data_range = np.array(lim)

    
    return new_tube, surf_bar, surf_bar_label

def change_surface_scalars(new_tube, surf_bar_label, scalar_name, colorbar_label=None,
                           lim=None, log10=False):
    """
    Change the surface scalars of an existing surface object, and change the associated
    colorbar and text objects. Note: This is hard coded to be easier for velocity surfaces.
    
    Parameters
    ----------
    new_tube: mayavi.modules.surface.Surface
        The Surface
    
    surf_bar_label: mayavi.modules.text.Text
        The text object

    scalar_name: string
        The name to set the scalar to.
    
    colorbar_label: string
        The string to set the colorbar tex object to, if None it will guess velocity ones.
    
    lim:[None, float, int or [min,max]]
        If lim is None symmetric limits will be created covering the whole dynamic range,
        if the input is a number symmetric limts will e made with that as a maximum,
        if a [min,max] pair is passed that will be the limits.
    
    log10: bool
        Scale the scalar with log scaling.
    
    Returns
    -------
    None
    """
    new_tube.parent.parent.point_scalars_name = scalar_name
    if scalar_name == 'vperp':
        title = 'Velocity Perpendicular\n    to Surface [km/s]'
    elif scalar_name == 'vpar':
        title = '  Velocity Parallel \nto Surface [km/s]'
    elif scalar_name == 'vphi':
        title = 'Velocity Azimuthally\n   to Surface [km/s]'
    else:
        title='unknown'
    new_tube.parent.scalar_lut_manager.data_name = ''
    surf_bar_label.text = title
    surf_poly = new_tube.parent.parent.outputs[0]
    
    #Parse lim
    if lim is None:
        lim = np.max([np.nanmax(surf_poly.point_data.scalars),
                      np.abs(np.nanmin(surf_poly.point_data.scalars))])
        lim = [-lim,lim]
    lim = _parse_limits(lim)
    new_tube.module_manager.scalar_lut_manager.data_range = np.array(lim)
    
    if log10:
        new_tube.module_manager.scalar_lut_manager.lut.scale = 'log10'
    if not log10:
        new_tube.module_manager.scalar_lut_manager.lut.scale = 'linear'
    