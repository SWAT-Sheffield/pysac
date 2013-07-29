# -*- coding: utf-8 -*-
"""
:Created on: Tue Oct  2 15:10:46 2012

:author: Stuart Mumford
"""
import numpy as np
from mayavi import mlab

def set_text(text_property):
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

def add_axes(ranges):
    axes1 = mlab.axes()
    axes1.axes.ranges = ranges
    axes1.axes.use_ranges = True
    axes1.axes.axis_label_text_property.font_family = 'times'
    axes1.axes.axis_label_text_property.bold = False
    axes1.axes.axis_label_text_property.italic = False
    axes1.axes.axis_label_text_property.color = (0,0,0)
    axes1.axes.label_format = '%3.1f'
    axes1.axes.number_of_labels = 5
    axes1.axes.property.color = (0,0,0)
    axes1.axes.property.line_width = 2.0
    axes1.axes.axis_title_text_property.font_family = 'times'
    axes1.axes.axis_title_text_property.color = (0,0,0)
    axes1.axes.axis_title_text_property.bold = False
    axes1.axes.axis_title_text_property.italic = False
    axes1.axes.x_label = "X [Mm]"
    axes1.axes.y_label = "Y [Mm]"
    axes1.axes.z_label = "Z [Mm]"
    o = mlab.outline(color=(0,0,0))
    o.actor.actor.property.line_width = 1.5
    return axes1, o

def add_cbar_label(cbar,title):
    position = cbar.scalar_bar_representation.position
    position2 = cbar.scalar_bar_representation.position2
    
    x = position[0] + position2[0] + 0.004
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
        raise Exception()
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

def draw_surface(surf_poly,cmap,lines=False,**colourbar_args):
    if not lines:
        surf_poly.lines = None
    new_tube = mlab.pipeline.surface(surf_poly)
    new_tube.module_manager.scalar_lut_manager.lut.table = cmap
    new_tube.parent.parent.point_scalars_name = 'vperp'
    
    surf_bar = add_colourbar(new_tube, [0.84, 0.35], [0.11,0.31],
                                  title='',**colourbar_args)
    surf_bar_label = add_cbar_label(surf_bar,'Velocity Perpendicular\n    to Surface [km/s]')
    #finite = np.isfinite(surf_poly.point_data.scalars)
    lim = np.max([np.nanmax(surf_poly.point_data.scalars),
                  np.abs(np.nanmin(surf_poly.point_data.scalars))])
    new_tube.module_manager.scalar_lut_manager.data_range = np.array([-lim,lim])
    return new_tube, surf_bar, surf_bar_label

def change_surface_scalars(new_tube,surf_bar_label,scalar_name,sym_lim=True,log10=False,lims=None):
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
    lim = np.max([np.nanmax(surf_poly.point_data.scalars),
                      np.abs(np.nanmin(surf_poly.point_data.scalars))])
    if sym_lim:
        new_tube.module_manager.scalar_lut_manager.data_range = np.array([-lim,lim])
    if lims:
        new_tube.module_manager.scalar_lut_manager.data_range = np.array(lims)
    if log10:
        new_tube.module_manager.scalar_lut_manager.lut.scale = 'log10'
    