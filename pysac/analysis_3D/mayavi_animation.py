# -*- coding: utf-8 -*-
"""
:Created on: Tue Oct  2 15:12:43 2012

:author: Stuart Mumford
"""
import os
import sys
import numpy as np
from mayavi import mlab

import pysac.io as sacio

import mayavi_plotting_functions as mpf
import tvtk_tube_functions as ttf
from pysac.plot.CustomColourmaps import get_mayavi_colourmap
from mayavi_cust_streamlines import sStreamline

driver = "Sarch"
post_amp = "A10"
period = "p240"
raw_amp = "A10"
tube_r = "r60"
exp_fac = "B0005"

if driver in ['Slog', 'Sarch', 'Suni']:
    out_dir = "/home/stuart/_PhD/VAC/Visualisations/output_project1/%s/%s_%s_%s_%s/frames/"%(driver,period,post_amp,tube_r,exp_fac)
else:
    out_dir = "/home/stuart/_PhD/VAC/Visualisations/output_project1/%s/%s_%s_%s/frames/"%(driver,period,post_amp,tube_r)

data_prefix = "/data/SAC Data/Processed"

if driver in ['Slog', 'Sarch', 'Suni']:
    data_dir = os.path.join(data_prefix,"%s/%s_%s_%s_%s/"%(driver,period,post_amp,tube_r,exp_fac))
else:
    data_dir = os.path.join(data_prefix,"%s/%s_%s_%s/"%(driver,period,post_amp,tube_r))
        
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

h5_dir = "/home/stuart/iceberg_fastdata/SAC Output"
#h5_dir = "/home/stuart/barcelona/"
#h5_dir = "/data/SAC Data/"
#print os.path.join(h5_dir,"3D_tube128_%s_%s_C_%s_600s_np020204.h5"%(driver,period,raw_amp))
f = sacio.SACdata(os.path.join(h5_dir,"3D_tube128_%s_%s_C_%s_%s_600s_np020204.h5"%(driver,period,raw_amp,exp_fac)))

def path_join(filename):
    return os.path.join(data_dir,filename)
    
def out_path_join(filename):
    return os.path.join(out_dir,filename)
    
cube_slice = np.s_[:,:,:-5]
x_slice = np.s_[:,:,:,:-5]

xmax,ymax,zmax = np.array(f.w_sac['rho'][cube_slice].shape)-1
#Get Secondary Varibles
va_f = f.get_va()
cs_f = f.get_cs()
thermal_p,mag_p = f.get_thermalp(beta=True)
beta_f = mag_p / thermal_p
#Convert B to Tesla
f.convert_B()

#==============================================================================
# Create Mayavi Scene
#==============================================================================
scene = mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),size=[1200,900])

#==============================================================================
# Create Mayavi datasets from HDF5 file
#==============================================================================
bfield = mlab.pipeline.vector_field(f.w_sac['b3'][cube_slice] * 1e3,
                                    f.w_sac['b2'][cube_slice] * 1e3, 
                                    f.w_sac['b1'][cube_slice] * 1e3,
                                    name="Magnetic Field")
bmag = mlab.pipeline.extract_vector_norm(bfield, name="Field line Normals")

vfield = mlab.pipeline.vector_field(f.w_sac['v3'][cube_slice] / 1e3,
                                    f.w_sac['v2'][cube_slice] / 1e3,
                                    f.w_sac['v1'][cube_slice] / 1e3,
                                    name="Velocity Field")
v_mag = mlab.pipeline.extract_vector_norm(vfield, name="Velocity Magnitude")
base_vel = mlab.pipeline.extract_vector_components(vfield,
                                                   name="Base Velocity")
density = mlab.pipeline.scalar_field(f.w_sac['rho'][cube_slice],
                                     name="Density")
valf = mlab.pipeline.scalar_field(va_f,name="Alven Speed")
cs = mlab.pipeline.scalar_field(cs_f)
beta = mlab.pipeline.scalar_field(beta_f)

#==============================================================================
# Add a timestamp
#==============================================================================
text = mlab.text(0.75,0.01,"t=%3.2f s"%f.header['t'])
mpf.set_text(text.property)
text.actor.text_scale_mode = 'none'
text.property.font_size = 24

#==============================================================================
# Add Axes (Has to be done after first object added)
#==============================================================================
#axes, outline = add_axes(np.array([0.0,2.0,0.0,2.0,0.0,1.6]))
x = np.array(f.x)
axes, outline = mpf.add_axes(np.array([np.min(x[x_slice][2,:,0,0])/1e6, 
                                   np.max(x[x_slice][2,:,0,0])/1e6, 
                                   np.min(x[x_slice][1,0,:,0])/1e6,
                                   np.max(x[x_slice][1,0,:,0])/1e6,
                                   np.min(x[x_slice][0,0,0,:])/1e6,
                                   np.max(x[x_slice][0,0,0,:])/1e6]))

#==============================================================================
# Add velocity cut plane at base, with colourmap
#==============================================================================
vvel = mlab.pipeline.vector_cut_plane(vfield)
vvel.implicit_plane.normal = [0,0,1]
vvel.implicit_plane.origin = [64.5,64.5,5.0]
vvel.implicit_plane.visible = False
vvel.glyph.glyph.scale_factor = 50
vvel.glyph.mask_points.maximum_number_of_points = 100000000 #vtk >6.0 bug fix
vvel.glyph.mask_input_points = True
vvel.glyph.mask_points.on_ratio = 10

div_vvel = get_mayavi_colourmap([1.0,0.0,0.0], [0.0,1.0,0.0], 0.0, 256, 
                                5, 0.0, 255)
                                
vvel.module_manager.vector_lut_manager.lut.table = div_vvel
vmag = np.sqrt(f.w_sac['v1'][cube_slice][:,:,5]**2 + 
               f.w_sac['v2'][cube_slice][:,:,5]**2 +
               f.w_sac['v3'][cube_slice][:,:,5]**2)
               
lim = np.max([np.max(vmag),np.abs(np.min(vmag))])

vvel.module_manager.vector_lut_manager.data_range = np.array([0,
                                                              lim/1e3])
vv_bar = mpf.add_colourbar(vvel,[0.84,0.03], [0.11,0.31],
                       '',lut_manager='vector')
vv_text = mpf.add_cbar_label(vv_bar,'Driver Plane Velocity\n          [km/s]')
                       
#==============================================================================
# Generate seeds for fieldlines and plot them in the image
#==============================================================================
xc = xmax/2
yc = ymax/2
seeds = [[xc,yc,zmax]]
for ti,r in enumerate(np.linspace(15, yc, 3)):
    for theta in np.linspace(0, 2 * np.pi, 4, endpoint=False):
        seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                      r * np.sin(theta + 0.5 * ti) + yc, zmax])
seeds = np.array(seeds)
seed_points = mlab.points3d(seeds[:,0], seeds[:,1], seeds[:,2],
                            color=(0.231, 0.298, 0.752), scale_mode='none',
                            scale_factor=1.8)

#==============================================================================
# Calculate and plot streamlines and add top panel field strength contours
#==============================================================================
#div_flines = get_mayavi_colourmap([0.231,0.298,0.752], [1.0,0.0126,0.0], 0.05,
                                  #256, 95, 0, 255)
#div_flines = get_mayavi_colourmap([0.231,0.298,0.752], [1.0,0.0126,0.0], 0.5, 256, 5, 0.5, 255)
field_lines = sStreamline(seed_points = np.array(seeds))
bmag.add_child(field_lines)
#field_lines.module_manager.scalar_lut_manager.lut.table = div_flines
field_lines.module_manager.scalar_lut_manager.lut_mode = 'winter'
field_lines.module_manager.scalar_lut_manager.reverse_lut = True
field_lines.module_manager.scalar_lut_manager.data_range = np.array([0.01,100])
field_lines.parent.scalar_lut_manager.lut.scale = 'log10'

flines_cbar = mpf.add_colourbar(field_lines, [0.84,0.68], [0.11,0.31],
              '')
vv_text = mpf.add_cbar_label(flines_cbar,'Magnetic Field Strength\n              [mT]')
#Contour the top panel
cntr = mlab.pipeline.contour_grid_plane(bmag, name="Top B Contours")
cntr.grid_plane.axis = 'z'
cntr.grid_plane.position = zmax-1
cntr.contour.minimum_contour = 0.01
cntr.contour.maximum_contour = 0.5
cntr.contour.number_of_contours = 5
#cntr.module_manager.scalar_lut_manager.lut.table = div_flines
#cntr.module_manager.scalar_lut_manager.reverse_lut = False
cntr.module_manager.scalar_lut_manager.lut_mode = 'winter'
surf_poly = ttf.read_step(path_join('Fieldline_surface_%s_%s_%s_00000.vtp'%(driver,period,post_amp)))

div_perp = get_mayavi_colourmap([0.,1.,1.], [1.,0.,1.], 0.5 ,256,
                             100, 0, 190)
new_tube, surf_bar, surf_bar_label = mpf.draw_surface(surf_poly,div_perp)

def set_view():
    mlab.view(153.1813884481075, 62.474202967954426, 417.66115037685745, 
              np.array([ 75.21308804,  16.09446165,  29.63298413]))
set_view()
#mlab.show()


def animate_scene(n):
    global new_tube, surf_bar_label
    """Animate scence to the specified frame number"""
#==============================================================================
#   Read in New Data
#==============================================================================
    t1 = f.header['t'] #dt needed for seed points
    f.readrecord(n)
    f.convert_B()
    t2 = f.header['t']
    text.text = "t=%3.2f s"%t2
#==============================================================================
#   Modify fieldline seeds with old velocity data
#==============================================================================
    new_seeds = ttf.move_seeds(field_lines.stream_tracer.source,vfield ,t2-t1)
#==============================================================================
#   Update Datasets
#==============================================================================
    bfield.set(vector_data = np.rollaxis(np.array([f.w_sac['b3'][cube_slice] * 1e3,
                                                   f.w_sac['b2'][cube_slice] * 1e3,
                                                   f.w_sac['b1'][cube_slice] * 1e3]),
                                                   0, 4))
    vfield.set(vector_data = np.rollaxis(np.array([f.w_sac['v3'][cube_slice] / 1e3,
                                                   f.w_sac['v2'][cube_slice] / 1e3,
                                                   f.w_sac['v1'][cube_slice] / 1e3]),
                                                   0, 4))
    
#==============================================================================
#   Update Colour bar x
#==============================================================================
    vmag = np.sqrt(f.w_sac['v1'][cube_slice][:,:,5]**2 +
                    f.w_sac['v2'][cube_slice][:,:,5]**2 + 
                    f.w_sac['v3'][cube_slice][:,:,5]**2)
    lim = np.max([np.max(vmag),np.abs(np.min(vmag))])
    vvel.module_manager.vector_lut_manager.data_range = np.array([0,
                                                                  lim/1e3])

#==============================================================================
#   Update Seeds and contours
#==============================================================================
    field_lines.seed_points = new_seeds
    field_lines.stream_tracer.update()
    seed_points.mlab_source.points = new_seeds
    cntr.contour.maximum_contour = 0.5 #This forces redraw of contours

#==============================================================================
#   Update the Surface
#==============================================================================
    surf_poly = ttf.read_step(path_join('Fieldline_surface_%s_%s_%s_%05i.vtp'%(driver,period,post_amp,n)))
    
    new_tube.parent.parent.remove()
    new_tube, surf_bar, surf_bar_label = mpf.draw_surface(surf_poly,div_perp)
    #Weird Colourbar hack
    field_lines.module_manager.scalar_lut_manager.data_range = np.array([0.01,100])
    mpf.change_surface_scalars(new_tube,surf_bar_label,'vphi')
    new_tube.module_manager.scalar_lut_manager.data_range = np.array([-0.1,0.1])    
    set_view()


def save_animation(step):
    import time
    t = time.clock()
    for i in range(0,f.num_records,step):
        t1 = time.clock()
        print i
        try:
            animate_scene(i)
        except Exception as e:
            print e
        mlab.savefig(out_path_join("Mayavi_%s_%s_%s_%s_vperp_%05i.png"%(driver,period,post_amp,tube_r,i)),size=[1200,900],magnification=1)
        mpf.change_surface_scalars(new_tube,surf_bar_label,'vpar')
        mlab.savefig(out_path_join("Mayavi_%s_%s_%s_%s_vpar_%05i.png"%(driver,period,post_amp,tube_r,i)),size=[1200,900],magnification=1)
        mpf.change_surface_scalars(new_tube,surf_bar_label,'vphi')
        mlab.savefig(out_path_join("Mayavi_%s_%s_%s_%s_vphi_%05i.png"%(driver,period,post_amp,tube_r,i)),size=[1200,900],magnification=1)
        print "step %i done in %f s\n"%(i,time.clock()-t1)+"_"*80+"\n"
    t2 = time.clock()
    print "All done in %f s %f min\n"%(t2-t,t2-t/60.)+"_"*80+"\n"
    from pushover import pushover
    pushover(title="Animation Complete",
             message="The animation script for driver:%s, period:%s and Amplitude:%s completed in %f seconds"%(driver,period,post_amp,t2-t))