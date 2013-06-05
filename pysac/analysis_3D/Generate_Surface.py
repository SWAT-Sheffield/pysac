# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 11:07:36 2012

@author: Stuart Mumford

This routine processes a whole SAC output file, in HDF5 format, and calculates
flux and all related varibles.

It saves all surface varibles to a vtk file for each time step and all line
vars to npy files.
"""

import time
import os
import sys

import numpy as np
import pysac.io as sacio

import process_utils as util
import tvtk_tube_functions as ttf

#==============================================================================
# Get paths right and set up problem
#==============================================================================
print sys.argv, len(sys.argv)
if len(sys.argv) == 1:
    driver = "Slog"
    post_amp = "A10"
    period = "p240"
    raw_amp = "A10"
    tube_r = "r30"
    exp_fac = "B05"

    data_dir = "/data/SAC Data/Processed/%s/vtk_flines_%s_%s/"%(driver,
                                                                period,
                                                                post_amp)
    h5_dir = "/home/stuart/iceberg_fastdata/SAC Output/"
    
else:
    if len(sys.argv) > 5:
        data_dir, h5_dir, [driver, post_amp, period,
                           raw_amp, tube_r, exp_fac] = util.read_argv(sys.argv)
    else:
        [driver, post_amp, period, 
         raw_amp, tube_r, exp_fac] = util.read_argv(sys.argv)
    
f = sacio.SACdata(os.path.join(
                    h5_dir,"3D_tube128_%s_%s_C_%s_%s_600s_np020204.h5"%(driver,
                                                                     period,
                                                                     raw_amp,
                                                                     exp_fac)))
def path_join(filename):
    return os.path.join(data_dir,filename)

#==============================================================================
# Set some Parameters then Read in HDF5 data
#==============================================================================
cube_slice = np.s_[:,:,:-5]
x_slice = np.s_[:,:,:,:-5]
xmax, ymax, zmax = np.array(f.w_sac['rho'][cube_slice].shape) - 1
#nlines is the number of fieldlines used in the surface
n_lines = 100
#the line is the fieldline to use as "the line"
line_n = 25

#Read in all data from hdf5 file
bfield, vfield, density, valf, cs, beta = util.get_hdf5(f, 
                                                        cube_slice, flux=True)

#Define domain parameters
domain = {'xmax':xmax,'ymax':ymax,'zmax':zmax}

#Make inital surface seeds (n=100, r=45)
surf_seeds = ttf.make_circle_seeds(n_lines, int(tube_r[1:]), **domain)

#Make surface using seeds, surf filter and contour
surf_field_lines, surface = ttf.create_flux_surface(bfield, surf_seeds)

#Make the PolyDataNormals object
poly_norms = ttf.make_poly_norms(surface)

#Interpolate the vfield to the surface
surface_vel_filter, surface_velocities = ttf.interpolate_vectors(vfield.outputs[0],
                                                              poly_norms.output)
surface_mag_filter, surface_bfield = ttf.interpolate_vectors(bfield.outputs[0],
                                                         poly_norms.output)

#Interpolate the vfield to the surface
surface_den_filter, surface_density = ttf.interpolate_vectors(density.outputs[0],
                                                              poly_norms.output)
surface_va_filter, surface_va = ttf.interpolate_vectors(valf.outputs[0],
                                                              poly_norms.output)
surface_cs_filter, surface_cs = ttf.interpolate_vectors(cs.outputs[0],
                                                              poly_norms.output)
surface_beta_filter, surface_beta = ttf.interpolate_vectors(beta.outputs[0],
                                                              poly_norms.output)
#Make the line
the_line = ttf.get_the_line(bfield, surf_seeds, line_n)

#==============================================================================

#Set up arrays to store all the varibles for all the timesteps
save_times = np.array([i[1].attrs['t'] for i in f.file.time_group.items()])

save_points = np.zeros([f.num_records, len(the_line.output.points),3])
save_vpar = np.zeros([f.num_records, len(the_line.output.points)])
save_vperp = np.zeros([f.num_records, len(the_line.output.points)])
save_vphi = np.zeros([f.num_records, len(the_line.output.points)])
save_index = np.zeros([f.num_records, len(the_line.output.points)])


Fpar_line = np.zeros([f.num_records, len(the_line.output.points)])
Fperp_line = np.zeros([f.num_records, len(the_line.output.points)])
Fphi_line = np.zeros([f.num_records, len(the_line.output.points)])

density_line = np.zeros([f.num_records, len(the_line.output.points)])
va_line = np.zeros([f.num_records, len(the_line.output.points)])
beta_line = np.zeros([f.num_records, len(the_line.output.points)])
cs_line = np.zeros([f.num_records, len(the_line.output.points)])

#==============================================================================

def save_step(n):
    global the_line, bfield, vfield, density, valf, cs, beta
    t1 = f.header['t'] #dt needed for seed points
    f.read_timestep(n)
    t2 = f.header['t']
    
    [bfield, vfield, density,
     valf, cs, beta] = util.process_next_step(f, cube_slice, bfield, vfield,
                                              density, valf, cs, beta)

    #Move fieldlines with current step velocities
    surf_seeds = ttf.move_seeds(surf_field_lines.source, vfield, t2-t1)
    surf_field_lines.source.points = surf_seeds
    
    #Update surface
    ttf.update_flux_surface(surf_seeds, surf_field_lines, surface)
    #Get velocities at surface
    surface_velocities = ttf.update_interpolated_vectors(False,
                                                         surface_vel_filter)
    #Get surface magfield
    surface_bfield = ttf.update_interpolated_vectors(False, surface_mag_filter)
    #Get vectors at surface
    normals, torsionals, parallels = ttf.get_surface_vectors(poly_norms,
                                                         surface_bfield)
    #Get velocity componets at surface
    vperp, vpar, vphi = ttf.get_surface_velocity_comp(surface_velocities,
                                                      normals, torsionals,
                                                      parallels)
    #Save to file
    ttf.write_step(path_join('Fieldline_surface_%s_%s_%s_%s_%s_%05i.vtp'%(driver, period, post_amp, tube_r, exp_fac, n)),
               surface,normals,parallels,torsionals,vperp,vpar,vphi)

    the_line = ttf.update_the_line(the_line, surf_seeds, line_n, len(the_line.output.points))

    #line varibles
    surf_line_index, surf_line_points = ttf.get_surface_indexes(surface.output,the_line)
    
    
    save_index[n] = surf_line_index
    save_points[n] = surf_line_points
    save_vpar[n] = vpar[surf_line_index]
    save_vperp[n] = vperp[surf_line_index]
    save_vphi[n] = vphi[surf_line_index]

    step_poly = surface.output

    
    surface_density = ttf.update_interpolated_scalars(step_poly, surface_den_filter)
    surface_va = ttf.update_interpolated_scalars(step_poly, surface_va_filter)
    surface_beta = ttf.update_interpolated_scalars(step_poly, surface_beta_filter)
    surface_cs = ttf.update_interpolated_scalars(step_poly, surface_cs_filter)
    
    #Calculate Fluxs
    vpar = ttf.get_data(step_poly,'vpar') * 1e3 #stored in km/s
    if vpar.shape != surface_density.shape or surface_density.shape != surface_cs.shape:
        if n > 450: #If we fuck up far enough in, save it anyway.
            save_all_data()
        import pdb; pdb.set_trace()
    Fpar = surface_density * vpar**2 * surface_cs
    
    vperp = ttf.get_data(step_poly,'vperp') * 1e3 #stored in km/s
    Fperp = surface_density * vperp**2 * np.sqrt(surface_cs**2 + surface_va**2)
    
    vphi = ttf.get_data(step_poly,'vphi') * 1e3 #stored in km/s
    Fphi = surface_density * vphi**2 * surface_va

    ttf.write_flux(path_join(
        "SurfaceFlux_AllVars_%s_%s_%s_%05i.vtp"%(driver, period, post_amp, n)),
                             step_poly, surface_density, surface_va,
                             surface_beta, surface_cs, Fpar, Fperp, Fphi)    
    
    Fpar_line[n] = Fpar[surf_line_index]
    Fperp_line[n] = Fperp[surf_line_index]
    Fphi_line[n] = Fphi[surf_line_index]
    
    density_line[n] = surface_density[surf_line_index]
    va_line[n] = surface_va[surf_line_index]
    beta_line[n] = surface_beta[surf_line_index]
    cs_line[n] = surface_cs[surf_line_index]
    
    print "Done Flux %i"%n

t = time.clock()
for n in xrange(0,f.num_records):
    t2 = time.clock()
    save_step(n)
    print "Step %i done in %f s\n"%(n,time.clock()-t2)+"_"*80
    
#==============================================================================
#Save final line arrays to file
#==============================================================================
def save_all_data():
    np.save(path_join("LineVar_%s_%s_%s_%s_%s_index.npy"%(driver, period,
                                                          post_amp, tube_r, exp_fac)),save_index)
    np.save(path_join("LineVar_%s_%s_%s_%s_%s_vphi.npy"%(driver, period,
                                                         post_amp, tube_r, exp_fac)),save_vphi)
    np.save(path_join("LineVar_%s_%s_%s_%s_%s_vperp.npy"%(driver, period,
                                                          post_amp, tube_r, exp_fac)),save_vperp)
    np.save(path_join("LineVar_%s_%s_%s_%s_%s_vpar.npy"%(driver, period,
                                                         post_amp, tube_r, exp_fac)),save_vpar)
    np.save(path_join("LineVar_%s_%s_%s_%s_%s_points.npy"%(driver, period,
                                                           post_amp, tube_r, exp_fac)),save_points)
    np.save(path_join("LineVar_%s_%s_%s_%s_%s_times.npy"%(driver, period,
                                                          post_amp, tube_r, exp_fac)),save_times)
    
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_Fpar.npy"%(driver, period,
                                                          post_amp, tube_r, exp_fac)),Fpar_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_Fperp.npy"%(driver, period,
                                                           post_amp, tube_r, exp_fac)),Fperp_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_Fphi.npy"%(driver, period,
                                                          post_amp, tube_r, exp_fac)),Fphi_line)
    
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_rho.npy"%(driver, period,
                                                         post_amp, tube_r, exp_fac)),density_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_va.npy"%(driver, period,
                                                        post_amp, tube_r, exp_fac)),va_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_cs.npy"%(driver, period,
                                                        post_amp, tube_r, exp_fac)),cs_line)
    np.save(path_join("LineFlux_%s_%s_%s_%s_%s_beta.npy"%(driver, period,
                                                          post_amp, tube_r, exp_fac)),beta_line)

save_all_data()
print "All Done in %f s"%(time.clock()-t)