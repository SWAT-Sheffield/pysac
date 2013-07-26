#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os
import sys
sys.path.append("/home/stuart/BitBucket/pySAC/pysac/analysis_3D/") 
import copy

import numpy as np
import pysac.io

import process_utils as util
import tvtk_tube_functions as ttf
from tvtk.api import tvtk
from mayavi import mlab
import h5py

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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

    data_dir = "/archive/processed/data/parallel_testing/mpi_test/%s/vtk_flines_%s_%s/"%(driver,
                                                                period,
                                                                post_amp)
    h5_dir = "/archive/simulations/driver_paper/"
    
else:
    if len(sys.argv) > 5:
        data_dir, h5_dir, [driver, post_amp, period,
                           raw_amp, tube_r, exp_fac] = util.read_argv(sys.argv)
    else:
        [driver, post_amp, period, 
         raw_amp, tube_r, exp_fac] = util.read_argv(sys.argv)
    
hdf5_fname = os.path.join(
                    h5_dir,"3D_tube128_%s_%s_C_%s_%s_600s_np020204.h5"%(driver,
                                                                     period,
                                                                     raw_amp,
                                                                     exp_fac))
f  = pysac.io.SACdata(hdf5_fname)
#Define a var to limit iterations, no limt = f.num_records
max_n = f.num_records
                                                                     
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

#==============================================================================
# Do Serial compute of seed points
# This is designed to have absolute minimal I/O over heads to reduce execution
# time. This will hopfully read the least data needed.
#==============================================================================
if rank == 0:
    seeds_slice = np.s_[:,:,-6]#slicing upto -5 is not the same as indexing -5
    h5file = h5py.File(hdf5_fname,'r')
    sac_group = h5file["SACdata"]
    time_group = sac_group['wseries']
    wstepname = time_group.keys()
    
    #Don't reader header define raw indcies
    m1_ = 1; m2_ = 2; m3_ = 3; h_ = 0; rhob_ = 9
    
    rho = time_group[wstepname[0]][h_][seeds_slice] +  time_group[wstepname[0]][rhob_][seeds_slice]
    v3 = time_group[wstepname[0]][m3_][seeds_slice] / rho
    v2 = time_group[wstepname[0]][m2_][seeds_slice] / rho
    v1 = time_group[wstepname[0]][m1_][seeds_slice] / rho
    
    vfield = mlab.pipeline.vector_field(v3 / 1e3,
                                        v2 / 1e3,
                                        v1 / 1e3,
                                            name="Velocity Field",figure=None)
    #Define domain parameters
    domain = {'xmax':xmax,'ymax':ymax,'zmax':0}
    #Create initial seed points in tvtk.PolyData
    surf_seeds_poly = ttf.make_circle_seeds(n_lines, int(tube_r[1:]), **domain)
    #Make surface using seeds, surf filter and contour
    #Extract all times
    save_times = np.array([i[1].attrs['t'] for i in f.file.time_group.items()])
    next_seeds = np.array(surf_seeds_poly.points)
    next_seeds[:,-1] = zmax
    surf_seeds = [next_seeds]
    for i in range(1,max_n):
        print i
        t1 = save_times[i-1]
        t2 = save_times[i]
        rho = time_group[wstepname[i]][h_][seeds_slice] +  time_group[wstepname[i]][rhob_][seeds_slice]
        v3 = time_group[wstepname[i]][m3_][seeds_slice] / rho
        v2 = time_group[wstepname[i]][m2_][seeds_slice] / rho
        v1 = time_group[wstepname[i]][m1_][seeds_slice] / rho
        vfield = mlab.pipeline.vector_field(v3 / 1e3,
                                            v2 / 1e3,
                                            v1 / 1e3,
                                                name="Velocity Field",figure=None)
        next_seeds = ttf.move_seeds(surf_seeds_poly, vfield, t2-t1)
        next_seeds[:,-1] = zmax
        surf_seeds.append(next_seeds)
else:
    surf_seeds = None

## tvtk objects are presumeably not pikleable, therefore cannot be trasmitted
## via MPI, create the tvtk object upon recieve.
surf_seeds = np.array(surf_seeds)
#Give all processes the surf_seeds
surf_seeds_arr = comm.bcast(surf_seeds, root=0)

surf_seeds = []
for seeds in surf_seeds_arr:
    pd = tvtk.PolyData()
    pd.points = seeds
    surf_seeds.append(pd)
#Divide up the time steps to each process
all_indices = np.arange(0,max_n,1,dtype=int)
chunks = np.array_split(all_indices, size)

rank_indices = comm.scatter(chunks, root=0)

print rank_indices

#==============================================================================
# Begin Parallel compute
#==============================================================================

#Read in all data from hdf5 file
bfield, vfield, density, valf, cs, beta = util.get_hdf5(f, 
                                                        cube_slice, flux=True)

#Make surface using seeds, surf filter and contour
surf_field_lines, surface = ttf.create_flux_surface(bfield, surf_seeds[0])

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
the_line = ttf.get_the_line(bfield, surf_seeds[0], line_n)

#==============================================================================
#Number of indicies for this process:
number_ind = len(rank_indices)
#Set up arrays to store all the varibles for all the timesteps
save_times = np.array([i[1].attrs['t'] for i in f.file.time_group.items()])

save_points = np.zeros([number_ind, len(the_line.output.points),3])
save_vpar = np.zeros([number_ind, len(the_line.output.points)])
save_vperp = np.zeros([number_ind, len(the_line.output.points)])
save_vphi = np.zeros([number_ind, len(the_line.output.points)])
save_index = np.zeros([number_ind, len(the_line.output.points)])


Fpar_line = np.zeros([number_ind, len(the_line.output.points)])
Fperp_line = np.zeros([number_ind, len(the_line.output.points)])
Fphi_line = np.zeros([number_ind, len(the_line.output.points)])

density_line = np.zeros([number_ind, len(the_line.output.points)])
va_line = np.zeros([number_ind, len(the_line.output.points)])
beta_line = np.zeros([number_ind, len(the_line.output.points)])
cs_line = np.zeros([number_ind, len(the_line.output.points)])

#==============================================================================

for i,n in enumerate(rank_indices):
    f.read_timestep(n)
    
    [bfield, vfield, density,
     valf, cs, beta] = util.process_next_step(f, cube_slice, bfield, vfield,
                                              density, valf, cs, beta)

    #Update surface
    ttf.update_flux_surface(surf_seeds[n], surf_field_lines, surface)
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

    the_line = ttf.update_the_line(the_line, surf_seeds[n].points, line_n, len(the_line.output.points))

    #line varibles
    surf_line_index, surf_line_points = ttf.get_surface_indexes(surface.output,the_line)
    
    
    save_index[i] = surf_line_index
    save_points[i] = surf_line_points
    save_vpar[i] = vpar[surf_line_index]
    save_vperp[i] = vperp[surf_line_index]
    save_vphi[i] = vphi[surf_line_index]

    step_poly = surface.output

    
    surface_density = ttf.update_interpolated_scalars(step_poly, surface_den_filter)
    surface_va = ttf.update_interpolated_scalars(step_poly, surface_va_filter)
    surface_beta = ttf.update_interpolated_scalars(step_poly, surface_beta_filter)
    surface_cs = ttf.update_interpolated_scalars(step_poly, surface_cs_filter)
    
    #Calculate Fluxs
    vpar = ttf.get_data(step_poly,'vpar') * 1e3 #stored in km/s
    Fpar = surface_density * vpar**2 * surface_cs
    
    vperp = ttf.get_data(step_poly,'vperp') * 1e3 #stored in km/s
    Fperp = surface_density * vperp**2 * np.sqrt(surface_cs**2 + surface_va**2)
    
    vphi = ttf.get_data(step_poly,'vphi') * 1e3 #stored in km/s
    Fphi = surface_density * vphi**2 * surface_va

    ttf.write_flux(path_join(
        "SurfaceFlux_AllVars_%s_%s_%s_%05i.vtp"%(driver, period, post_amp, n)),
                             step_poly, surface_density, surface_va,
                             surface_beta, surface_cs, Fpar, Fperp, Fphi)    
    
    Fpar_line[i] = Fpar[surf_line_index]
    Fperp_line[i] = Fperp[surf_line_index]
    Fphi_line[i] = Fphi[surf_line_index]
    
    density_line[i] = surface_density[surf_line_index]
    va_line[i] = surface_va[surf_line_index]
    beta_line[i] = surface_beta[surf_line_index]
    cs_line[i] = surface_cs[surf_line_index]
    
    print "Done Flux %i"%n

#==============================================================================
# Final rank 0 stuff
#==============================================================================
#Gather the data
save_points_r0 = comm.gather(save_points, root=0)
save_vpar_r0 = comm.gather(save_vpar, root=0)
save_vperp_r0 = comm.gather(save_vperp, root=0)
save_vphi_r0 = comm.gather(save_vphi, root=0)
save_index_r0= comm.gather(save_index, root=0)

Fpar_line_r0 = comm.gather(Fpar_line, root=0)
Fperp_line_r0 = comm.gather(Fperp_line, root=0)
Fphi_line_r0 = comm.gather(Fphi_line, root=0)

density_line_r0 = comm.gather(density_line, root=0)
va_line_r0 = comm.gather(va_line, root=0)
beta_line_r0 = comm.gather(beta_line, root=0)
cs_line_r0 = comm.gather(cs_line, root=0)

if rank == 0:
    save_points = np.concatenate(save_points_r0,axis=0)
    save_vpar = np.concatenate(save_vpar_r0,axis=0)
    save_vperp = np.concatenate(save_vperp_r0,axis=0)
    save_vphi = np.concatenate(save_vphi_r0,axis=0)
    save_index = np.concatenate(save_index_r0,axis=0)
    
    Fpar_line = np.concatenate(Fpar_line_r0,axis=0)
    Fperp_line = np.concatenate(Fperp_line_r0,axis=0)
    Fphi_line = np.concatenate(Fphi_line_r0,axis=0)
    
    density_line = np.concatenate(density_line_r0,axis=0)
    va_line = np.concatenate(va_line_r0,axis=0)
    beta_line = np.concatenate(beta_line_r0,axis=0)
    cs_line = np.concatenate(cs_line_r0,axis=0)

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