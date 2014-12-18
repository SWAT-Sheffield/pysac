# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 17:56:58 2014

@author: sm1fg
"""
import numpy as np
import pysac.io.gdf_writer as gdf
import h5py
import astropy.units as u
##============================================================================
## Save a file!!!
##============================================================================
""" For large data arrays this has been producing overfull memory so moved to 
dedicated serial function, which should be parallel and moved ahead of gather
for plotting -- maybe handle plotting from hdf5 also

This file is potentially large - recommended to mkdir hdf5 in /data/${USER}
and add symlink to ${HOME} to avoid exceeding quota.
"""
def save_SACvariables(
              model,
              filename,
              rho,
              Bx,
              By,
              Bz,
              energy,
              logical_pars,
              physical_constants,
              scales,
              coords,
              Nxyz
             ):

    """ 
    """
    rank = 0
    if logical_pars['l_mpi']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        gather_vars = [
                       rho,Bx,By,Bz,energy
                      ]
        concat_vars = []
        for var in gather_vars:
            concat_vars.append(comm.gather(var, root=0))
    
        if rank == 0:
            out_vars = []
            for cvar in concat_vars:
                out_vars.append(np.concatenate(cvar, axis=0))
    
            rho,Bx,By,Bz,energy = out_vars
    if rank == 0:
    
        grid_dimensions = [[Nxyz[0], Nxyz[1], Nxyz[2]]]
        left_edge =  coords['xmin']*scales['length'],\
                     coords['ymin']*scales['length'],\
                     coords['zmin']*scales['length']
        right_edge = coords['xmax']*scales['length'],\
                     coords['ymax']*scales['length'],\
                     coords['zmax']*scales['length']
        g0_SI = physical_constants['gravity']*scales['velocity']/scales['time']

        dummy = np.zeros(rho.shape)
        simulation_parameters = gdf.SimulationParameters([
                            ['boundary_conditions', np.zeros(6) + 2],
                            ['cosmological_simulation', 0          ],
                            ['current_iteration', 0                ],
                            ['current_time', 0.0                   ],
                            ['dimensionality', 3                   ],
                            ['domain_dimensions', grid_dimensions  ],
                            ['domain_left_edge', left_edge         ],
                            ['domain_right_edge', right_edge       ],
                            ['eta', 0.0                            ],
                            ['field_ordering', 0                   ],
                            ['gamma', 1.66666667                   ],
                            ['gravity0', 0.0                       ],
                            ['gravity1', 0.0                       ],
                            ['gravity2', g0_SI                     ],
                            ['nu', 0.0                             ],
                            ['num_ghost_zones', 0                  ],
                            ['refine_by', 0                        ],
                            ['unique_identifier', 'sacgdf2014'     ]
                            ])
    
        gdf_file = gdf.create_file(h5py.File(filename,'w'), simulation_parameters, grid_dimensions)
    
        gdf.write_field_u(gdf_file, rho*scales['density']*u.Unit('kg/m^3'),
                              'density_bg',
                              'Background Density' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('kg/m^3'), 
                              'density_pert',
                              'Perturbation Density' 
                              )
        gdf.write_field_u(gdf_file, energy*scales['energy']*u.Unit('Pa'), 
                              'internal_energy_bg',
                              'Background Internal Energy' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('Pa'), 
                              'internal_energy_pert',
                              'Perturbation Internal Energy' 
                              )
        gdf.write_field_u(gdf_file, Bx*scales['magnetic']*u.Unit('Tesla'), 
                              'mag_field_x_bg',
                              'x Component of Background Magnetic Field' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('Tesla'), 
                              'mag_field_x_pert',
                              'x Component of Pertubation Magnetic Field' 
                              )
        gdf.write_field_u(gdf_file, By*scales['magnetic']*u.Unit('Tesla'), 
                              'mag_field_y_bg',
                              'y Component of Background Magnetic Field' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('Tesla'), 
                              'mag_field_y_pert',
                              'y Component of Pertubation Magnetic Field' 
                              )
        gdf.write_field_u(gdf_file, Bz*scales['magnetic']*u.Unit('Tesla'), 
                              'mag_field_z_bg',
                              'z Component of Background Magnetic Field' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('Tesla'), 
                              'mag_field_z_pert',
                              'z Component of Pertubation Magnetic Field' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('m/s'), 
                              'velocity_x',
                              'x Component of Velocity' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('m/s'), 
                              'velocity_y',
                              'y Component of Velocity' 
                              )
        gdf.write_field_u(gdf_file, dummy*u.Unit('m/s'), 
                                'velocity_z',
                                'z Component of Velocity' 
                                )
    
        gdf_file.close()

