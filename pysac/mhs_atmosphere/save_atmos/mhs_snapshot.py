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

    """ Save the background variables for a SAC model in hdf5 (gdf default)
    format after collating the data from mpi sub processes if necessary.
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

        grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
        left_edge =  np.array([coords['xmin']*scales['length'],
                               coords['ymin']*scales['length'],
                               coords['zmin']*scales['length']])
        right_edge = np.array([coords['xmax']*scales['length'],
                               coords['ymax']*scales['length'],
                               coords['zmax']*scales['length']])
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
        gdf.write_field_u(gdf_file,
                              energy*scales['energy density']*u.Unit('Pa'),
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
#============================================================================

def save_SACsources(
                    model,
                    sourcesfile,
                    Fx,
                    Fy,
                    logical_pars,
                    physical_constants,
                    scales,
                    coords,
                    Nxyz
                   ):
    """ Save the balancing forces for a SAC model with multiple flux tubes in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    rank = 0
    if logical_pars['l_mpi']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        gather_vars = [
                       Fx,Fy
                      ]
        concat_vars = []
        for var in gather_vars:
            concat_vars.append(comm.gather(var, root=0))

        if rank == 0:
            out_vars = []
            for cvar in concat_vars:
                out_vars.append(np.concatenate(cvar, axis=0))

            Fx,Fy = out_vars
    if rank == 0:

        grid_dimensions = [[Nxyz[0], Nxyz[1], Nxyz[2]]]
        left_edge =  coords['xmin']*scales['length'],\
                     coords['ymin']*scales['length'],\
                     coords['zmin']*scales['length']
        right_edge = coords['xmax']*scales['length'],\
                     coords['ymax']*scales['length'],\
                     coords['zmax']*scales['length']
        g0_SI = physical_constants['gravity']*scales['velocity']/scales['time']

        dummy = np.zeros(Fx.shape)
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

        gdf_file = gdf.create_file(h5py.File(sourcesfile,'w'), simulation_parameters, grid_dimensions)

        gdf.write_field_u(gdf_file,
                              Fx*scales['force density']*u.Unit('Newton/m^3'),
                              'balancing_force_x_bg',
                              'x Component of Background Balancing Force'
                              )
        gdf.write_field_u(gdf_file,
                              Fy*scales['force density']*u.Unit('Newton/m^3'),
                              'balancing_force_y_bg',
                              'y Component of Background Balancing Force'
                              )
#        gdf.write_field_u(gdf_file, rho*scales['density']*u.Unit('kg/m^3'),
#                              'density_bg',
#                              'Background Density'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('kg/m^3'),
#                              'density_pert',
#                              'Perturbation Density'
#                              )
#        gdf.write_field_u(gdf_file, energy*scales['energy density']*u.Unit('Pa'),
#                              'internal_energy_bg',
#                              'Background Internal Energy'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('Pa'),
#                              'internal_energy_pert',
#                              'Perturbation Internal Energy'
#                              )
#        gdf.write_field_u(gdf_file, Bx*scales['magnetic']*u.Unit('Tesla'),
#                              'mag_field_x_bg',
#                              'x Component of Background Magnetic Field'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('Tesla'),
#                              'mag_field_x_pert',
#                              'x Component of Pertubation Magnetic Field'
#                              )
#        gdf.write_field_u(gdf_file, By*scales['magnetic']*u.Unit('Tesla'),
#                              'mag_field_y_bg',
#                              'y Component of Background Magnetic Field'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('Tesla'),
#                              'mag_field_y_pert',
#                              'y Component of Pertubation Magnetic Field'
#                              )
#        gdf.write_field_u(gdf_file, Bz*scales['magnetic']*u.Unit('Tesla'),
#                              'mag_field_z_bg',
#                              'z Component of Background Magnetic Field'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('Tesla'),
#                              'mag_field_z_pert',
#                              'z Component of Pertubation Magnetic Field'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('m/s'),
#                              'velocity_x',
#                              'x Component of Velocity'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('m/s'),
#                              'velocity_y',
#                              'y Component of Velocity'
#                              )
#        gdf.write_field_u(gdf_file, dummy*u.Unit('m/s'),
#                                'velocity_z',
#                                'z Component of Velocity'
#                                )

        gdf_file.close()
#============================================================================

def save_auxilliary1D(
                    model,
                    auxfile,
                    pressure_m,
                    rho_m,
                    temperature,
                    pbeta,
                    alfven,
                    cspeed,
                    dxB2,
                    dyB2,
                    val,
                    mtw,
                    pressure_Z,
                    rho_Z,
                    Rgas_Z,
                    logical_pars,
                    physical_constants,
                    scales,
                    coords,
                    Nxyz
                   ):
    """ Save auxilliary variables for use in plotting background setup in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    rank = 0
    if logical_pars['l_mpi']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        gather_vars = [
                       pressure_m,
                       rho_m,
                       temperature,
                       pbeta,
                       alfven,
                       cspeed,
                       dxB2,
                       dyB2
                      ]
        concat_vars = []
        for var in gather_vars:
            concat_vars.append(comm.gather(var, root=0))

        if rank == 0:
            out_vars = []
            for cvar in concat_vars:
                out_vars.append(np.concatenate(cvar, axis=0))

            pressure_m, rho_m, temperature, pbeta,\
                alfven, cspeed, dxB2, dyB2 = out_vars
    if rank == 0:
#
        grid_dimensions = [[Nxyz[0], Nxyz[1], Nxyz[2]]]
        left_edge =  coords['xmin']*scales['length'],\
                     coords['ymin']*scales['length'],\
                     coords['zmin']*scales['length']
        right_edge = coords['xmax']*scales['length'],\
                     coords['ymax']*scales['length'],\
                     coords['zmax']*scales['length']
        g0_SI = physical_constants['gravity']*scales['velocity']/scales['time']

        dummy3D = np.zeros(pressure_m.shape)
        dummy1D = np.zeros(pressure_Z.shape)
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

        gdf_file = gdf.create_file(h5py.File(auxfile,'w'), simulation_parameters, grid_dimensions)

        gdf.write_field_u(gdf_file,
                              pressure_Z*scales['energy density']*u.Unit('Pa'),
                              '1D_plasma_pressure',
                              'Background pressure Z-profile'
                              )
        gdf.write_field_u(gdf_file,
                              rho_Z*scales['density']*u.Unit('kg/m^3'),
                              '1D_plasma_density',
                              'Background density Z-profile'
                              )
        gdf.write_field_u(gdf_file,
                              Rgas_Z*scales['velocity']**2/
                                     scales['temperature']*u.Unit('J/(kg*K)'),
                              '1D_ideal_gas_constant',
                              'Background R_gas Z-profile'
                              )
        gdf.write_field_u(gdf_file,
                              pressure_m*scales['energy density']*u.Unit('Pa'),
                              '3D_plasma_pressure_balance',
                              'Background magneto-pressure balance'
                              )
        gdf.write_field_u(gdf_file,
                              rho_m*scales['density']*u.Unit('kg/m^3'),
                              '3D_plasma_density_balance',
                              'Background magneto-density balance'
                              )
        gdf.write_field_u(gdf_file,
                              temperature*scales['temperature']*u.Unit('K'),
                              '3D_temperature',
                              'Background temperature'
                              )
        gdf.write_field_u(gdf_file,
                              pbeta*u.Unit(''),
                              '3D_plasma_beta',
                              'Background plasma beta'
                              )
        gdf.write_field_u(gdf_file,
                              alfven*scales['velocity']*u.Unit('m/s'),
                              '3D_Alfven_speed',
                              'Background Alfven speed'
                              )
        gdf.write_field_u(gdf_file,
                              cspeed*scales['velocity']*u.Unit('m/s'),
                              '3D_sound_speed',
                              'Background sound speed'
                              )
        gdf.write_field_u(gdf_file,
                              dxB2*scales['energy density']*u.Unit('Pa'),
                              '3D_mag_tension_x',
                              'x-component bg magnetic tension'
                              )
        gdf.write_field_u(gdf_file,
                              dyB2*scales['energy density']*u.Unit('Pa'),
                              '3D_mag_tension_y',
                              'y-component bg magnetic tension'
                              )
        gdf.write_field_u(gdf_file,
                              val[:,0]*u.Unit('m'),
                              'val3c_z',
                              'Height of VAL data'
                              )
        gdf.write_field_u(gdf_file,
                              val[:,1]*u.Unit('kg/m^3'),
                              'val3c_density',
                              'VAL plasma density data'
                              )
        gdf.write_field_u(gdf_file,
                              val[:,2]*u.Unit('Pa'),
                              'val3c_pressure',
                              'VAL plasma pressure data'
                              )
        gdf.write_field_u(gdf_file,
                              val[:,3]*u.Unit('K'),
                              'val3c_temperature',
                              'VAL temperature data'
                              )
        gdf.write_field_u(gdf_file,
                              mtw[:,0]*u.Unit('m'),
                              'mtw_z',
                              'Height of MTW data'
                              )
        gdf.write_field_u(gdf_file,
                              mtw[:,2]*u.Unit('Pa'),
                              'mtw_pressure',
                              'MTW plasma pressure data'
                              )
        gdf.write_field_u(gdf_file,
                              mtw[:,1]*u.Unit('K'),
                              'mtw_temperature',
                              'MTW temperature data'
                              )


        gdf_file.close()
