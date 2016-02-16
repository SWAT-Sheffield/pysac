# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 17:56:58 2014

@author: sm1fg
"""
import numpy as np
import pysac.io.gdf_writer as gdf
import h5py
import astropy.units as u
#import h5py as h5

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
              filename,
              rho,
              Bx,
              By,
              Bz,
              energy,
              option_pars,
              physical_constants,
              coords,
              Nxyz
             ):

    """ Save the background variables for a SAC model in hdf5 (gdf default)
    format after collating the data from mpi sub processes if necessary.
    """
    rank = 0
    if option_pars['l_mpi']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        gather_vars = [
                       rho.to(u.Unit('kg/m**3')),
                       Bx.to(u.Unit('T')),
                       By.to(u.Unit('T')),
                       Bz.to(u.Unit('T')),
                       energy.to(u.Unit('kg/(m s**2)'))
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
        print'writing',filename
        print'SAC background atmosphere'

        grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
        left_edge =  u.Quantity([coords['xmin'],
                                 coords['ymin'],
                                 coords['zmin']]).to(u.m)
        right_edge = u.Quantity([coords['xmax'],
                                 coords['ymax'],
                                 coords['zmax']]).to(u.m)
        g0 = physical_constants['gravity']
        gamma = 5./3.
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
                            ['gamma', gamma                   ],
                            ['gravity0', 0.0                       ],
                            ['gravity1', 0.0                       ],
                            ['gravity2', g0                        ],
                            ['nu', 0.0                             ],
                            ['num_ghost_zones', 0                  ],
                            ['refine_by', 0                        ],
                            ['unique_identifier', 'sacgdf2014'     ]
                            ])

        gdf_file = gdf.create_file(h5py.File(filename,'w'), simulation_parameters, grid_dimensions)

        gdf.write_field(gdf_file, rho,
                              'density_bg',
                              'Background Density'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('kg/m**3'),
                              'density_pert',
                              'Perturbation Density'
                              )
        gdf.write_field(gdf_file,
                              energy,
                              'internal_energy_bg',
                              'Background Internal Energy'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('kg/(m s**2)'),
                              'internal_energy_pert',
                              'Perturbation Internal Energy'
                              )
        gdf.write_field(gdf_file, Bx,
                              'mag_field_x_bg',
                              'x Component of Background Magnetic Field'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('T'),
                              'mag_field_x_pert',
                              'x Component of Pertubation Magnetic Field'
                              )
        gdf.write_field(gdf_file, By,
                              'mag_field_y_bg',
                              'y Component of Background Magnetic Field'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('T'),
                              'mag_field_y_pert',
                              'y Component of Pertubation Magnetic Field'
                              )
        gdf.write_field(gdf_file, Bz,
                              'mag_field_z_bg',
                              'z Component of Background Magnetic Field'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('T'),
                              'mag_field_z_pert',
                              'z Component of Pertubation Magnetic Field'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('m/s'),
                              'velocity_x',
                              'x Component of Velocity'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('m/s'),
                              'velocity_y',
                              'y Component of Velocity'
                              )
        gdf.write_field(gdf_file, dummy*u.Unit('m/s'),
                                'velocity_z',
                                'z Component of Velocity'
                              )

        gdf_file.close()


#============================================================================

def save_SACsources(
                    sourcesfile,
                    Fx,
                    Fy,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz
                   ):
    """ Save the balancing forces for a SAC model with multiple flux tubes in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    rank = 0
    if option_pars['l_mpi']:
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
        print'writing',sourcesfile
        print'SAC background source terms'

        grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
        left_edge =  u.Quantity([coords['xmin'],
                                 coords['ymin'],
                                 coords['zmin']]).to(u.m)
        right_edge = u.Quantity([coords['xmax'],
                                 coords['ymax'],
                                 coords['zmax']]).to(u.m)
        g0 = physical_constants['gravity']

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
                            ['gamma', physical_constants['gamma']  ],
                            ['gravity0', 0.0                       ],
                            ['gravity1', 0.0                       ],
                            ['gravity2', physical_constants['gravity']],
                            ['nu', 0.0                             ],
                            ['num_ghost_zones', 0                  ],
                            ['refine_by', 0                        ],
                            ['unique_identifier', 'sacgdf2014'     ]
                            ])

        gdf_file = gdf.create_file(h5py.File(sourcesfile,'w'), simulation_parameters, grid_dimensions)

        gdf.write_field(gdf_file,
                              Fx,
                              'balancing_force_x_bg',
                              'x Component of Background Balancing Force'
                              )
        gdf.write_field(gdf_file,
                              Fy,
                              'balancing_force_y_bg',
                              'y Component of Background Balancing Force'
                              )
        gdf_file.close()

#============================================================================

def save_auxilliary3D(
                    auxfile,
                    pressure_m,
                    rho_m,
                    temperature,
                    pbeta,
                    alfven,
                    cspeed,
                    dxB2,
                    dyB2,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz
                   ):
    """ Save auxilliary variables for use in plotting background setup in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    rank = 0
    if option_pars['l_mpi']:
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
        print'writing',auxfile
        print'non-SAC 3D auxilliary data for plotting'

        grid_dimensions = [Nxyz[0], Nxyz[1], Nxyz[2]]
        left_edge =  u.Quantity([coords['xmin'],
                                 coords['ymin'],
                                 coords['zmin']]).to(u.m)
        right_edge = u.Quantity([coords['xmax'],
                                 coords['ymax'],
                                 coords['zmax']]).to(u.m)
        

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
                            ['gamma', physical_constants['gamma']  ],
                            ['gravity0', 0.0                       ],
                            ['gravity1', 0.0                       ],
                            ['gravity2', physical_constants['gravity']],
                            ['nu', 0.0                             ],
                            ['num_ghost_zones', 0                  ],
                            ['refine_by', 0                        ],
                            ['unique_identifier', 'sacgdf2014'     ]
                            ])

        gdf_file = gdf.create_file(h5py.File(auxfile,'w'), simulation_parameters, grid_dimensions)

        gdf.write_field(gdf_file,
                              pressure_m,
                              'pressure_mhs',
                              'Background magneto-pressure balance'
                              )
        gdf.write_field(gdf_file,
                              rho_m,
                              'density_mhs',
                              'Background magneto-density balance'
                              )
        gdf.write_field(gdf_file,
                              temperature.to(u.K),
                              'temperature',
                              'Background temperature'
                              )
        gdf.write_field(gdf_file,
                              pbeta,
                              'plasma_beta',
                              'Background plasma beta'
                              )
        gdf.write_field(gdf_file,
                              alfven,
                              'alfven_speed',
                              'Background Alfven speed'
                              )
        gdf.write_field(gdf_file,
                              cspeed,
                              'sound_speed',
                              'Background sound speed'
                              )
        gdf.write_field(gdf_file,
                              dxB2,
                              'mag_tension_x',
                              'x-component background magnetic tension'
                              )
        gdf.write_field(gdf_file,
                              dyB2,
                              'mag_tension_y',
                              'y-component background magnetic tension'
                              )
        gdf_file.close()

#============================================================================

def save_auxilliary1D(
                    auxfile,
                    pressure_Z,
                    rho_Z,
                    Rgas_Z,
                    option_pars,
                    physical_constants,
                    coords,
                    Nxyz
                   ):
    """ Save auxilliary variables for use in plotting background setup in
    hdf5 (gdf default) format after collating the data from mpi sub processes
    if necessary.
    """
    rank = 0
    if option_pars['l_mpi']:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

    if rank == 0:
        print'writing',auxfile
        print'non-SAC 1D auxilliary data for plotting'

        grid_dimensions = [2, 2, Nxyz[2]] #dims > 1 to be read by yt
        left_edge =  u.Quantity([coords['xmin'],
                                 coords['ymin'],
                                 coords['zmin']]).to(u.m)
        right_edge = u.Quantity([coords['xmax'],
                                 coords['ymax'],
                                 coords['zmax']]).to(u.m)

        pressureHS = u.Quantity(np.zeros(grid_dimensions),
                                unit=pressure_Z.to("Pa").unit)
        rhoHS = u.Quantity(np.zeros(grid_dimensions),
                           unit=rho_Z.to('kg/m**3').unit)
        RgasHS = u.Quantity(np.zeros(grid_dimensions),
                            unit=Rgas_Z.to('m**2/(K s**2)').unit)
        pressureHS[:] = pressure_Z
        rhoHS[:] = rho_Z
        RgasHS[:] = Rgas_Z
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
                            ['gamma', physical_constants['gamma']  ],
                            ['gravity0', 0.0                       ],
                            ['gravity1', 0.0                       ],
                            ['gravity2', physical_constants['gravity']],
                            ['nu', 0.0                             ],
                            ['num_ghost_zones', 0                  ],
                            ['refine_by', 0                        ],
                            ['unique_identifier', 'sacgdf2014'     ]
                            ])

#        import pdb; pdb.set_trace()
        gdf_file = gdf.create_file(h5py.File(auxfile,'w'), simulation_parameters, grid_dimensions)

        gdf.write_field(gdf_file,
                              pressureHS,
                              'pressure_HS',
                              'Background 1D hydrostatic-pressure'
                              )
        gdf.write_field(gdf_file,
                              rhoHS,
                              'density_HS',
                              'Background 1D hydrostatic-density'
                              )
        gdf.write_field(gdf_file,
                              RgasHS,
                              'ideal_gas_constant_HS',
                              'Background 1D hydrostatic-R_gas'
                              )
        gdf_file.close()
    

#=============================================================================
