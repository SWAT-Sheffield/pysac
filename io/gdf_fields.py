# -*- coding: utf-8 -*-
"""
This file defines a dictionary that holds the gdf parameters for a 3D MHD variable set in SAC
"""

sac_gdf_output = {
      'density_pert':{
                      'field_title': 'density_pert',
                      'field_name': 'pertubation density',
                      'field_units': 'kg m^{-3}',
                      'field_to_cgs': 1e-3,
                      'staggering': 0
                      },
      'velocity_x':{
                      'field_title': 'velocity_x',
                      'field_name': 'x component velocity',
                      'field_units': 'ms^{-1}',
                      'field_to_cgs': 100,
                      'staggering': 0
                      },
      'velocity_y':{
                      'field_title': 'velocity_y',
                      'field_name': 'y component velocity',
                      'field_units': 'ms^{-1}',
                      'field_to_cgs': 100,
                      'staggering': 0
                      },
      'velocity_z':{
                      'field_title': 'velocity_z',
                      'field_name': 'z component velocity',
                      'field_units': 'ms^{-1}',
                      'field_to_cgs': 100,
                      'staggering': 0
                      },
      'internal_energy_pert':{
                      'field_title': 'internal_energy_pert',
                      'field_name': 'pertubation internal energy',
                      'field_units': 'Pa',
                      'field_to_cgs': 10,
                      'staggering': 0
                      },
      'mag_field_x_pert':{
                      'field_title': 'mag_field_x_pert',
                      'field_name': 'pertubation magnetic field x component',
                      'field_units': 'T',
                      'field_to_cgs': 1e4,
                      'staggering': 0
                      },
      'mag_field_y_pert':{
                      'field_title': 'mag_field_y_pert',
                      'field_name': 'pertubation magnetic field y component',
                      'field_units': 'T',
                      'field_to_cgs': 1e4,
                      'staggering': 0
                      },
      'mag_field_z_pert':{
                      'field_title': 'mag_field_z_pert',
                      'field_name': 'pertubation magnetic field z component',
                      'field_units': 'T',
                      'field_to_cgs': 1e4,
                      'staggering': 0
                      },
      'internal_energy_bg':{
                      'field_title': 'internal_energy_bg',
                      'field_name': 'background internal energy',
                      'field_units': 'Pa',
                      'field_to_cgs': 10,
                      'staggering': 0
                      },
      'density_bg':{
                      'field_title': 'density_bg',
                      'field_name': 'background density',
                      'field_units': 'kg m^{-3}',
                      'field_to_cgs': 1e-3,
                      'staggering': 0
                      },
      'mag_field_x_bg':{
                      'field_title': 'mag_field_x_bg',
                      'field_name': 'background magnetic field x component',
                      'field_units': 'T',
                      'field_to_cgs': 1e4,
                      'staggering': 0
                      },
      'mag_field_y_bg':{
                      'field_title': 'mag_field_y_bg',
                      'field_name': 'background magnetic field y component',
                      'field_units': 'T',
                      'field_to_cgs': 1e4,
                      'staggering': 0
                      },
      'mag_field_z_bg':{
                      'field_title': 'mag_field_z_bg',
                      'field_name': 'background magnetic field z component',
                      'field_units': 'T',
                      'field_to_cgs': 1e4,
                      'staggering': 0
                      },
     }
