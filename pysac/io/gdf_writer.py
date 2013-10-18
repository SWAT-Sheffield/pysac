# -*- coding: utf-8 -*-

import numpy as np
import h5py

def write_gdf(gdf_path, header, x, fields, data_author=None, data_comment=None):
    """
    Write a gdf file from a vac-data header a w array and an x array
    
    gdf files should be written in a x,y,z order, please swap them before 
    calling this function!!
    
    Parameters
    ----------
    gdf_path: string
        Filename to save out
    
    header: dict
        A 'VACdata' like header
    
    x: np.ndarray
        The x array described in header, in x,y,z order
    
    fields: dict
        A dictionary with a key for each varname in header, and a value containing
        a dictionary with each of 
        {field: np.ndarray,
             The array for this field, in x,y,z (should be in SI not code units)
        field_value: string
             The save name for the field in yt nomencleture i.e. denisty_pert
        field_name: string,
             A descriptive name for the field, spaces allowed
        field_to_cgs: float,
             The conversion factor to cgs units
        field_units: string,
             A unit string
        staggering: integer
             0: cell-centred, 1: face-centred, 2: vertex-centred}
    
    data_author: (optional) string
        Author to write to file
    
    data_comment: (optional) string
        A comment to write to file
    
    Notes
    -----
    GDF is defined here: https://bitbucket.org/yt_analysis/grid_data_format/
    """
    
    # Create and open the file with h5py
    f = h5py.File(gdf_path, "w")

    # "gridded_data_format" group
    g = f.create_group("gridded_data_format")
    g.attrs["data_software"] = "SAC"
    if data_author is not None:
        g.attrs["data_author"] = data_author
    if data_comment is not None:
        g.attrs["data_comment"] = data_comment

    # "simulation_parameters" group
    g = f.create_group("simulation_parameters")
    g.attrs["refine_by"] = 0
    g.attrs["dimensionality"] = header['ndim']
    g.attrs["domain_dimensions"] = header['nx']
    g.attrs["current_time"] = header['t']
    g.attrs["domain_left_edge"] = [x[0][0,0,0], x[1][0,0,0], x[2][0,0,0]]
    g.attrs["domain_right_edge"] = [x[0][-1,-1,-1], x[1][-1,-1,-1], x[2][-1,-1,-1]]
    g.attrs["unique_identifier"] = np.random.randint(1e15)
    g.attrs["cosmological_simulation"] = 0
    #TODO: hmmmm
    g.attrs["num_ghost_zones"] = 0
    g.attrs["field_ordering"] = 0
    # @todo: not yet supported by yt.
    g.attrs["boundary_conditions"] = np.array([0, 0, 0, 0, 0, 0], 'int32')

    # "field_types" group
    g = f.create_group("field_types")

    # "particle_types" group
    g = f.create_group("particle_types")

    # root datasets -- info about the grids
    f["grid_dimensions"] = header['nx']
    f["grid_left_index"] = [0]*header['nx']
    f["grid_level"] = 0
    # @todo: Fill with proper values
    f["grid_parent_id"] = 0
    f["grid_particle_count"] = 0
    
    # "data" group -- where we should spend the most time
    d = f.create_group("data")
    
    gr = f.create_group("grid_%10i"%0)
    for afield in fields:
        gr[afield['field_value']] = afield['field']
        
        fv = f['field_types'].create_group(afield['field_value'])
        fv.attrs['field_name'] = afield['field_name']
        fv.attrs['field_to_cgs'] = afield['field_to_cgs']
        fv.arrts['field_units'] = afield['field_units']
        fv.attrs['staggering'] = afield['staggering']
    
    f.close()