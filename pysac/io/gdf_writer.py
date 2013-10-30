# -*- coding: utf-8 -*-

import astropy.units as u
import numpy as np
import h5py

def convert_w_3D(w, w_):
    """
    This function converts a 4D w array to a field dictionary.
    
    It is hardcoded to expect the following varibles:
    
    varnames: ['h', 'm1', 'm2', 'm3', 'e', 'b1', 'b2', 'b3', 'eb', 'rhob', 'bg1', 'bg2', 'bg3']
    where magnticfield is in NON-scaled units of B i.e. Tesla
    and everything else is also in base SI, m, Pa, kg /m^3 etc.
    
    Paramters
    ---------
    w: np.ndarray
        The w array in Si units
    
    w_: dict
        A conversion between varnames and w indices
    """
    
    rhot = w[w_['h']] + w[w_['rhob']]
    
    sac_gdf_output = {}
    
    sac_gdf_output['density_pert'] = {
                    'field': w[w_['h']] * u.Unit('kg/m^3'),
                    'field_name': 'pertubation density',       
                    'staggering': 0
                    }
    sac_gdf_output['density_bg'] = {
                    'field': w[w_['rhob']] * u.Unit('kg/m^3'),
                    'field_name':'background density',       
                    'staggering':0
                    }
    sac_gdf_output['velocity_x'] = {
                    'field': w[w_['m2']]/rhot * u.Unit('m/s'),
                    'field_name': 'velocity x component',       
                    'staggering': 0
                    }
    sac_gdf_output['velocity_y'] = {
                    'field': w[w_['m3']]/rhot * u.Unit('m/s'),
                    'field_name': 'velocity y component',     
                    'staggering': 0
                    }
    sac_gdf_output['velocity_z'] = {
                    'field': w[w_['m1']]/rhot * u.Unit('m/s'),
                    'field_name':'velocity z component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_x_pert'] = {
                    'field':u.Quantity(w[w_['b2']], unit=u.T),
                    'field_name':'pertubation magnetic field x component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_y_pert'] = {
                    'field': w[w_['b3']] * u.T,
                    'field_name':'pertubation magnetic field y component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_z_pert'] = {
                    'field': w[w_['b1']] * u.T,
                    'field_name': 'pertubation magnetic field z component',        
                    'staggering': 0
                    }
    sac_gdf_output['mag_field_x_bg'] = {
                    'field': w[w_['bg2']] * u.T,
                    'field_name':'background magnetic field x component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_y_bg'] = {
                    'field': w[w_['bg3']] * u.T,
                    'field_name': 'background magnetic field y component',        
                    'staggering': 0
                    }
    sac_gdf_output['mag_field_z_bg'] = {
                    'field': w[w_['bg1']] * u.T,
                    'field_name': 'background magnetic field z component',        
                    'staggering': 0
                    }
    sac_gdf_output['internal_energy_pert'] = {
                    'field': w[w_['e']] * u.Pa,
                    'field_name':'pertubation internal energy',        
                    'staggering':0
                    }
    sac_gdf_output['internal_energy_bg'] = {
                    'field': w[w_['eb']] * u.Pa,
                    'field_name':'background internal energy',        
                    'staggering':0
                    }

    return sac_gdf_output

def convert_w_2D(w, w_):
    """
    This function converts a 3D w array to a field dictionary.
    
    It is hardcoded to expect the following varibles:
    
    varnames: ['h', 'm1', 'm2', 'e', 'b1', 'b2', 'eb', 'rhob', 'bg1', 'bg2']
    where magnticfield is in NON-scaled units of B i.e. Tesla
    and everything else is also in base SI, m, Pa, kg /m^3 etc.
    
    Paramters
    ---------
    w: np.ndarray
        The w array in Si units
    
    w_: dict
        A conversion between varnames and w indices
    """
    
    rhot = w[w_['h']] + w[w_['rhob']]
    
    sac_gdf_output = {}
    
    sac_gdf_output['density_pert'] = {
                    'field': w[w_['h']] * u.Unit('kg/m^3'),
                    'field_name': 'pertubation density',       
                    'staggering': 0
                    }
    sac_gdf_output['density_bg'] = {
                    'field': w[w_['rhob']] * u.Unit('kg/m^3'),
                    'field_name':'background density',       
                    'staggering':0
                    }
    sac_gdf_output['velocity_x'] = {
                    'field': w[w_['m2']]/rhot * u.Unit('m/s'),
                    'field_name': 'velocity x component',       
                    'staggering': 0
                    }
    sac_gdf_output['velocity_z'] = {
                    'field': w[w_['m1']]/rhot * u.Unit('m/s'),
                    'field_name':'velocity z component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_x_pert'] = {
                    'field':u.Quantity(w[w_['b2']], unit=u.T),
                    'field_name':'pertubation magnetic field x component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_z_pert'] = {
                    'field': w[w_['b1']] * u.T,
                    'field_name': 'pertubation magnetic field z component',        
                    'staggering': 0
                    }
    sac_gdf_output['mag_field_x_bg'] = {
                    'field': w[w_['bg2']] * u.T,
                    'field_name':'background magnetic field x component',        
                    'staggering':0
                    }
    sac_gdf_output['mag_field_z_bg'] = {
                    'field': w[w_['bg1']] * u.T,
                    'field_name': 'background magnetic field z component',        
                    'staggering': 0
                    }
    sac_gdf_output['internal_energy_pert'] = {
                    'field': w[w_['e']] * u.Pa,
                    'field_name':'pertubation internal energy',        
                    'staggering':0
                    }
    sac_gdf_output['internal_energy_bg'] = {
                    'field': w[w_['eb']] * u.Pa,
                    'field_name':'background internal energy',        
                    'staggering':0
                    }

    return sac_gdf_output

def write_gdf(gdf_path, header, x, fields, arr_slice=np.s_[:],
              data_author=None, data_comment=None):
    """
    Write a gdf file from a vac-data header a w array and an x array
    
    gdf files should be written in a x,y,z order, please swap them before 
    calling this function!!
    
    Parameters
    ----------
    gdf_path: string or h5py instance
        Filename to save out
    
    header: dict
        A 'VACdata' like header
    
    x: astropy.units.Quantity
        The x array described in header, in x,y,z order
    
    fields: dict
        A dictionary with a key for each varname in header, and a value containing
        a dictionary with each of 
        {field: astropy.units.Quantity,
             The array for this field, in x,y,z (should be in SI not code units)
        field_name: string,
             A descriptive name for the field, spaces allowed
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
    if isinstance(gdf_path, h5py.File):
        f = gdf_path
    else:
        f = h5py.File(gdf_path, "w")
        
    domain_left_edge = [x[0][0,0,0].to(u.cm).value,
                        x[1][0,0,0].to(u.cm).value, 
                        x[2][0,0,0].to(u.cm).value]
    domain_right_edge = [x[0][-1,-1,-1].to(u.cm).value,
                         x[1][-1,-1,-1].to(u.cm).value,
                         x[2][-1,-1,-1].to(u.cm).value]

    create_file(f, header, domain_left_edge=domain_left_edge,
                domain_right_edge=domain_right_edge, data_author=data_author, 
                data_comment=data_comment)
    
    #Write the x array
    f['x'] = x

    for field_title,afield in fields.items():
       write_field(f, afield['field'], field_title, afield['field_name'],
                   field_shape=header['nx'], arr_slice=arr_slice, 
                   staggering=afield['staggering'])
    
    f.close()

def create_file(f, header, domain_left_edge=[], domain_right_edge=[],
                data_author="", data_comment=""):
    """
    Do all the structral creation of a gdf file
    
    gdf files should be written in a x,y,z order, please swap them before 
    calling this function!!
    
    Parameters
    ----------
    gdf_path: string or h5py instance
        Filename to save out
    
    header: dict
        A 'VACdata' like header
    
    domain_left_edge
    
    domain_right_edge
    
    data_author: (optional) string
        Author to write to file
    
    data_comment: (optional) string
        A comment to write to file
    
    Returns
    -------
    h5py.File instance
    
    Notes
    -----
    GDF is defined here: https://bitbucket.org/yt_analysis/grid_data_format/
    """
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
    g.attrs["domain_left_edge"] = domain_left_edge
    g.attrs["domain_right_edge"] = domain_right_edge
    g.attrs["unique_identifier"] = np.random.randint(1e15)
    g.attrs["cosmological_simulation"] = 0
    #TODO: hmmmm
    g.attrs["num_ghost_zones"] = 0
    g.attrs["field_ordering"] = 0
    # @todo: not yet supported by yt.
    g.attrs["boundary_conditions"] = np.array([0, 0, 0, 0, 0, 0], 'int32')
    
    #Write some VAC info just in case
    if 'final t' in header.keys():
        g.attrs['simulation_run_time'] = header['final t']
    g.attrs["eqpar"] = header['eqpar']
    index = next((i for i in xrange(len(header['varnames']))
                    if not(header['varnames'][i] in ["x","y","z"])),header['ndim'])
    g.attrs["eqpar_names"] = header['varnames'][index+header['nw']-1:]
    
    # "field_types" group
    g = f.create_group("field_types")

    # "particle_types" group
    g = f.create_group("particle_types")

    # root datasets -- info about the grids
    f["grid_dimensions"] = np.array([header['nx']]) #needs to be 1XN
    f["grid_left_index"] = np.array([[0]*header['ndim']]) #needs to be 1XN
    f["grid_level"] = np.zeros(1)
    # @todo: Fill with proper values
    f["grid_parent_id"] = np.zeros(1)
    f["grid_particle_count"] = np.zeros((1,1))

    # "data" group -- where we should spend the most time
    d = f.create_group("data")
    gr = d.create_group("grid_%010i"%0)
    
    return f

def write_field(gdf_file, data, field_title, field_name, field_shape=None, arr_slice=np.s_[:], staggering=0):
    gr = gdf_file["/data/grid_%010i"%0]
    field = data.si
    
    if not field_shape:
        field_shape = field.shape
        
    gr.create_dataset(field_title, field_shape, dtype='d')
    
    fv = gdf_file['field_types'].create_group(field_title)
    fv.attrs['field_name'] = field_name
    fv.attrs['field_to_cgs'] = field.unit.to_system(u.cgs)[0].scale
    fv.attrs['field_units'] = np.string_(field.unit.to_string("latex").strip('$'))
    fv.attrs['staggering'] = staggering
   
    gr[field_title][arr_slice] = np.array(field)