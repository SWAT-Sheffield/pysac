"""
Routines for the writing of GDF files
"""

import copy

import numpy as np
import astropy.units as u
from astropy.utils import deprecated
import yt

import h5py
from h5py import h5s, h5p, h5fd

__all__ = ['write_field', 'write_field_raw', 'create_file', 'SimulationParameters']


class SimulationParameters(dict):
    def __init__(self, *args, **kwargs):
        self.REQUIRED_KEYS = ('refine_by', 'dimensionality', 'domain_dimensions',
                              'current_time', 'domain_left_edge', 'domain_right_edge',
                              'unique_identifier', 'cosmological_simulation',
                              'num_ghost_zones', 'field_ordering',
                              'boundary_conditions')

        for key in self.REQUIRED_KEYS:
            self[key] = None

        self['refine_by'] = 0
        self['cosmological_simulation'] = 0
        self['unique_identifier'] = 'sacgdf2014'
        super(SimulationParameters, self).__init__(*args, **kwargs)

    def __str__(self):

        head = """
Simulation Parameters Object
----------------------------

Required Attributes:

"""

        for key in self.REQUIRED_KEYS:
            head += key+': {{{0}}}\n'.format(key)

        others = copy.copy(self)
        for key in self.REQUIRED_KEYS:
            others.pop(key)

        if len(others):
            tail = """
Other Attributes:

"""
            for key in others:
                tail += key+': {{{0}}}\n'.format(key)

        else:
            tail = ''

        string = head.format(**self) + tail.format(**others)

        return string

    def __repr__(self):
        return str(self)


def create_file(f, simulation_parameters, grid_dimensions,
                data_author=None, data_comment=None):
    """
    Do all the structral creation of a gdf file.

    gdf files should be written in a x,y,z order, please swap them before
    calling this function!!

    Parameters
    ----------
    gdf_path: string or h5py instance
        Filename to save out

    simulation_parameters: dict
        Key value pairs for attributes to be written to the
        simulation_parameters group.

    grid_dimensions

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
    if isinstance(f, basestring):
        f = h5py.File(f, 'a')

    # "gridded_data_format" group
    g = f.create_group("gridded_data_format")
    g.attrs["data_software"] = "Sheffield Advanced Code"
    g.attrs["data_software_version"] = "pySAC0.2"
    if data_author is not None:
        g.attrs["data_author"] = data_author
    if data_comment is not None:
        g.attrs["data_comment"] = data_comment

    # "simulation_parameters" group
    g = f.create_group("simulation_parameters")
    for key, value in simulation_parameters.items():
        g.attrs[key] = value


    # "field_types" group
    g = f.create_group("field_types")

    # "particle_types" group
    g = f.create_group("particle_types")

    # "dataset_units" group
    g = f.create_group("dataset_units")
    BASE = [("length_unit", "m"), ("mass_unit", "kg"), ("time_unit", "s"),
            ("velocity_unit", "m / s"), ("magnetic_unit", "T")]
    for base_name, base_unit in BASE:
        f["dataset_units"][base_name] = 1.0
        f["dataset_units"][base_name].attrs["unit"] = base_unit

    # root datasets -- info about the grids
    f["grid_dimensions"] = np.reshape(grid_dimensions, (1, 3)) #needs to be 1XN
    f["grid_left_index"] = np.zeros((1,3)) #needs to be 1XN
    f["grid_level"] = np.zeros(1)
    # @todo: Fill with proper values
    f["grid_parent_id"] = np.zeros(1)
    f["grid_particle_count"] = np.zeros((1,1))

    # "data" group -- where we should spend the most time
    d = f.create_group("data")
    gr = d.create_group("grid_%010i"%0)

    return f

def write_field(gdf_file, data, field_title, field_name, field_shape=None,
                arr_slice=np.s_[:], staggering=0,
                collective=False, api='high'):
    """
    Write a field to an existing gdf file.

    Parameters
    ----------

    gdf_file: h5py.File
        Open, writeable gdf file

    data: astropy.units.Quantity
        The data to be written

    field_title: str
        The name of the field dataset

    field_name: str
        The long name for the field

    field_shape: (optional) tuple
        The shape of the whole dataset, if not specified use data.shape.

    arr_slice: (optional) np.s_
        The slice of the whole dataset to write

    staggering: (optional) int
        The 'staggering' of the gdf field

    """
    gr = gdf_file["/data/grid_%010i"%0]

    if not isinstance(data, u.Quantity):
        raise TypeError("data must be an astropy Quantity")

    # Make sure we are in SI land
    field = data.si

    if not field_shape:
        field_shape = field.shape

    dset = gr.create_dataset(field_title, field_shape, dtype='d')

    fv = gdf_file['field_types'].create_group(field_title)
    fv.attrs['field_name'] = field_name
    fv.attrs['field_to_cgs'] = field.unit.to_system(u.cgs)[0].scale
    fv.attrs['field_units'] = np.string_(field.unit.to_string("latex").strip('$'))
    fv.attrs['staggering'] = staggering

    ytunit = np.string_(yt.YTQuantity.from_astropy(field).units)
    gdf_file["dataset_units"][field_title] = 1.0
    gdf_file["dataset_units"][field_title].attrs["unit"] = ytunit

#    gr[field_title][arr_slice] = np.array(field)
    if api == 'high':
        _write_dset_high(dset, field, arr_slice, collective=collective)
    elif api == 'low':
        _write_dset_low(dset, field, arr_slice, collective=collective)
    else:
        raise ValueError("Please specifiy 'high' or 'low'")

@deprecated('0.3', "Writing GDF files from non-quantity arrays is no longer supported")
def write_field_raw(gdf_file, field, field_shape=None, arr_slice=np.s_[:],
                collective=False, api='high'):
    """
    Write a field to an existing gdf file.

    Parameters
    ----------

    gdf_file: h5py.File
        Open, writeable gdf file

    field: dict
        A dict containing the following keys:
            'field': ndarray
            'field_title': string
            'field_name': string
            'field_units': string
            'field_to_cgs': float
            'staggering': 0 or 1

    arr_slice: (optional) np.s_
        The slice of the whole dataset to write

    staggering: (optional) int
        The 'staggering' of the gdf field

    collective: bool
        Use MPI collective write

    api: str
        'high' or 'low' signifiyng the h5py API to use for write. Used for
        benchmarking.

    """
    gr = gdf_file["/data/grid_%010i"%0]

    if not field_shape:
        field_shape = field['field'].shape

    dset = gr.create_dataset(field['field_title'], field_shape, dtype='d')

    fv = gdf_file['field_types'].create_group(field['field_title'])
    fv.attrs['field_name'] = field['field_name']
    fv.attrs['field_to_cgs'] = field['field_to_cgs']
    fv.attrs['field_units'] = field['field_units']
    fv.attrs['staggering'] = field['staggering']

    if api == 'high':
        _write_dset_high(dset, field['field'], arr_slice, collective=collective)
    elif api == 'low':
        _write_dset_low(dset, field['field'], arr_slice, collective=collective)
    else:
        raise ValueError("Please specifiy 'high' or 'low'")


def _write_dset_high(dset, data, arr_slice,collective=False):
    if collective:
        with dset.collective:
            dset[arr_slice] = np.ascontiguousarray(data)
    else:
        dset[arr_slice] = np.ascontiguousarray(data)

def _write_dset_low(dset, data, arr_slice, collective=False):
    memory_space = h5s.create_simple(data.shape)
    file_space = dset.id.get_space()

    s = (arr_slice[0].start,arr_slice[1].start,arr_slice[2].start)
    e = (arr_slice[0].stop,arr_slice[1].stop,arr_slice[2].stop)

    count = tuple([ee - ss for ss,ee in zip(s,e)])

    file_space.select_hyperslab(s, count)

    if collective:
        dxpl = h5p.create(h5p.DATASET_XFER)
        dxpl.set_dxpl_mpio(h5fd.MPIO_COLLECTIVE)
    else:
        dxpl = None

    dset.id.write(memory_space, file_space,
            np.ascontiguousarray(data),dxpl=dxpl)
