# -*- coding: utf-8 -*-
"""
:Created on: Tue Oct  2 12:25:10 2012

:author: Stuart Mumford
"""

import numpy as np
from tvtk.api import tvtk
import vtk_read_functions as vrf

def vector_field(u, v, w, transpose=True, copy_data=True, field_name=None):
    """
    Take 3 numpy arrays and return a tvtk ImageData object
    
    It is unknown if this will work for non-cube data
    
    Parameters
    ----------
    u, v, w: np.ndarray
        Components of a vector field
    
    transpose: bool
        Transpose the arrays before conversion. Normally needed to convert from
        numpy's z,y,x order to vtk's x,y,z
    
    copy_data: bool
        Copy the numpy data in memory. If False, the numpy array must not be 
        deleted while the ImageData object exists, as it is using the pointer 
        to the numpy.ndata array.
    
    Returns
    -------
    vfield: tvtk.ImageData
    """
    vector_data = np.c_[u.ravel(), w.ravel(),  v.ravel()].ravel()
    data = vector_data.reshape([u.shape[0] , w.shape[1], v.shape[2], 3])
    
    imd = tvtk.ImageData()
    dims = list(data.shape)
    if len(dims) == 3:
        dims.insert(2, 1)
        data = np.reshape(data, dims)
    imd.dimensions = tuple(dims[:-1])
    imd.extent = 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1
    imd.update_extent = 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1
    sz = np.size(data)
    if transpose:
        data_t = np.transpose(data, (2, 1, 0, 3))
    else:
        data_t = data
    imd.point_data.vectors = np.reshape(data_t, (sz/3, 3))
    imd.point_data.vectors.name = field_name
    imd.update()
    
    return imd

def scalar_field(u, ravel=True, copy_data=True):
    """
    Take a numpy array and return a tvtk ImageData object
    
    Parameters
    ----------
    u: np.ndarray
        A scalar field
    
    ravel: bool
        Re-order the numpy array to be FORTRAN contiguious to be in x,y,z 
        order.
    
    copy_data: bool
        Copy the numpy data in memory. If False, the numpy array must not be 
        deleted while the ImageData object exists, as it is using the pointer 
        to the numpy.ndata array.
    
    Returns
    -------
    sdata: tvtk.ImageData
    """
    if ravel:
        sd = u.ravel().reshape(u.shape,order='F')
    else:
        sd = u
    
    sdata = vrf.array_to_vtk_image(sd, copy_data, data_type='scalar')
    
    return tvtk.to_tvtk(sdata)

def move_seeds(seeds,vfield,dt):
    """ Move a list of seeds based on a velocity field
    
    TODO:
    WARNING: THIS IS HARD CODED FOR GRID SIZE!
    
    Parameters
    ----------
    seeds: tvtk.PolyData
        Old seed points
    
    vfield: mayavi.sources.array_source.ArraySource object
        The velocity field
    
    dt: float
        The time step betweent the current and the previous step.
    
    Returns
    -------
    
    seeds_arr: ndarray
        New Seed points    
    """
    v_seed = tvtk.ProbeFilter(source=vfield,input=seeds)
    v_seed.update()
    int_vels = np.array(v_seed.output.point_data.vectors)[:,:2]/15.625
    seed_arr = np.array(seeds.points)
    seed_arr[:,:2] += int_vels * dt
    
    #seeds.points = seed_arr
    return seed_arr

def make_circle_seeds(n,r,**domain):
    """ Generate an array of n seeds evenly spaced in a circle at radius r
    
    Parameters
    ----------
    n:  integer
        Number of Seeds to Create
    
    r:  float
        Radius of the Circle in grid points
    
    **domain: Dict
        kwargs specifiying the properties of the domain.
    
    Returns
    -------
    surf_seeds: tvtk.PolyData
        vtkPolyData containg point data with the seed locations.
        Needs: xmax, ymax, zmax
    """
    xc = domain['xmax']/2
    yc = domain['ymax']/2
    ti = 0
    
    surf_seeds = []
    for theta in np.linspace(0, 2 * np.pi, n, endpoint=False):
        surf_seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                      r * np.sin(theta + 0.5 * ti) + yc, domain['zmax']])
                      
    surf_seeds_arr = np.array(surf_seeds)
    surf_seeds = tvtk.PolyData()
    surf_seeds.points = surf_seeds_arr
    return surf_seeds

def create_flux_surface(bfield,surf_seeds):
    #Make a streamline instance with the bfield
    surf_field_lines = tvtk.StreamTracer()
    surf_field_lines.input = bfield

    surf_field_lines.source = surf_seeds
    surf_field_lines.integrator = tvtk.RungeKutta4()
    surf_field_lines.maximum_propagation = 1000
    surf_field_lines.integration_direction = 'backward'
    surf_field_lines.update()

    #Create surface from 'parallel' lines
    surface = tvtk.RuledSurfaceFilter()
    surface.input = surf_field_lines.output
    surface.close_surface = True
    surface.pass_lines = True
    surface.offset = 0
    surface.distance_factor = 30
    surface.ruled_mode = 'point_walk'
#    surface.ruled_mode = 'resample'
#    surface.resolution = (10,1)
    surface.update()
    
    return surf_field_lines, surface

def update_flux_surface(surf_seeds, surf_field_lines, surface):
    #import pdb; pdb.set_trace()
    surf_field_lines.update()
    surface.update()   


def make_poly_norms(poly_data):
    poly_norms = tvtk.PolyDataNormals()
    poly_norms.input = poly_data
    poly_norms.compute_point_normals = True
    poly_norms.flip_normals = False
    poly_norms.update()
    
    return poly_norms

def norms_sanity_check(poly_norms):
    norm1 = poly_norms.output.point_data.normals[1000]
    norm_sanity = np.dot(norm1,
                         np.array(poly_norms.input.points.get_point(1000))-
                         np.array([63,63,poly_norms.input.points.get_point(1000)[2]]))         
    if norm_sanity < 0:
        poly_norms.flip_normals = not(poly_norms.flip_normals)
        poly_norms.update()
        
        passfail = False
    else:
        passfail = True
    return passfail, poly_norms

def get_surface_vectors(poly_norms, surf_bfield):
    """ Calculate the vector normal, vertically parallel and Azimuthally around
    the surface cont"""
    # Update the Normals
    poly_norms.update()

    passfail, poly_norms = norms_sanity_check(poly_norms)
#    print "pass norm check?", passfail
    
    normals = np.array(poly_norms.output.point_data.normals)
    
    
    parallels = surf_bfield / np.sqrt(np.sum(surf_bfield**2,axis=1))[:, np.newaxis]
    
    torsionals = np.cross(normals,parallels)
    torsionals /= np.sqrt(np.sum(torsionals**2,axis=1))[:, np.newaxis]
    
    return normals, torsionals, parallels

def interpolate_scalars(image_data, poly_data):
    """ Interpolate a imagedata scalars to a set points in polydata"""
    surface_probe_filter = tvtk.ProbeFilter(source=image_data,
                                          input=poly_data)
    surface_probe_filter.update()
    
    # Calculate Vperp, Vpar, Vphi
    surface_scalars = np.array(surface_probe_filter.output.point_data.scalars)
    return surface_probe_filter, surface_scalars
    
def interpolate_vectors(image_data, poly_data):
    """ Interpolate a imagedata vectors to a set points in polydata"""
    surface_probe_filter = tvtk.ProbeFilter(source=image_data,
                                          input=poly_data)
    surface_probe_filter.update()
    
    # Calculate Vperp, Vpar, Vphi
    surface_vectors = np.array(surface_probe_filter.output.point_data.vectors)
    return surface_probe_filter, surface_vectors

def update_interpolated_vectors(poly_data, surface_probe_filter):
    if poly_data:
        surface_probe_filter.input = poly_data
    surface_probe_filter.update()
    # Calculate Vperp, Vpar, Vphi
    surface_vectors = np.array(surface_probe_filter.output.point_data.vectors)
    
    return surface_vectors

def update_interpolated_scalars(poly_data, surface_probe_filter):
    if poly_data:
        surface_probe_filter.input = poly_data
    surface_probe_filter.update()
    # Calculate Vperp, Vpar, Vphi
    surface_scalars = np.array(surface_probe_filter.output.point_data.scalars)
    
    return surface_scalars
   
def get_surface_velocity_comp(surface_velocities, normals, torsionals, parallels):
    vperp = np.zeros(len(surface_velocities))
    vpar = np.zeros(len(surface_velocities))
    vphi = np.zeros(len(surface_velocities))
    
    for i in xrange(0,len(surface_velocities)):
        vperp[i] = np.dot(normals[i],surface_velocities[i])
        vpar[i] = np.dot(parallels[i],surface_velocities[i])
        vphi[i] = np.dot(torsionals[i],surface_velocities[i])
    
    return vperp, vpar, vphi

def get_the_line(bfield, surf_seeds, n):
    """Generate the vertical line on the surface"""
    the_line = tvtk.StreamTracer(input=bfield,
                                 source=tvtk.PolyData(
                                 points=np.array(
                                    [surf_seeds.points.get_point(n),[0,0,0]]
                                                )
                                                     )
                                 )

    the_line.integrator = tvtk.RungeKutta4()
    the_line.maximum_propagation = 1000
    the_line.integration_direction = 'backward'
    the_line.update()
    
    return the_line
    
def update_the_line(the_line, surf_seeds, seed, length):
    """ Updates the TD line at each time step, while making sure the length is fixed"""
    the_line.source.points = np.array([surf_seeds.get_point(seed), [0.0, 0.0, 0.0]])
    the_line.update()
    
    N = len(the_line.output.points)
    if N < length:
        print len(the_line.output.points)
        for i in list(the_line.output.points)[N-1:]:
            the_line.output.points.append(i)
    if N > length:
        print len(the_line.output.points)
        the_line.output.points = list(the_line.output.points)[:length]
    return the_line
    
def get_surface_indexes(surf_poly,the_line):
    point_locator = tvtk.PointLocator(data_set=surf_poly)
    
    surf_line_index = [] 
    surf_line_points = []
    for point in the_line.output.points:
        surf_line_index.append(point_locator.find_closest_point(point))
        surf_line_points.append(surf_poly.points.get_point(surf_line_index[-1]))
    
    return surf_line_index, surf_line_points

def write_step(file_name,surface,normals,parallels,torsionals,vperp,vpar,vphi):
    pd_par = tvtk.PointData(scalars=vpar,vectors=parallels)
    pd_par.scalars.name = "vpar"
    pd_par.vectors.name = "par"
    
    pd_perp = tvtk.PointData(scalars=vperp,vectors=normals)
    pd_perp.scalars.name = "vperp"
    pd_perp.vectors.name = "perp"
    
    pd_phi = tvtk.PointData(scalars=vphi,vectors=torsionals)
    pd_phi.scalars.name = "vphi"
    pd_phi.vectors.name = "phi"
    
    poly_out = surface.output
    poly_out.point_data.add_array(pd_par.scalars)
    poly_out.point_data.add_array(pd_par.vectors)
    
    poly_out.point_data.add_array(pd_perp.scalars)
    poly_out.point_data.add_array(pd_perp.vectors)
    
    poly_out.point_data.add_array(pd_phi.scalars)
    poly_out.point_data.add_array(pd_phi.vectors)
    
    w = tvtk.XMLPolyDataWriter(input=poly_out,file_name=file_name)
    w.write()

def write_flux(file_name, surface, surface_density, surface_va, surface_beta,
               surface_cs, Fpar, Fperp, Fphi):
    pd_density = tvtk.PointData(scalars=surface_density)
    pd_density.scalars.name = "surface_density"
    
    pd_va = tvtk.PointData(scalars=surface_va)
    pd_va.scalars.name = "surface_va"
    
    pd_beta = tvtk.PointData(scalars=surface_beta)
    pd_beta.scalars.name = "surface_beta"
    
    pd_cs = tvtk.PointData(scalars=surface_cs)
    pd_cs.scalars.name = "surface_cs"
    
    pd_Fpar = tvtk.PointData(scalars=Fpar)
    pd_Fpar.scalars.name = "Fpar"
    
    pd_Fperp = tvtk.PointData(scalars=Fperp)
    pd_Fperp.scalars.name = "Fperp"
    
    pd_Fphi = tvtk.PointData(scalars=Fphi)
    pd_Fphi.scalars.name = "Fphi"
    
    poly_out = surface
    poly_out.point_data.add_array(pd_density.scalars)
    poly_out.point_data.add_array(pd_va.scalars)
    poly_out.point_data.add_array(pd_beta.scalars)
    poly_out.point_data.add_array(pd_cs.scalars)
    poly_out.point_data.add_array(pd_Fpar.scalars)
    poly_out.point_data.add_array(pd_Fperp.scalars)
    poly_out.point_data.add_array(pd_Fphi.scalars)
    
    w = tvtk.XMLPolyDataWriter(input=poly_out,file_name=file_name)
    w.write()
    
def write_wave_flux(file_name, surface_poly, parallels, normals, torsionals,
                    Fwpar, Fwperp, Fwphi):

    pd_Fwpar = tvtk.PointData(scalars=Fwpar, vectors=parallels)
    pd_Fwpar.scalars.name = "Fwpar"
    pd_Fwpar.vectors.name = "par"
    
    pd_Fwperp = tvtk.PointData(scalars=Fwperp, vectors=normals)
    pd_Fwperp.scalars.name = "Fwperp"
    pd_Fwperp.vectors.name = "perp"
    
    pd_Fwphi = tvtk.PointData(scalars=Fwphi, vectors=torsionals)
    pd_Fwphi.scalars.name = "Fwphi"
    pd_Fwphi.vectors.name = "phi"
    
    poly_out = surface_poly
    poly_out.point_data.add_array(pd_Fwpar.scalars)
    poly_out.point_data.add_array(pd_Fwperp.scalars)
    poly_out.point_data.add_array(pd_Fwphi.scalars)
    
    w = tvtk.XMLPolyDataWriter(input=poly_out,file_name=file_name)
    w.write()

def read_step(filename):
    """ Read back in a saved surface file"""
    r = tvtk.XMLPolyDataReader(file_name=filename)
    r.update()
    return r.output

def get_data(poly_out,name):
    names = {}
    #Extract varibles from file
    for i in xrange(0,poly_out.point_data.number_of_arrays):
        names.update({poly_out.point_data.get_array_name(i):i})
    
    data = np.array(poly_out.point_data.get_array(names[name]))
    
    return data
    