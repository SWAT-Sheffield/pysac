##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import operator

import numpy

from vtk import vtkDataArray, vtkImageData, vtkImageExport, vtkImageImport, vtkLookupTable
import vtk.util.numpy_support
from vtk.util import vtkConstants

def array_to_vtk_image(array, copy_data, data_type="scalar"):
    """ Create an ``vtkImage`` matching the contents and type of given array. If
        ``copy_data`` is ``True``, then the data of the will be copied. Otherwise the
        data will be shared, and the array **must not** be destroyed before the 
        ``vtkImage``. ``data_type`` specifies how the array should be interpreted : 
        either as a n-array of scalars (``data_type="scalar"``) or as an n-1 array
        of vectors (``data_type="vector"``).
    """
    
    if data_type not in ["scalar", "vector"] :
        raise Exception("Unknown data_type: {0}".format(repr(data_type)))
    
    if data_type == "scalar" :
        ndim = array.ndim
    elif data_type == "vector" :
        ndim = array.ndim-1
    
    if ndim > 3 :
        raise Exception(
            "Cannot convert a {0} array of dimension {1}".format(data_type, 
                                                                 array.ndim))
    
    importer = vtkImageImport()
    
    if numpy.iscomplexobj(array) :
        # Get the first element of the array
        element = array.flat.next()
        scalar_type = vtk.util.numpy_support.get_vtk_array_type(element.real.dtype)
    else :
        scalar_type = vtk.util.numpy_support.get_vtk_array_type(array.dtype)
    importer.SetDataScalarType(scalar_type)
    
    if data_type == "scalar" :
        number_of_components = 1
    elif data_type == "vector" :
        number_of_components = array.shape[ndim]
    if numpy.iscomplexobj(array) :
        number_of_components *= 2
    importer.SetNumberOfScalarComponents(number_of_components)
    
    extent = 6*[0]
    extent[1] = array.shape[ndim-1]-1
    if ndim >= 2 :
        extent[3] = array.shape[ndim-2]-1
    if ndim >= 3 :
        extent[5] = array.shape[ndim-3]-1
    importer.SetDataExtent(extent)
    importer.SetWholeExtent(extent)
    
    size = array.itemsize*reduce(operator.mul, array.shape, 1)
    if copy_data :
        importer.CopyImportVoidPointer(array, size)
    else :
        importer.SetImportVoidPointer(array, size)
    
    importer.Update()

    return importer.GetOutput()

def vtk_image_to_array(vtk_image) :
    """ Create an ``numpy.ndarray`` matching the contents and type of given image. 
        If the number of scalars components in the image is greater than 1, then
        the ndarray will be 4D, otherwise it will be 3D. 
    """
    
    exporter = vtkImageExport()
    exporter.SetInput(vtk_image)
    
    # Create the destination array
    extent = vtk_image.GetWholeExtent()
    shape = [extent[5]-extent[4]+1,
             extent[3]-extent[2]+1,
             extent[1]-extent[0]+1]
    if vtk_image.GetNumberOfScalarComponents() > 1:
        shape += [vtk_image.GetNumberOfScalarComponents()]
    dtype = vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType())
    array = numpy.zeros(shape, dtype=dtype)
    
    exporter.Export(array)
    
    return array