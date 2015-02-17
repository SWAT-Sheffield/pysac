# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:13:23 2015

@author: sm1fg
"""

import numpy as np


def make_1d_slices(ds, var_label, oneD_arrays):
    """extract 1D arrays from gdf data for plotting
    """
    var_ = ds.index.grids[0][var_label]

#    oneD_arrays.create[var_label]
    oneD_arrays.update({var_label:{}})
    oneD_arrays[var_label]['axis'] = var_[ds.domain_dimensions[0]/2,
                                          ds.domain_dimensions[1]/2,:]
    oneD_arrays[var_label]['edge'] = var_[0,        0        ,:]
    oneD_arrays[var_label]['mean'] = np.mean(np.mean(var_, axis=0), axis=0)
    oneD_arrays[var_label]['min' ] = np.min(np.min(var_, axis=0), axis=0)
    oneD_arrays[var_label]['max' ] = np.max(np.max(var_, axis=0), axis=0)
    oneD_arrays[var_label]['Z'   ] = np.linspace(ds.domain_left_edge[2],
                                                 ds.domain_right_edge[2],
                                                 ds.domain_dimensions[2])

    return oneD_arrays

