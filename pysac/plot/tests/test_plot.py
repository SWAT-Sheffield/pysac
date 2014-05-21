# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 14:46:23 2014

@author: stuart
"""
import numpy as np
import matplotlib as mpl

import pysac.plot

def test_FixedCentre():
    fc = pysac.plot.FixedCentre([0,1])

def test_jetblack():
    assert isinstance(pysac.plot.jetblack, mpl.colors.LinearSegmentedColormap)

def test_red_temp():
    assert isinstance(pysac.plot.red_temp, mpl.colors.LinearSegmentedColormap)

def test_jetblack_mayavi():
    assert isinstance(pysac.plot.jetblack_mayavi, np.ndarray)
    assert pysac.plot.jetblack_mayavi.shape == (255,4)

def test_red_temp_mayavi():
    assert isinstance(pysac.plot.red_temp_mayavi, np.ndarray)
    assert pysac.plot.red_temp_mayavi.shape == (256,4)