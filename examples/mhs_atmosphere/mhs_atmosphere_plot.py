# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: stuart
"""
import os
import glob

import yt

model = 'spruit'
datadir = os.path.expanduser('~/mhs_atmosphere/'+model+'/')

files = glob.glob(datadir+'/*')
files.sort()

print(files)


ds = yt.load(files[0])

# uncomment for axis swapping for normal='y'
ds.coordinates.x_axis = {0: 2, 1: 0, 2: 1, 'x': 2, 'y': 0, 'z': 1}
ds.coordinates.y_axis = {0: 1, 1: 2, 2: 0, 'x': 1, 'y': 2, 'z': 0}

slc = yt.SlicePlot(ds, normal='y', fields='density_bg')
slc.save('~/yt.png')