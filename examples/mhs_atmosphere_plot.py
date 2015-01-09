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

slc = yt.SlicePlot(ds, fields='density_bg', normal='x')
slc.save('~/yt.png')