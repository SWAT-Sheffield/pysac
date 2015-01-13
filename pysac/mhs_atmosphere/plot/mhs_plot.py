# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 14:55:17 2015

@author: sm1fg
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import commands
#import pysac.io as sacio
import yt
import h5py
import astropy.units as u
import pysac.io.gdf_writer as gdf
import pysac.mhs_atmosphere as atm
#==============================================================================
#check whether mpi is required and the number of procs = size
#==============================================================================
try:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    l_mpi = True
    l_mpi = l_mpi and (size != 1)
except ImportError:
    l_mpi = False
    rank = 0
    size = 1
#==============================================================================
# Read in saved data for plotting
#==============================================================================
# Set path and file names to read
model = 'spruit'
homedir = os.environ['HOME']
#standard set of logical switches
option_pars = atm.set_options(model, l_mpi, l_gdf=True)
datadir = os.path.expanduser(homedir+'/mhs_atmosphere/'+model+'/')
filename = datadir + model + option_pars['suffix']
if not os.path.exists(datadir+model):
    os.makedirs(datadir+model)
sourcefile = datadir + model + 'sources' + option_pars['suffix']
auxfile = datadir + model + 'aux' + option_pars['suffix']
test = yt.load(filename)
