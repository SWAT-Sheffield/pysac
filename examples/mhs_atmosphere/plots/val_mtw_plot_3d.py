# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: fred
"""
import os
import numpy as np
import pysac.yt as sacyt
import pysac.mhs_atmosphere as atm
from pysac.mhs_atmosphere.parameters.model_pars import paper1 as model_pars
import astropy.units as u
#papers = ['paper1','paper2a','paper2b','paper2c','paper2d','mfe_setup'

l_mpi=False
paper = model_pars['model']
datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+paper+'/')
figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+paper+'/')
if not os.path.exists(figsdir):
    os.makedirs(figsdir)
#open all gdf files in the model directory
file_ = datadir+paper+'.gdf'
ds = sacyt.SACGDFDataset(file_)
figname = figsdir+paper+'_3dplot.eps'
fkeys = ['density','alfven_speed',
         'mag_field_x','mag_field_y','mag_field_z']

option_pars = atm.set_options(model_pars, l_mpi, l_gdf=True)
coords = atm.get_coords(model_pars['Nxyz'], u.Quantity(model_pars['xyz']))
xc, yc, S = atm.get_flux_tubes(model_pars, coords, option_pars)
maxr = ds.domain_right_edge[1].in_units('Mm').value
maxz = ds.domain_right_edge[2].in_units('Mm').value
nrad, nangle = 7, 4
seeds = [[0.,0.,maxz]]
for ix in xc:
    for ti,r in enumerate(maxr*(1.01-np.exp(-(np.linspace(0, 2, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            seeds.append([r * np.cos(theta + r + 0.5 * ti) + xc[0].value,
                          r * np.sin(theta + r + 0.5 * ti) + yc[0].value, maxz])
    seeds = np.array(seeds)

if 'paper1' in file_:
    figxy = [500,550]
    view = (-45.0, 90.0, 22.0, np.array([0,0,4.25]))   
if 'paper2a' in file_:
    figxy = [500,750]
    view = (-45.0, 90.0, 20.0, np.array([0,0,3.75]))
    xc = 0
    yc = 0
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = ds.domain_right_edge[2].in_units('Mm').value
    nrad, nangle = 7 , 6
    seeds = [[xc,yc,maxr]]
    for ti,r in enumerate(maxr*(1-np.exp(-(np.linspace(0, 1, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                          r * np.sin(theta + 0.5 * ti) + yc, maxz])
    seeds = np.array(seeds)
if 'paper2b' in file_:
    figxy = [500,350]
    view = (-45., 90., 20., np.array([0,0,3.75]))
    xc = 0
    yc = 0
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = ds.domain_right_edge[2].in_units('Mm').value
    nrad, nangle = 7 , 6
    seeds = [[xc,yc,maxr]]
    for ti,r in enumerate(maxr*(1-np.exp(-(np.linspace(0, 1, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                          r * np.sin(theta + 0.5 * ti) + yc, maxz])
    seeds = np.array(seeds)
if 'paper2c' in file_:
    figxy = [500,750]
    view = (-45., 90., 20., np.array([0,0,3.75]))
    xc = 0
    yc = 0
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = ds.domain_right_edge[2].in_units('Mm').value
    nrad, nangle = 7 , 6
    seeds = [[xc,yc,maxr]]
    for ti,r in enumerate(maxr*(1-np.exp(-(np.linspace(0, 1, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                          r * np.sin(theta + 0.5 * ti) + yc, maxz])
    seeds = np.array(seeds)
if 'paper2d' in file_:
    figxy = [500,750]
    view = (-45., 90., 20., np.array([0,0,3.75]))
    xc = 0
    yc = 0
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = ds.domain_right_edge[2].in_units('Mm').value
    nrad, nangle = 7 , 6
    seeds = [[xc,yc,maxr]]
    for ti,r in enumerate(maxr*(1-np.exp(-(np.linspace(0, 1, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                          r * np.sin(theta + 0.5 * ti) + yc, maxz])
    seeds = np.array(seeds)
if 'mfe_setup' in file_:
    figxy = [500,750]
    view = (-45., 90., 20., np.array([0,0,3.75]))
    xc = 0
    yc = 0
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = ds.domain_right_edge[2].in_units('Mm').value
    nrad, nangle = 7 , 6
    seeds = [[xc,yc,maxr]]
    for ti,r in enumerate(maxr*(1-np.exp(-(np.linspace(0, 1, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            seeds.append([r * np.cos(theta + 0.5 * ti) + xc,
                          r * np.sin(theta + 0.5 * ti) + yc, maxz])
    seeds = np.array(seeds)
atm.make_3d_plot(ds, figname, fields=fkeys, figxy=figxy, view=view, 
                 seeds=seeds)

