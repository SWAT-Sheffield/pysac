# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 12:52:31 2015

@author: fred
"""
import os
import numpy as np
import pysac.yt as sacyt
import pysac.mhs_atmosphere as atm
from pysac.mhs_atmosphere.parameters.model_pars import paper2d as model_pars
import astropy.units as u

l_mpi=False
paper = model_pars['model']
datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+paper+'/')
figsdir = os.path.expanduser('~/Documents/mhs_atmosphere/figs/'+paper+'/')
if not os.path.exists(figsdir):
    os.makedirs(figsdir)
file_ = datadir+paper+'.gdf'
ds = sacyt.SACGDFDataset(file_)
figname = figsdir+paper+'_3dplot.jpg'
fkeys = ['thermal_pressure','plasma_beta',
         'mag_field_x','mag_field_y','mag_field_z']

option_pars = atm.set_options(model_pars, l_mpi, l_gdf=True)
coords = atm.get_coords(model_pars['Nxyz'], u.Quantity(model_pars['xyz']))
xc, yc, S = atm.get_flux_tubes(model_pars, coords, option_pars)
if 'paper1' in file_:
    figxy = [1000,2000]
    view = (-55.0, 90.0, 21.0, np.array([0,0,4.25]))   
    nrad, nangle = 5, 4
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = [0.35]#ds.domain_right_edge[2].in_units('Mm').value * 0.125
if 'paper2a' in file_:
    figxy = [1100,2000]
    view = (-60.0, 90.0, 20.0, np.array([0,0,4.25]))
    nrad, nangle = 5, 4
    maxr = ds.domain_right_edge[0].in_units('Mm').value
    maxz = [0.5]#ds.domain_right_edge[2].in_units('Mm').value * 0.125
if 'paper2b' in file_:
    figxy = [1200,2000]
    view = (-45., 90., 7., np.array([0,0,1.35]))
    nrad, nangle = 5, 3
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = [0.275]#ds.domain_right_edge[2].in_units('Mm').value * 0.125
if 'paper2c' in file_:
    figxy = [2000,2000]
    view = (-45., 90., 12.5, np.array([0,0,1.35]))
    nrad, nangle = 2, 3
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    maxz = [0.15,ds.domain_right_edge[2].in_units('Mm').value*0.5]
if 'paper2d' in file_:
    figxy = [2000,2000]
    view = (-45., 90., 9.5, np.array([0,0,1.35]))
    maxr = ds.domain_right_edge[0].in_units('Mm').value# * 3.
    maxz = [0.15,ds.domain_right_edge[2].in_units('Mm').value*0.5]
    nrad, nangle = 2 , 3
if 'mfe_setup' in file_:
    figxy = [3200,2000]
    view = (-45., 90., 6.5, np.array([0,0,0.95]))
    maxr = ds.domain_right_edge[1].in_units('Mm').value
    nrad, nangle = 5 , 5
    maxz = 1.05#ds.domain_right_edge[2].in_units('Mm').value * 0.125
seeds = [[xc[0][0].value,yc[0][0].value,maxz[0]]]
for ix in range(0,xc.size):
    for ti,r in enumerate(maxr*(1.01-np.exp(-(np.linspace(0, 1, nrad))**2))):
        for theta in np.linspace(0, 2 * np.pi, nangle, endpoint=False):
            for iz in range(0,len(maxz)):
                seeds.append([r * np.cos(theta + r + 0.5 * ti +
                          np.random.uniform()) + xc[ix].value, 
                          r * np.sin(theta + r + 0.5 * ti + 
                          np.random.uniform()) + yc[ix].value, maxz[iz]])
seeds = np.array(seeds)
# ensure all seeds are located inside the computational domain
nix = np.where(seeds.T[0] >   maxr * .99)
if nix[0].size > 0:
    seeds.T[0][nix] =   maxr * .99
nix = np.where(seeds.T[0] < - maxr * .99)
if nix[0].size > 0:
    seeds.T[0][nix] = - maxr * .99
nix = np.where(seeds.T[1] >   maxr * .99)
if nix[0].size > 0:
    seeds.T[1][nix] =   maxr * .99
nix = np.where(seeds.T[1] < - maxr * .99)
if nix[0].size > 0:
    seeds.T[1][nix] = - maxr * .99
atm.make_3d_plot(ds, figname, fields=fkeys, figxy=figxy, view=view, 
                 seeds=seeds)
