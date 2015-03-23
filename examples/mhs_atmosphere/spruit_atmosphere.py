# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 13:55:17 2014

@author: sm1fg

This is the main module to construct a magnetohydrostatic solar atmosphere,
given a specified magnetic network of self-similar magnetic flux tubes and
save the output to gdf format.

To select an existing configuration change the import as model_pars, set Nxyz,
xyz_SI and any other special parameters, then execute mhs_atmopshere.

To add new configurations:
add the model options to set_options in parameters/options.py;
add options required in parameters/model_pars.py;
add alternative empirical data sets to hs_model/;
add alternativ table than interploate_atmosphere in hs_model/hs_atmosphere.py;
add option to get_flux_tubes in mhs_model/flux_tubes.py

If an alternative formulation of the flux tube is required add options to
construct_magnetic_field and construct_pairwise_field in
mhs_model/flux_tubes.py

Plotting options are included in plot/mhs_plot.py
"""

import os
import numpy as np
import pysac.mhs_atmosphere as atm
import astropy.units as u
from pysac.mhs_atmosphere.parameters.model_pars import spruit as model_pars
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
#set up model parameters
#==============================================================================
local_procs=1
#standard set of logical switches
option_pars = atm.set_options(model_pars, l_mpi, l_gdf=True)
#standard conversion to dimensionless units and physical constants
scales, physical_constants = \
    atm.get_parameters()
# select the option in the next line
option_pars['l_const'] = True
# Alfven speed constant along the axis of the flux tube
if option_pars['l_const']:
    option_pars['l_B0_quadz'] = True
    model_pars['chrom_scale'] *= 5e1
    model_pars['p0'] *= 1.5e1
    physical_constants['gravity'] *= 1.
    model_pars['radial_scale'] *= 1.
# Alfven speed proportional to sqrt(Z) along the axis of the flux tube
elif option_pars['l_sqrt']:
    option_pars['l_B0_rootz'] = True
    model_pars['chrom_scale'] *= 5.65e-3
    model_pars['p0'] *= 1.
    physical_constants['gravity'] *= 7.5e3
    model_pars['radial_scale'] *= 0.7
# Alfven speed proportional to Z along the axis of the flux tube
elif option_pars['l_linear']:
    option_pars['l_B0_rootz'] = True
    model_pars['chrom_scale'] *= 0.062
    model_pars['p0'] *= 3e2
    physical_constants['gravity'] *= 8e3
    model_pars['radial_scale'] *= 1.
# Alfven speed proportional to Z^2 along the axis of the flux tube
elif option_pars['l_square']:
    option_pars['l_B0_rootz'] = True
    model_pars['chrom_scale'] *= 1.65
    model_pars['p0'] *= 2e4
    physical_constants['gravity'] *= 5e4
    model_pars['radial_scale'] *= 1.
# Alfven speed not defined along the axis of the flux tube
else:
    option_pars['l_B0_rootz'] = True
    model_pars['chrom_scale'] *= 1.
    model_pars['p0'] *= 1.


#obtain code coordinates and model parameters in astropy units
coords = atm.get_coords(model_pars['Nxyz'], u.Quantity(model_pars['xyz']))

#==============================================================================
#calculate 1d hydrostatic balance from empirical density profile
#==============================================================================
pressure_Z, rho_Z, Rgas_Z = atm.get_spruit_hs(coords['Z'],
                                                  model_pars,
                                                  physical_constants,
                                                  option_pars
                                                  )
#==============================================================================
# load flux tube footpoint parameters
#==============================================================================
# axial location and value of Bz at each footpoint
xi, yi, Si = atm.get_flux_tubes(
                                model_pars,
                                coords,
                                option_pars
                               )
#==============================================================================
# split domain into processes if mpi
#==============================================================================
ax, ay, az = np.mgrid[coords['xmin']:coords['xmax']:1j*model_pars['Nxyz'][0],
                      coords['ymin']:coords['ymax']:1j*model_pars['Nxyz'][1],
                      coords['zmin']:coords['zmax']:1j*model_pars['Nxyz'][2]]

# split the grid between processes for mpi
if l_mpi:
    x_chunks = np.array_split(ax, size, axis=0)
    y_chunks = np.array_split(ay, size, axis=0)
    z_chunks = np.array_split(az, size, axis=0)

    x = comm.scatter(x_chunks, root=0)
    y = comm.scatter(y_chunks, root=0)
    z = comm.scatter(z_chunks, root=0)
else:
    x, y, z = ax, ay, az

x = u.Quantity(x, unit=coords['xmin'].unit)
y = u.Quantity(y, unit=coords['ymin'].unit)
z = u.Quantity(z, unit=coords['zmin'].unit)
#==============================================================================
# initialize zero arrays in which to add magnetic field and mhs adjustments
#==============================================================================
Bx   = u.Quantity(np.zeros(x.shape), unit=u.T)  # magnetic x-component
By   = u.Quantity(np.zeros(x.shape), unit=u.T)  # magnetic y-component
Bz   = u.Quantity(np.zeros(x.shape), unit=u.T)  # magnetic z-component
pressure_m = u.Quantity(np.zeros(x.shape), unit=u.Pa) # magneto-hydrostatic adjustment to pressure
rho_m = u.Quantity(np.zeros(x.shape), unit=u.kg/u.m**3)      # magneto-hydrostatic adjustment to density
# initialize zero arrays in which to add balancing forces and magnetic tension
Fx   = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)  # balancing force x-component
Fy   = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)  # balancing force y-component
# total tension force for comparison with residual balancing force
Btensx  = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)
Btensy  = u.Quantity(np.zeros(x.shape), unit=u.N/u.m**3)
#==============================================================================
#calculate the magnetic field and pressure/density balancing expressions
#==============================================================================
for i in range(0,model_pars['nftubes']):
    for j in range(i,model_pars['nftubes']):
        if rank == 0:
            print'calculating ij-pair:',i,j
        if i == j:
            pressure_mi, rho_mi, Bxi, Byi ,Bzi, B2x, B2y =\
                atm.construct_magnetic_field(
                                             x, y, z,
                                             xi[i], yi[i], Si[i],
                                             model_pars, option_pars,
                                             physical_constants,
                                             scales
                                            )
            Bx, By, Bz = Bxi+Bx, Byi+By ,Bzi+Bz
            Btensx += B2x
            Btensy += B2y
            pressure_m += pressure_mi
            rho_m += rho_mi
        else:
            pressure_mi, rho_mi, Fxi, Fyi, B2x, B2y =\
                atm.construct_pairwise_field(
                                             x, y, z,
                                             xi[i], yi[i],
                                             xi[j], yi[j], Si[i], Si[j],
                                             model_pars,
                                             option_pars,
                                             physical_constants,
                                             scales
                                            )
            pressure_m += pressure_mi
            rho_m += rho_mi
            Fx   += Fxi
            Fy   += Fyi
            Btensx += B2x
            Btensy += B2y

# clear some memory
del pressure_mi, rho_mi, Bxi, Byi ,Bzi, B2x, B2y
#==============================================================================
# Construct 3D hs arrays and then add the mhs adjustments to obtain atmosphere
#==============================================================================
# select the 1D array spanning the local mpi process; the add/sub of dz to
# ensure all indices are used, but only once
indz = np.where(coords['Z'] >= z.min()-0.1*coords['dz']) and \
       np.where(coords['Z'] <= z.max()+0.1*coords['dz'])
pressure_z, rho_z, Rgas_z = pressure_Z[indz], rho_Z[indz], Rgas_Z[indz]
# local proc 3D mhs arrays
pressure, rho = atm.mhs_3D_profile(z,
                                   pressure_z,
                                   rho_z,
                                   pressure_m,
                                   rho_m
                                  )
magp = (Bx**2 + By**2 + Bz**2)/(2.*physical_constants['mu0'])
if rank ==0:
    print'max B corona = ',magp[:,:,-1].max().decompose()
energy = atm.get_internal_energy(pressure,
                                                  magp,
                                                  physical_constants)
#============================================================================
# Save data for SAC and plotting
#============================================================================
# set up data directory and file names
# may be worthwhile locating on /data if files are large
datadir = os.path.expanduser('~/Documents/mhs_atmosphere/'+model_pars['model']+'/')
filename = datadir + model_pars['model'] + option_pars['suffix']
if not os.path.exists(datadir):
    os.makedirs(datadir)
sourcefile = datadir + model_pars['model'] + '_sources' + option_pars['suffix']
aux3D = datadir + model_pars['model'] + '_3Daux' + option_pars['suffix']
aux1D = datadir + model_pars['model'] + '_1Daux' + option_pars['suffix']

# save the variables for the initialisation of a SAC simulation
atm.save_SACvariables(
              filename,
              rho,
              Bx,
              By,
              Bz,
              energy,
              option_pars,
              physical_constants,
              coords,
              model_pars['Nxyz']
             )
# save the balancing forces as the background source terms for SAC simulation
atm.save_SACsources(
              sourcefile,
              Fx,
              Fy,
              option_pars,
              physical_constants,
              coords,
              model_pars['Nxyz']
             )
# save auxilliary variable and 1D profiles for plotting and analysis
Rgas = u.Quantity(np.zeros(x.shape), unit=Rgas_z.unit)
Rgas[:] = Rgas_z
temperature = pressure/rho/Rgas
if not option_pars['l_hdonly']:
    inan = np.where(magp <=1e-7*pressure.min())
    magpbeta = magp
    magpbeta[inan] = 1e-7*pressure.min()  # low pressure floor to avoid NaN
    pbeta  = pressure/magpbeta
else:
    pbeta  = magp+1.0    #dummy to avoid NaN
alfven = np.sqrt(2.*magp/rho)
if rank == 0:
    print'Alfven speed Z.min to Z.max =',\
    alfven[model_pars['Nxyz'][0]/2,model_pars['Nxyz'][1]/2, 0].decompose(),\
    alfven[model_pars['Nxyz'][0]/2,model_pars['Nxyz'][1]/2,-1].decompose()
cspeed = np.sqrt(physical_constants['gamma']*pressure/rho)
atm.save_auxilliary3D(
              aux3D,
              pressure_m,
              rho_m,
              temperature,
              pbeta,
              alfven,
              cspeed,
              Btensx,
              Btensy,
              option_pars,
              physical_constants,
              coords,
              model_pars['Nxyz']
             )
atm.save_auxilliary1D(
              aux1D,
              pressure_Z,
              rho_Z,
              Rgas_Z,
              option_pars,
              physical_constants,
              coords,
              model_pars['Nxyz']
             )

