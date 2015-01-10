# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 13:55:17 2014

@author: sm1fg

This is the main module to construct a magnetohydrostatic solar atmosphere,
given a specified magnetic network of self-similar magnetic flux tubes and
save the output to gdf format.

To select an existing configuration set the model label, Nxyz, xyz_SI and
any other special parameters, then execute mhs_atmopshere.

To add new configurations:
add the model label options to get_logical in parameters/logical.py;
add options to get_model in parameters/model_pars.py;
add alternative empirical data sets to hs_model/;
add option to interpolate_atmosphere in hs_model/hs_atmosphere.py;
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
from pysac.mhs_atmosphere.parameters.model_pars import mfe_setup as model_pars
model = model_pars['model']

local_procs=1
#standard set of logical switches
logical_pars = atm.get_logical(model, l_mpi, l_SI=True, l_gdf=True)

#standard conversion to dimensionless units and physical constants
scales, physical_constants = \
    atm.get_parameters(model, l_mpi, logical_pars, size)

#if 1D or 2D set unused dimensions to 0, and unrequired xyz limits to 1.
Nxyz = [128,128,128] # 3D grid
xyz = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,35*u.km,1.6*u.Mm] # xyz limits SI/CGS units

#obtain code coordinates and model parameters in code units
coords = atm.get_coords(Nxyz, u.Quantity(xyz))

#from pysac.mhs_atmosphere.hs_model import VALIIIc_data as VAL
#from pysac.mhs_atmosphere.hs_model import MTWcorona_data as MTW

#filenames = [VAL, MTW]
# uncomment and switch to l_const/l_sqrt/l_linear/l_square as required
logical_pars['l_const'] = True
#interpolate the hs 1D profiles from empirical data source[s]
empirical_data = atm.read_VAL3c_MTW(mu=physical_constants['mu'])

table = \
    atm.interpolate_atmosphere(empirical_data,
                               coords['Zext']
                              )
pressure_1d, temperature_1d, rho_1d, muofT_1d = \
    table['p'], table['T'], table['rho'], table['mu']
#==============================================================================
#calculate 1d hydrostatic balance from empirical density profile
#==============================================================================
# the hs pressure balance is enhanced by pressure equivalent to the residual mean
# coronal magnetic pressure contribution once the magnetic field is applied
magp_meanz = np.ones(len(coords['Z'])) * u.one
magp_meanz *= model_pars['pBplus']**2/(2*physical_constants['mu0'])

pressure_Z, rho_Z, Rgas_Z = atm.vertical_profile(
                                                 coords['Zext'],
                                                 coords['Z'],
                                                 table,
                                                 magp_meanz,
                                                 physical_constants,
                                                 coords['dz'],
                                                 )

#==============================================================================
# load flux tube footpoint parameters
#==============================================================================
# axial location and value of Bz at each footpoint
xi, yi, Si = atm.get_flux_tubes(
                                model,
                                model_pars,
                                coords,
                                logical_pars
                               )
#==============================================================================
# split domain into processes if mpi
#==============================================================================
ax, ay, az = np.mgrid[coords['xmin']:coords['xmax']:1j*Nxyz[0],
                      coords['ymin']:coords['ymax']:1j*Nxyz[1],
                      coords['zmin']:coords['zmax']:1j*Nxyz[2]]

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
                                             model_pars, logical_pars,
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
                                             logical_pars,
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
# select the 1D array spanning the local mpi process
indz = np.where(coords['Z'] >= z.min()) and np.where(coords['Z'] <= z.max())
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
datadir = os.path.expanduser('~/mhs_atmosphere/'+model+'/')
filename = datadir + model + logical_pars['suffix']
if not os.path.exists(datadir):
    os.makedirs(datadir)
sourcefile = datadir + model + '_sources' + logical_pars['suffix']
auxfile = datadir + model + '_aux' + logical_pars['suffix']

# save the variables for the initialisation of a SAC simulation
atm.save_SACvariables(model,
              filename,
              rho,
              Bx,
              By,
              Bz,
              energy,
              logical_pars,
              physical_constants,
              coords,
              Nxyz
             )
# save the balancing forces as the background source terms for SAC simulation
atm.save_SACsources(model,
              sourcefile,
              Fx,
              Fy,
              logical_pars,
              physical_constants,
              coords,
              Nxyz
             )
# save auxilliary variable and 1D profiles for plotting and analysis
Rgas = np.zeros(x.shape)
Rgas[:] = Rgas_z
temperature = pressure/rho/Rgas
if not logical_pars['l_hdonly']:
    inan = np.where(magp <=1e-7*pressure.min())
    magpbeta = magp
    magpbeta[inan] = 1e-7*pressure.min()  # low pressure floor to avoid NaN
    pbeta  = pressure/magpbeta
else:
    pbeta  = magp+1.0    #dummy to avoid NaN
alfven = np.sqrt(2.*physical_constants['mu0']*magp/rho)
cspeed = np.sqrt(physical_constants['gamma']*pressure/rho)
atm.save_auxilliary1D(
              model,
              auxfile,
              pressure_m,
              rho_m,
              temperature,
              pbeta,
              alfven,
              cspeed,
              Btensx,
              Btensy,
              pressure_Z,
              rho_Z,
              Rgas_Z,
              logical_pars,
              physical_constants,
              coords,
              Nxyz
             )

