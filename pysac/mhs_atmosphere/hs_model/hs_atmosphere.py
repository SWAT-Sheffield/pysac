# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 11:37:39 2014

@author: sm1fg

    Generate a 1D non-magnetic atmosphere vector based on an empirical model
    based on observational data, or specify an analytical hydrostatic
    equilibrium atmosphere.

"""
import numpy as np
from scipy.interpolate import UnivariateSpline

import astropy.table
from astropy.table import Table
import astropy.units as u
from astropy.constants import k_B, m_p

__all__ = ['read_VAL3c_MTW', 'interpolate_atmosphere', 'get_spruit_hs', 'vertical_profile']

#============================================================================
# Read in and interpolate HD atmosphere
#============================================================================

def read_VAL3c_MTW(VAL_file=None, MTW_file=None, mu=0.6):
    """
    Read in the data from Table 12 in Vernazza (1981) and combine with
    McWhirter (1975).

    Parameters
    ----------
    VAL_file : string
        The data file for the VAL3c atmosphere, defaults to
        `pysac.mhs_atmosphere.hs_model.VALIIIc_data`

    MTW_file : string
        The data file for the McWhirter atmosphere, defaults to
        `pysac.mhs_atmosphere.hs_model.MTWcorona_data`

    mu : float
        The mean molecular weight ratio for the corona. defaults to 0.6.

    Returns
    -------
    data : `astropy.table.Table`
        The combined data, sorted by Z.
    """
    from . import VALIIIc_data, MTWcorona_data
    if not VAL_file:
        VAL_file = VALIIIc_data
    if not MTW_file:
        MTW_file = MTWcorona_data

    VAL3c = Table.read(VAL_file, format='ascii', comment='#')
    VAL3c['Z'].unit = u.km
    VAL3c['rho'].unit = u.Unit('g cm-3')
    VAL3c['p'].unit = u.Unit('dyne/cm^2')
    VAL3c['T'].unit = u.K
    VAL3c['n_i'].unit = u.one/u.cm**3
    VAL3c['n_e'].unit = u.one/u.cm**3

    # Calculate the mean molecular weight ratio
    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e']/VAL3c['n_i'])
#    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e'].quantity/VAL3c['n_i'].quantity)

    MTW = Table.read(MTW_file, format='ascii', comment='#')
    MTW['Z'].unit = u.km
    MTW['T'].unit = u.K
    MTW['p'].unit = u.Unit('dyne/cm^2')
    MTW['rho'] = (MTW['p'] / k_B / MTW['T'] * m_p * mu).to('g cm-3')

    MTW['mu'] = mu

    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
#    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
    data.sort('Z')

    return data

def interpolate_atmosphere(data, Z):
    """ This module generates a 1d array for the model plasma preesure, plasma
    density, temperature and mean molecular weight.
    """

    hdata = np.array(data['Z'])
    # interpolate total pressure, temperature and density profiles
    pdata_f = UnivariateSpline(hdata,np.array(np.log(data['p'])),k=1, s=0.25)
    Tdata_f = UnivariateSpline(hdata,np.array(np.log(data['T'])),k=1, s=0.25)
    rdata_f = UnivariateSpline(hdata,np.array(np.log(data['rho'])),k=1, s=0.25)
    #s=0.0 to ensure all points are strictly used for ionisation state
    muofT_f = UnivariateSpline(hdata,np.array(np.log(data['mu'])),k=1, s=0.0)

    outdata = Table()
    outdata['Z'] = Z
    outdata['p'] = np.exp(pdata_f(Z)) * data['p'].unit
    outdata['T'] = np.exp(Tdata_f(Z)) * data['T'].unit
    outdata['rho'] = np.exp(rdata_f(Z)) * data['rho'].unit
    outdata['mu'] = np.exp(muofT_f(Z)) * u.one


    return outdata


#----------------------------------------------------------------------------
# a simpler exponential atmosphere to test Spruit's analytical result
#----------------------------------------------------------------------------
def get_spruit_hs(
                   Z,
                   model_pars,
                   physical_constants,
                   option_pars
                 ):
    """ photospheric values of pressure and density are taken from VAL3c.
        Four options are available to select Alfven speed along the flux tube
        axis to be:
        constant, increase as the square root of Z, increase linearly and
        increase as the square 0f Z. We apply Bz~exp(-2z/chrom_scale) hence
        for Alfven speed \sqrt(B^2/rho) constant rho~exp(-4z/chrom_scale)...
        These are approximate due to the effect on density of the non-zero
        magnetic tension force.
        For HS equilibrium dp/dz = rho g.
    FAG: the hs balance needs to be rechecked, especially for last three
         options. But code compiles correctly
    """
    p0 = 117200.0 * u.dyne/u.cm**2
    r0 = 2.727e-07 * u.g/u.cm**3
    g0 = physical_constants['gravity']
    if option_pars['l_const']:
        pressure_Z = p0\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rho_Z = -4*p0/(g0*model_pars['chrom_scale'])\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
    elif option_pars['l_sqrt']:
        pressure_Z = p0/(1+Z/model_pars['chrom_scale'])**0.5\
                     *np.exp(-4.0*Z**0.5/model_pars['chrom_scale'])
        rho_Z = -  p0 (2 * g0 * model_pars['chrom_scale']**2 * (1 +
                          Z/model_pars['chrom_scale'])**1.5
                  ) * np.exp(-4.*Z/model_pars['chrom_scale'])
    elif option_pars['l_linear']:
        pressure_Z = p0/(1+Z/model_pars['chrom_scale'])**1\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rho_Z = -4*p0/g0/model_pars['chrom_scale']/(
                  1+Z/model_pars['chrom_scale'])**1\
                * np.exp(-4. * Z / model_pars['chrom_scale'])
    elif option_pars['l_square']:
        pressure_Z = p0/(1+Z/model_pars['chrom_scale'])**2\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rho_Z = -4*p0/g0/model_pars['chrom_scale']/(1+Z/model_pars['chrom_scale'])**2\
                     *np.exp(-4.*Z/model_pars['chrom_scale'])
    else:
        raise ValueError("in hs_model.hs_atmosphere.get_spruit_hs set \
                  option_pars True for axial Alfven speed Z dependence")
    rtest = -4*p0/g0/model_pars['chrom_scale']
    #to compare the derived density from hs-balance with VAL3c value:
    print'VAL rho(0) = ',r0.decompose(),' vs spruit rho(0) = ',rtest.decompose()
    Rgas_Z = u.Quantity(np.ones(Z.size), u.one)
    Rgas_Z *= physical_constants['boltzmann']/\
                physical_constants['proton_mass']/physical_constants['mu']

    return pressure_Z, rho_Z, Rgas_Z

#============================================================================
# Construct 3D hydrostatic profiles and include the magneto adjustments
#============================================================================
def vertical_profile(Zint, Z,
                     table,
                     magp,
                     physical_constants, dz,
                    ):
    """Return the vertical profiles for thermal pressure and density in 1D.
       Integrate in reverse from the corona to the photosphere to remove
       sensitivity to larger chromospheric gradients."""

    g0 = physical_constants['gravity']
    Rgas_v = u.Quantity(np.ones(Z.size), u.one)
    Rgas_v *= physical_constants['boltzmann']/\
                physical_constants['proton_mass']/table['mu'][4:-4]
    rdata   = u.Quantity(table['rho'], copy=True)
    #it may be necessary to enhance the density near the footpoint emprically
#    rdata[:] = rdata_i[:]*(1. + 1.5 * np.exp(-Zint[:]/0.00005)) #+ 1.*rdata_i.min()#+ 0.1*rdata_i[4]*np.exp(-Zint/Z[50])
    rdata_v = rdata[4:-4].copy()
    # inverted SAC 4th order derivative scheme to minimise numerical error
    """evaluate upper boundary pressure from equation of state + enhancement,
       magp, which will be replaced by the mean magnetic pressure in the
       corona, then integrate from inner next pressure
    """
    table_T = u.Quantity(table['T'])
    linp_1 = table_T[-5]*rdata_v[-1]*Rgas_v[-1] + magp[-1]
    linp = u.Quantity(np.ones(len(Z)), unit=linp_1.unit)
    linp[-1] = linp_1
    linp[-2] = linp[-1] + magp[-2] - magp[-1] - ( 9.*g0*rdata[-1]*dz+19.*g0*rdata[-2]*dz- 5.*g0*rdata[-3]*dz+g0*rdata[-4]*dz)/24.

    for i in range(2,Z.size):
        linp[-i-1] = (144.*linp[-i]+18.*linp[-i+1]
                  -102.*(g0*rdata[-i]  )*dz
                  - 84.*(g0*rdata[-i-1])*dz
                  +  6.*(g0*rdata[-i-2])*dz
                  )/162. + magp[-i-1] - magp[-i]

    thermalp_v = linp
    return thermalp_v, rdata_v, Rgas_v
