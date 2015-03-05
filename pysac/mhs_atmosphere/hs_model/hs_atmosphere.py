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

def read_VAL3c_MTW(VAL_file=None, MTW_file=None, mu=0.602):
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
    MTW['p'].unit = u.Unit('dyne cm-2')
    MTW['rho'] = (MTW['p'] / k_B / MTW['T'] * m_p * mu).to('g cm-3')

    MTW['mu'] = mu

    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
#    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
    data.sort('Z')

    return data

#============================================================================
# interpolate the empirical data onto a Z array
#============================================================================
def interpolate_atmosphere(data, Z):
    """ This module generates a 1d array for the model plasma preesure, plasma
    density, temperature and mean molecular weight.
    """

    hdata = np.array(u.Quantity(data['Z']).to(u.m))
    # interpolate total pressure, temperature and density profiles
    pdata_f = UnivariateSpline(hdata,np.array(np.log(data['p'])),k=1, s=0.25)
    Tdata_f = UnivariateSpline(hdata,np.array(np.log(data['T'])),k=1, s=0.25)
    rdata_f = UnivariateSpline(hdata,np.array(np.log(data['rho'])),k=1, s=0.25)
    #s=0.0 to ensure all points are strictly used for ionisation state
    muofT_f = UnivariateSpline(hdata,np.array(np.log(data['mu'])),k=1, s=0.0)

    outdata = Table()
    outdata['Z'] = Z
    outdata['p'] = np.exp(pdata_f(Z.to(u.m))) * data['p'].unit
    outdata['T'] = np.exp(Tdata_f(Z.to(u.m))) * data['T'].unit
    outdata['rho'] = np.exp(rdata_f(Z.to(u.m))) * data['rho'].unit
    outdata['mu'] = np.exp(muofT_f(Z.to(u.m))) * u.one

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
        For HS equilibrium dp/dz = rho g., so cannot be isothermal?
    """
    p0 = model_pars['p0']
    r0 = 2.727e-07 * u.g/u.cm**3
    g0 = physical_constants['gravity']
    if option_pars['l_const']:
        pressure_Z = p0 * model_pars['chrom_scale']**3 /\
                            (model_pars['chrom_scale'] + Z)**3
        rho_Z = -p0 / g0 * 3. * model_pars['chrom_scale']**3/\
                            (model_pars['chrom_scale'] + Z)**4
        rtest = -p0 / g0 * 3. / model_pars['chrom_scale']
        model_pars['model'] += '_const'
    elif option_pars['l_sqrt']:
        pressure_Z = p0 *     model_pars['chrom_scale']**0.5/\
                             (model_pars['chrom_scale'] + Z)**0.5
        rho_Z = -0.5/g0 * p0 * model_pars['chrom_scale']**0.5/\
                             (model_pars['chrom_scale'] + Z)**1.5
        rtest = -0.5/g0 * p0 / model_pars['chrom_scale']
        model_pars['model'] += '_sqrt'
    elif option_pars['l_linear']:
        pressure_Z = p0 *     model_pars['chrom_scale']**1.5/\
                             (model_pars['chrom_scale'] + Z)**1.5
        rho_Z = -1.5/g0 * p0 * model_pars['chrom_scale']**1.5/\
                             (model_pars['chrom_scale'] + Z)**2.5
        rtest = -1.5/g0 * p0 / model_pars['chrom_scale']
        model_pars['model'] += '_linear'
    elif option_pars['l_square']:
        pressure_Z = p0 *     model_pars['chrom_scale']**3.5/\
                             (model_pars['chrom_scale'] + Z)**3.5
        rho_Z = -3.5/g0 * p0 * model_pars['chrom_scale']**3.5/\
                             (model_pars['chrom_scale'] + Z)**4.5
        rtest = -3.5/g0 * p0 / model_pars['chrom_scale']
        model_pars['model'] += '_square'
    else:
        raise ValueError("in hs_model.hs_atmosphere.get_spruit_hs set \
                  option_pars True for axial Alfven speed Z dependence")
    #to compare the derived density from hs-balance with VAL3c value:
    print'VAL rho(0) = ',r0.decompose(),' vs spruit rho(0) = ',rtest.decompose()
    Rgas_Z = u.Quantity(np.ones(Z.size), u.one)
    Rgas_Z *= physical_constants['boltzmann']/\
                physical_constants['proton_mass']/physical_constants['mu']

    return pressure_Z, rho_Z, Rgas_Z

#============================================================================
# Construct 3D hydrostatic profiles and include the magneto adjustments
#============================================================================
def vertical_profile(Z,
                     table,
                     magp0,
                     physical_constants, dz
                    ):
    """Return the vertical profiles for thermal pressure and density in 1D.
       Integrate in reverse from the corona to the photosphere to remove
       sensitivity to larger chromospheric gradients."""
    g0 = physical_constants['gravity'].to('m s-2')
    Rgas = u.Quantity(np.ones(table['Z'].size), u.one)
    Rgas *= (physical_constants['boltzmann']/\
                physical_constants['proton_mass']/table['mu']).to('m2 K-1 s-2')
    Rgas_Z  = Rgas[4:-4].copy()
    rdata   = u.Quantity(table['rho'], copy=True).to('kg m-3')
    rdata_Z = rdata[4:-4].copy()
    magp = magp0.to('kg m-1 s-2')
    # inverted SAC 4th order derivative scheme to minimise numerical error
    """evaluate upper boundary pressure from equation of state + enhancement,
       magp, which will be replaced by the mean magnetic pressure in the
       corona, then integrate from inner next pressure
    """
    table_T = u.Quantity(table['T'])
    linp_1 = table_T[-1]*rdata[-1]*Rgas[-1] + magp[-1]
    linp = u.Quantity(np.ones(len(Z)), unit=linp_1.unit)
    linp[-1] = table_T[-5]*rdata[-5]*Rgas[-5] + magp[-1]

#    for i in range(1,Z.size):
#        linp[-i-1] = (144.*linp[-i]+18.*linp[-i+1]
#                  -102.*(g0*rdata[-i-4]  )*dz
#                  - 84.*(g0*rdata[-i-5])*dz
#                  +  6.*(g0*rdata[-i-6])*dz
#                  )/162. + magp[-i-1] - magp[-i-0]
    for i in range(1,Z.size):
        linp[-i-1] = (1152.*linp[-i]
                  + 35.*(g0*rdata[-i-7])*dz
                  -112.*(g0*rdata[-i-6])*dz
                  -384.*(g0*rdata[-i-5])*dz
                  -784.*(g0*rdata[-i-4])*dz
                  + 77.*(g0*rdata[-i-3])*dz
                  )/1152. + magp[-i-1] - magp[-i-0]
    thermalp_Z = linp
    return thermalp_Z, rdata_Z, Rgas_Z
