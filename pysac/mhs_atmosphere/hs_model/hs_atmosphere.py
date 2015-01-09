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
    VAL3c['mu'] = 4.0/(3*0.74+1+VAL3c['n_e'].quantity/VAL3c['n_i'].quantity)

    MTW = Table.read(MTW_file, format='ascii', comment='#')
    MTW['Z'].unit = u.km
    MTW['T'].unit = u.K
    MTW['p'].unit = u.Unit('dyne/cm^2')
    MTW['rho'] = (MTW['p'] / k_B / MTW['T'] * m_p * mu).to('g cm-3')

    MTW['mu'] = mu

    data = astropy.table.vstack([VAL3c, MTW], join_type='inner')
    data.sort('Z')

    return data

def interpolate_atmosphere(data, Z):
    """ This module generates a 1d array for the model plasma preesure, plasma
    density, temperature and mean molecular weight.
    """

    hdata = data['Z'].quantity.value
    # interpolate total pressure, temperature and density profiles
    pdata_f = UnivariateSpline(hdata,np.array(np.log(data['p'])),k=1, s=0.25) * data['p'].unit
    Tdata_f = UnivariateSpline(hdata,np.array(np.log(data['T'])),k=1, s=0.25) * data['T'].unit
    rdata_f = UnivariateSpline(hdata,np.array(np.log(data['rho'])),k=1, s=0.25) * data['rho'].unit
    #s=0.0 to ensure all points are strictly used for ionisation state
    muofT_f = UnivariateSpline(hdata,np.array(np.log(data['mu'])),k=1, s=0.0) * u.one

    outdata = Table()
    outdata['Z'] = Z
    outdata['p'] = np.exp(pdata_f(Z))
    outdata['T'] = np.exp(Tdata_f(Z))
    outdata['rho'] = np.exp(rdata_f(Z))
    outdata['mu'] = np.exp(muofT_f(Z))


    return outdata


#----------------------------------------------------------------------------
# a simpler exponential atmosphere to test Spruit's analytical result
#----------------------------------------------------------------------------
def get_spruit_hs(
                   filenames,
                   Z,
                   scales,
                   model_pars,
                   physical_constants,
                   logical_pars,
                   plot
                 ):
    """ Reference data is collected to satisfy the global call for source data,
        but only the photospheric values of pressure and density are required
        for Spruit.
        Four options are available to select Alfven speed along the flux tube
        axis to be:
        constant, increase as the square root of Z, increase linearly and
        increase as the square 0f Z. These are approximate due to the effect on
        density of the non-zero magnetic tension force.
    """
    logical_pars['l_atmos_val3c_mtw'] = True
    logical_pars['l_spruit'] = False
    pdata_, Tdata_i, rdata_, muofT_i, [val ,mtw] = \
        interpolate_atmosphere(
                           filenames,
                           Z,
                           scales,
                           model_pars,
                           physical_constants,
                           logical_pars,
                           plot=True
                          )
    if logical_pars['l_const']:
        pdata_i = 1.1*pdata_[0]\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rdata_i = 0.5*rdata_[0]\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
    elif logical_pars['l_sqrt']:
        pdata_i = 1.1*pdata_[0]/(1+Z/model_pars['chrom_scale'])**0.25\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rdata_i = 0.5*rdata_[0]/(1+Z/model_pars['chrom_scale'])**0.25\
                     *np.exp(-4.*Z/model_pars['chrom_scale'])
    elif logical_pars['l_linear']:
        pdata_i = 1.1*pdata_[0]/(1+Z/model_pars['chrom_scale'])**1\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rdata_i = 0.5*rdata_[0]/(1+Z/model_pars['chrom_scale'])**1\
                     *np.exp(-4.*Z/model_pars['chrom_scale'])
    elif logical_pars['l_square']:
        pdata_i = 1.1*pdata_[0]/(1+Z/model_pars['chrom_scale'])**2\
                     *np.exp(-4.0*Z/model_pars['chrom_scale'])
        rdata_i = 0.5*rdata_[0]/(1+Z/model_pars['chrom_scale'])**2\
                     *np.exp(-4.*Z/model_pars['chrom_scale'])
    else:
        import sys
        sys.exit('in hs_model.hs_atmosphere.get_spruit_hs set \
                  logical_pars True for axial Alfven speed Z dependence')
    muofT_i[:] = physical_constants['mu']
    Tdata_i = pdata_i/rdata_i*physical_constants['mu']\
                             *physical_constants['proton_mass']\
                             /physical_constants['boltzmann']
    logical_pars['l_atmos_val3c_mtw'] = False
    logical_pars['l_spruit'] = True

    return pdata_i, Tdata_i, rdata_i, muofT_i, val, mtw

#============================================================================
# Construct 3D hydrostatic profiles and include the magneto adjustments
#============================================================================
def vertical_profile(Zint, Z,
                     pdata_i, rdata_i, Tdata_i, muofT_i,
                     magp,
                     physical_constants, dz,
                    ):
    """Return the vertical profiles for thermal pressure and density in 1D.
       Integrate in reverse from the corona to the photosphere to remove
       sensitivity to larger chromospheric gradients."""

    g0 = physical_constants['gravity']
    Rgas_v =np.zeros(Z.size)
    Rgas_v[:] = physical_constants['boltzmann']/\
                physical_constants['proton_mass']/muofT_i[4:-4]
    rdata_v = np.zeros(Z.size)
    rdata   = np.zeros(Zint.size)
    rdata[:] = rdata_i[:]
    #it may be necessary to enhance the density near the footpoint emprically
#    rdata[:] = rdata_i[:]*(1. + 1.5 * np.exp(-Zint[:]/0.00005)) #+ 1.*rdata_i.min()#+ 0.1*rdata_i[4]*np.exp(-Zint/Z[50])
    rdata_v[:] = rdata[4:-4]
    # inverted SAC 4th order derivative scheme to minimise numerical error
    linp=np.zeros(Z.size)
    """evaluate upper boundary pressure from equation of state + enhancement,
       magp, which will be replaced by the mean magnetic pressure in the
       corona, then integrate from inner next pressure
    """
    linp[-1] = Tdata_i[-5]*rdata_v[-1]*Rgas_v[-1] + magp[-1]
    linp[-2] = linp[-1] + magp[-2] - magp[-1] - ( 9.*g0*rdata[-1]*dz
                                                +19.*g0*rdata[-2]*dz
                                                - 5.*g0*rdata[-3]*dz
                                                +    g0*rdata[-4]*dz
                                            )/24.

    for i in range(2,Z.size):
        linp[-i-1] = (144.*linp[-i]+18.*linp[-i+1]
                  -102.*(g0*rdata[-i]  )*dz
                  - 84.*(g0*rdata[-i-1])*dz
                  +  6.*(g0*rdata[-i-2])*dz
                  )/162. + magp[-i-1] - magp[-i]

    thermalp_v = linp
    return thermalp_v, rdata_v, Rgas_v
