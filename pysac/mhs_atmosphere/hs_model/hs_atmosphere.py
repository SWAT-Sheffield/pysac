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

#============================================================================
# Read in and interpolate HD atmosphere
#============================================================================
def interpolate_atmosphere(filenames, 
                           Z, 
                           scales, 
                           model_pars, 
                           physicalconstants, 
                           logical_pars, 
                           plot=False):

    """ This module generates a 1d array for the model plasma preesure, plasma
    density, temperature and mean molecular weight. 
    """
#----------------------------------------------------------------------------
    if logical_pars['l_atmos_val3c_mtw']:
        """
    "filenames": a list of strings containing the paths of model data 
    for this module we combine data from Table 12 in Vernazza 1981 with
    McWhirter 1975 and convert units to SI
    In Vernazza data is:
        H(km) | \rho (g cm^-3) | Ptot (dyn cm^-2) | T (K) | Pgas / Ptot
    and in McWhirter 1975 data is:  
        H(km) | T (K)          | Ptot (dyn cm^-2) |
    The data is interpolated onto an array in Z applicable to the model
    resolution
    """
        VAL3c = np.loadtxt(filenames[0])[::-1,:]
        VAL3c[:,0] *= 1e3 #km -> m
        VAL3c[:,1] *= 1e3 #g/cm^3 -> kg/m^3
        VAL3c[:,2] /= 1e1 #dyne/cm^2 -> Pascals(N/m^2)
        muTV=4.0/(3*0.74+1+VAL3c[:,6]/VAL3c[:,5])
        nzv=VAL3c[:,0].size
        MTW =np.loadtxt(filenames[1])[::-1,:]
        MTW[:,0] *= 1e3 #km -> m
        MTW[:,1] *= 1   #Kelvin
        MTW[:,2] /= 1e1 #dyne/cm^2 -> Pascals(N/m^2)
        nzm=MTW[:,0].size
        #Combine both sets into single array
        hdata=np.zeros(nzv+nzm)
        hdata[0:nzv]=VAL3c[:,0]
        hdata[nzv:nzv+nzm]=MTW[:,0]
        hdata /= scales['length']
    
        pdata=np.zeros(nzv+nzm)
        pdata[0:nzv]=VAL3c[:,2]
        pdata[nzv:nzv+nzm]=MTW[:,2]
        pdata /= scales['energy density']
    
        Tdata=np.zeros(nzv+nzm)
        Tdata[0:nzv]=VAL3c[:,3]
        Tdata[nzv:nzv+nzm]=MTW[:,1] 
        Tdata /= scales['temperature']
    
        rdata=np.zeros(nzv+nzm)
        rdata[0:nzv]=VAL3c[:,1]
        rdata /= scales['density']   
        #MTW[:,2], kB, mp and mu are in code units so no rescale needed
        rdata[nzv:nzv+nzm] = MTW[:,2]/physicalconstants['boltzmann']/MTW[:,1] \
                        *physicalconstants['proton_mass']*physicalconstants['mu']
    
        muofT = np.zeros(nzv+nzm) # mean molecular weight
        muofT[0:nzv] = muTV
        muofT[nzv:nzv+nzm] = physicalconstants['mu']
        # interpolate total pressure, temperature and density profiles
        pdata_f = UnivariateSpline(hdata,np.log(pdata),k=1, s=0.25)
        Tdata_f = UnivariateSpline(hdata,np.log(Tdata),k=1, s=0.25)
        rdata_f = UnivariateSpline(hdata,np.log(rdata),k=1, s=0.25)
        #s=0.0 to ensure all points are strictly used for ionisation state
        muofT_f = UnivariateSpline(hdata,np.log(muofT),k=1, s=0.0)
        pdata_i = np.exp(pdata_f(Z))
        Tdata_i = np.exp(Tdata_f(Z))
        rdata_i = np.exp(rdata_f(Z))
        muofT_i = np.exp(muofT_f(Z))

        source_data = [VAL3c,MTW]                          
        if plot:
            return pdata_i, Tdata_i, rdata_i, muofT_i, source_data
        else:        
            return pdata_i, Tdata_i, rdata_i, muofT_i
#----------------------------------------------------------------------------
# a simpler exponential atmosphere to test Spruit's analytical result
#----------------------------------------------------------------------------
    if logical_pars['l_spruit']:
        pdata_i, Tdata_i, rdata_i, muofT_i, VAL3c, MTW\
               = get_spruit_hs(
                               filenames, 
                               Z, 
                               scales, 
                               model_pars,
                               physicalconstants,
                               logical_pars,
                               plot
                              )
        source_data = [VAL3c,MTW]                          
        if plot:
            return pdata_i, Tdata_i, rdata_i, muofT_i, source_data
        else:        
            return pdata_i, Tdata_i, rdata_i, muofT_i



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
                     scales
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
