# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 17:45:48 2014

@author: sm1fg
"""

import numpy as np
import os
import astropy.units as u

hmi_model = {'photo_scale': 0.6*u.Mm,       #scale height for photosphere
             'chrom_scale': 0.1*u.Mm,      #scale height for chromosphere
             'corona_scale': 2.5e3*u.Mm,      #scale height for the corona
             'coratio': 0.03*u.one,  #footpoint portion scaling as corona 
             'model': 'hmi_model',
             'phratio': 0.15*u.one,  #footpoint portion scaling as photosphere
             'pixel': 0.36562475*u.Mm,      #(HMI pixel)
             'radial_scale': 0.10979002*u.Mm, #=> FWHM = half pixel
             'nftubes': 1,
             'B_corona': 0.*u.T,
             'pBplus': 1e-3*u.T}
hmi_model['chratio'] = 1*u.one - hmi_model['coratio'] - hmi_model['phratio']

mfe_setup = {'photo_scale': 0.60*u.Mm,
             'chrom_scale': 0.4*u.Mm,
             'corona_scale': 0.25*u.Mm,  #scale height for the corona
             'coratio': 0.0*u.one,
             'model': 'mfe_setup',
             'phratio': 0.0*u.one,
             'pixel': 0.36562475*u.Mm,  #(HMI pixel)
             'radial_scale': 0.03938*u.Mm,
             #'radial_scale': 0.14*u.Mm,
             'nftubes': 1,
             #'B_corona': 4.85e-4*u.T,
             'B_corona': 5.5e-4*u.T,
             'pBplus': 12.0e-4*u.T}
mfe_setup['chratio'] = 1*u.one - mfe_setup['coratio'] - mfe_setup['phratio']
#if 1D or 2D set unused dimensions to 0, and unrequired xyz limits to 1.
mfe_setup['Nxyz'] = [128,128,128] # 3D grid
mfe_setup['Nxyz'] = [129,129,128] # 3D grid
mfe_setup['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.6641221e-2*u.Mm,1.5877863*u.Mm] #grid size

spruit = {'photo_scale': 1.5*u.Mm,
          'chrom_scale': 0.5*u.Mm,
          'corona_scale': 100*u.Mm,      #scale height for the corona
          'coratio': 0.0*u.one,
          'model': 'spruit',
          'phratio': 0.0*u.one,
          'pixel': 0.1*u.Mm,              
          'radial_scale': 0.075*u.Mm,
          'nftubes': 1,
          'p0': 117200.0 * u.dyne/u.cm**2,
          'B_corona': 0.*u.T,
          'pBplus': 4.250e-4*u.T}
spruit['chratio'] = 1*u.one - spruit['coratio'] - spruit['phratio']
spruit['Nxyz'] = [128,128,256] # 3D grid
spruit['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.0*u.Mm,25.5*u.Mm] #grid size


paper1 = {'photo_scale': 0.6*u.Mm,
          'chrom_scale': 0.1*u.Mm,
          'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
          'coratio': 0.03*u.one,
          'model': 'paper1',
          'phratio': 0.0*u.one,
          'pixel': 0.36562475*u.Mm,              #(HMI pixel)
          'radial_scale': 0.10979002*u.Mm,
          'nftubes': 1,
          'B_corona': 9.2e-4*u.T,
          'pBplus': 1e-3*u.T}
paper1['chratio'] = 1*u.one - paper1['coratio'] - paper1['phratio']
paper1['Nxyz'] = [128,128,432] # 3D grid
paper1['xyz']  = [-1.27*u.Mm,1.27*u.Mm,-1.27*u.Mm,1.27*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size
# uncomment to produce comaparable data to mfe_setup
#paper1['Nxyz'] = [127,128,128] # 3D grid
#paper1['xyz']  = [-1*u.Mm,1*u.Mm,-1*u.Mm,1*u.Mm,3.5e-2*u.Mm,1.6*u.Mm] #grid size

paper2a = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2a',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 4,
           'B_corona': 8.4e-4*u.T,
           'pBplus': 1e-3*u.T}
paper2a['chratio'] = 1*u.one - paper2a['coratio'] - paper2a['phratio']
paper2a['Nxyz'] = [160,80,432] # 3D grid
paper2a['xyz']  = [-1.59*u.Mm,1.59*u.Mm,-0.79*u.Mm,0.79*u.Mm,0.*u.Mm,8.62*u.Mm] #grid size

paper2b = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 2.5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2b',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 4,
           'B_corona': 8.2e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2b['chratio'] = 1*u.one - paper2b['coratio'] - paper2b['phratio']
paper2b['Nxyz'] = [50,50,140] # 3D grid
paper2b['xyz']  = [-0.49*u.Mm,0.49*u.Mm,-0.49*u.Mm,0.49*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2c = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2c',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 15,
           'B_corona': 5.95e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2c['chratio'] = 1*u.one - paper2c['coratio'] - paper2c['phratio']
paper2c['Nxyz'] = [224,224,140] # 3D grid
paper2c['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-2.23*u.Mm,2.23*u.Mm,0*u.Mm,2.78*u.Mm] #grid size

paper2d = {'photo_scale': 0.6*u.Mm,
           'chrom_scale': 0.1*u.Mm,
           'corona_scale': 5e3*u.Mm,         #scale height for the corona
           'coratio': 0.03*u.one,
           'model': 'paper2d',
           'phratio': 0.0*u.one,
           'pixel': 0.36562475*u.Mm,              #(HMI pixel)
           'radial_scale': 0.10979002*u.Mm,
           'nftubes': 18,
           'B_corona': 5.95e-4*u.T,
           'pBplus': 1.0e-3*u.T}
paper2d['chratio'] = 1*u.one - paper2d['coratio'] - paper2d['phratio']
paper2d['Nxyz'] = [224,224,140] # 3D grid
paper2d['xyz']  = [-2.23*u.Mm,2.23*u.Mm,-0.79*u.Mm,0.79*u.Mm,0*u.Mm,2.78*u.Mm] #grid size


def get_coords(Nxyz, xyz):
    """
    get_coords returns a non-dimensional dictionary describing the domain
    coordinates.
    """
    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = u.Quantity(np.linspace(xyz[4].value, xyz[5].value, Nxyz[2]), unit=xyz.unit)
    Zext = u.Quantity(np.linspace(Z.min().value-4.*dz.value, Z.max().value+4.*dz.value, Nxyz[2]+8), unit=Z.unit)
    coords = {
              'dx':(xyz[1]-xyz[0])/(Nxyz[0]-1),
              'dy':(xyz[3]-xyz[2])/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xyz[0],
              'xmax':xyz[1],
              'ymin':xyz[2],
              'ymax':xyz[3],
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords
#-----------------------------------------------------------------------------
#
def get_hmi_map(
                model_pars, option_pars,
                indx, 
                dataset = 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits', 
                sunpydir = os.path.expanduser('~/sunpy/data/'),
                figsdir = os.path.expanduser('~/figs/hmi/'),
                l_newdata = False
               ):
    """ indx is 4 integers: lower and upper indices each of x,y coordinates 
#    dataset of the form 'hmi_m_45s_2014_07_06_00_00_45_tai_magnetogram_fits'
#    """
    from scipy.interpolate import RectBivariateSpline
    from sunpy.net import vso
    import sunpy.map
    client = vso.VSOClient()
    results = client.query(vso.attrs.Time("2014/07/05 23:59:50",
                                          "2014/07/05 23:59:55"), 
                           vso.attrs.Instrument('HMI'),
                           vso.attrs.Physobs('LOS_magnetic_field'))
    #print results.show()                       

    if l_newdata:
        if not os.path.exists(sunpydir):
            raise ValueError("in get_hmi_map set 'sunpy' dir for vso data\n"+ 
        "for large files you may want link to local drive rather than network")
        client.get(results).wait(progress=True)
    if not os.path.exists(figsdir):
        os.makedirs(figsdir)

    hmi_map = sunpy.map.Map(sunpydir+dataset)
    #hmi_map = hmi_map.rotate()
    #hmi_map.peek()
    s = hmi_map.data[indx[0]:indx[1],indx[2]:indx[3]] #units of Gauss Bz
    s *= u.G
    nx = s.shape[0]
    ny = s.shape[1]
    nx2, ny2 = 2*nx, 2*ny # size of interpolant 
    #pixel size in arc seconds
    dx, dy = hmi_map.scale.items()[0][1],hmi_map.scale.items()[1][1] 
    x, y = np.mgrid[
             hmi_map.xrange[0]+indx[0]*dx:hmi_map.xrange[0]+indx[1]*dx:1j*nx2,
             hmi_map.xrange[0]+indx[2]*dy:hmi_map.xrange[0]+indx[3]*dy:1j*ny2
                     ]
    #arrays to interpolate s from/to
    fx =   u.Quantity(np.linspace(x.min().value,x.max().value,nx), unit=x.unit)
    fy =   u.Quantity(np.linspace(y.min().value,y.max().value,ny), unit=y.unit)
    xnew = u.Quantity(np.linspace(x.min().value,x.max().value,nx2), unit=x.unit)
    ynew = u.Quantity(np.linspace(y.min().value,y.max().value,ny2), unit=y.unit)
    f  = RectBivariateSpline(fx,fy,s.to(u.T))
    #The initial model assumes a relatively small region, so a linear 
    #Cartesian map is applied here. Consideration may be required if larger
    #regions are of interest, where curvature or orientation near the lim
    #of the surface is significant. 
    s_int  = f(xnew,ynew) #interpolate s and convert units to Tesla
    s_int /= 4. # rescale s as extra pixels will sum over FWHM
    x_int  = x  * 7.25e5 * u.m    #convert units to metres
    y_int  = y  * 7.25e5 * u.m
    dx_int = dx * 7.25e5 * u.m
    dy_int = dy * 7.25e5 * u.m 
    FWHM  = 0.5*(dx_SI+dy_SI)
    smax  = max(abs(s.min()),abs(s.max())) # set symmetric plot scale
    cmin  = -smax*1e-4
    cmax  =  smax*1e-4
#    
#    filename = 'hmi_map'
#    import loop_plots as mhs
#    mhs.plot_hmi(
#             s*1e-4,x_SI.min(),x_SI.max(),y_SI.min(),y_SI.max(),
#             cmin,cmax,filename,savedir,annotate = '(a)'
#            )
#    filename = 'hmi_2x2_map'
#    mhs.plot_hmi(
#             s_SI*4,x_SI.min(),x_SI.max(),y_SI.min(),y_SI.max(),
#             cmin,cmax,filename,savedir,annotate = '(a)'
#            )
#
#    return s_SI, x_SI, y_SI, nx2, ny2, dx_SI, dy_SI, cmin, cmax, FWHM
    dz=(xyz[5]-xyz[4])/(Nxyz[2]-1)
    Z    = u.Quantity(np.linspace(xyz[4].value, xyz[5].value, Nxyz[2]), unit=xyz.unit)
    Zext = u.Quantity(np.linspace(Z.min().value-4.*dz.value, Z.max().value+4.*dz.value, Nxyz[2]+8), unit=Z.unit)
    coords = {
              'dx':(xyz[1]-xyz[0])/(Nxyz[0]-1),
              'dy':(xyz[3]-xyz[2])/(Nxyz[1]-1),
              'dz':(xyz[5]-xyz[4])/(Nxyz[2]-1),
              'xmin':xyz[0],
              'xmax':xyz[1],
              'ymin':xyz[2],
              'ymax':xyz[3],
              'zmin':xyz[4],
              'zmax':xyz[5],
              'Z':Z,
              'Zext':Zext
             }

    return coords

