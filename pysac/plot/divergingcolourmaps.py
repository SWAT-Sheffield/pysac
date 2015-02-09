# -*- coding: utf-8 -*-
"""
:Created on: Tue Jun 26 10:57:21 2012

:author: Stuart Mumford

Create Diverging colormaps like (Moreland 2005)
"""
from __future__ import division
import numpy as np
from colormath.color_objects import BaseRGBColor as RGBColor, LabColor

__all__ = ['get_mpl_colourmap','get_mayavi_colourmap']


""" Colour Conversion Routines """
def Lab_MSH(L,a,b):
    M = np.sqrt(L**2 + a**2 + b**2)
    s = (M > 0.001) * np.arccos(L/M)
    h = (s > 0.001) * np.arctan2(b,a)
    return M, s, h

def MSH_Lab(M,s,h):
    L = M * np.cos(s)
    a = M * np.sin(s) * np.cos(h)
    b = M * np.sin(s) * np.sin(h)
    return L,a,b

def RGB_Msh(R,G,B):
    colour = RGBColor(R*255,G*255,B*255,illuminant='d65')
    LAB = colour.convert_to('Lab')
    L,a,b = LAB.lab_l,LAB.lab_a,LAB.lab_b
    return Lab_MSH(L,a,b)

def Msh_RGB(M,s,h):
    L,a,b = MSH_Lab(M,s,h)
    colour = LabColor(L,a,b,illuminant='d65')
    RGB = colour.convert_to('RGB')
    return RGB.rgb_r/255.,RGB.rgb_g/255.,RGB.rgb_b/255.

def AdjustHue(M1,s1,h1,M2):
    if M1 >= M2:
        return h1
    else:
        hSpin = (s1 * np.sqrt(M2**2 - M1**2))/(M1*np.sin(s1))
        if h1 > -(np.pi/3.):
            return h1 + hSpin
        else:
            return h1 - hSpin


""" Msh colour space interpolation """
def InterpolateColour(rgb1,rgb2,interp,n,Mw,Sw):
    M1,s1,h1 = RGB_Msh(rgb1[0],rgb1[1],rgb1[2])
    M2,s2,h2 = RGB_Msh(rgb2[0],rgb2[1],rgb2[2])

    h1u = h1
    h2u = h2
    #if s1 < 0.05 and s2 > 0.05:
        #h2u = AdjustHue(M2,s2,h2,Mw)
    #elif s2 < 0.05 and s1 > 0.05:
        #h1u = AdjustHue(M1,s1,h1,Mw)

    Mout = np.r_[M1:Mw:(n*interp)*1j,Mw:M2:(n*(1 - interp))*1j]#,np.zeros([n/2])+M2]
    sout = np.r_[s1:Sw:(n*interp)*1j,Sw:s2:(n*(1 - interp))*1j]

    hout = np.r_[h1:h1u:(n*interp)*1j,h2u:h2:(n*(1 - interp))*1j]

    Msh_out = np.array([Mout,sout,hout])
    return Msh_out



""" Get colour map routines """
def get_mpl_colourmap(name,rgb1,rgb2,interp,n,Mw,Sw):
    Msh = InterpolateColour(rgb1,rgb2,interp,n,Mw,Sw)
    RGB_out = np.zeros((n,3))
    for i in xrange(len(Msh[0,:])):
        r,g,b = Msh_RGB(Msh[0,i],Msh[1,i],Msh[2,i])
        RGB_out[i,:] = np.array([r,g,b])
    import matplotlib as mpl
    linmap = mpl.colors.ListedColormap(RGB_out,name=name)
    return linmap

def get_mayavi_colourmap(rgb1,rgb2,interp,n,Mw,Sw,alpha):
    Msh = InterpolateColour(rgb1,rgb2,interp,n,Mw,Sw)
    RGB_out = np.zeros((256,4))
    for i in xrange(len(Msh[0,:])):
        r,g,b = Msh_RGB(Msh[0,i],Msh[1,i],Msh[2,i])
        RGB_out[i,:3] = np.array([r,g,b])*255
    RGB_out[:,-1] = alpha
    return RGB_out

def lsm_to_mayavi(lsm):
    """
    Convert a `matplotlib.colors.LinearSegmentedColormap` to a mayavi array.
    """
    r = np.array(lsm._segmentdata['red'])[:,0]
    g = np.array(lsm._segmentdata['green'])[:,0]
    b = np.array(lsm._segmentdata['blue'])[:,0]

    # rescale to 0-255
    r, g, b = r*255, g*255, b*255

    RGB = np.zeros([len(r),4], dtype=np.uint8)+255
    RGB[:,0] = r
    RGB[:,1] = g
    RGB[:,2] = b

    return RGB