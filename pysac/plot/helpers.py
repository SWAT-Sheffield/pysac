# -*- coding: utf-8 -*-
"""
:Created on: Mon Jul 23 14:00:50 2012

:author: Stuart Mumford

This is my custom colourmaps module
"""

import numpy as np
import matplotlib as mpl
from scipy import ndimage

from divergingcolourmaps import *

__all__ = ['FixedCentre',
           'jetblack', 'jetblack_mayavi', 'red_temp', 'red_temp_mayavi',
           'fieldlines']

class FixedCentre(mpl.colors.Normalize):
    """
    Normalise with a Fixed Centre
    """
    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale(result)
        vmin, vmax = self.vmin, self.vmax
        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0.5)   # Or should it be all masked?  Or 0.5?
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            # ma division is very slow; we can take a shortcut
            resdat = result.data
            resdat[(result.data == 0.).nonzero()] = 0.0
            resdat[(result.data < 0.).nonzero()] /= (-vmin * 2.)
            resdat[(result.data > 0.).nonzero()] /= (vmax * 2.)
            result = np.ma.array(resdat+0.5, mask=result.mask, copy=False)
        if is_scalar:
            result = result[0]
        return result

jbdict = {'red': [(0.0,0.0,0.16),
                 (0.25,0.0,0.0),
                 (0.45,0.0,0.0),
                 (0.5,0.0,0.0),
                 (0.55,0.5,0.5),
                 (0.75,1.0,1.0),
                 (1.0,0.80,0.0)],
         'green': [(0.0,0.0,1.0),
                   (0.25,0.5,0.5),
                   (0.45,0.0,0.0),
                   (0.5,0.0,0.0),
                   (0.55,0.0,0.0),
                   (0.75,0.58,0.58),
                   (1.0,1.0,0.0)],
         'blue': [(0.0,0.0,0.80),
                  (0.25,1.0,1.0),
                  (0.45,0.5,0.5),
                  (0.5,0.0,0.0),
                  (0.55,0.0,0.0),
                  (0.75,0.0,0.0),
                  (1.0,0.16,0.0)]}

jetblack = mpl.colors.LinearSegmentedColormap('jetblack',jbdict,N=2048)

red_temp_dict = {
                    'red': ((0.0,0.0,0.0),
                            (0.6901,1.0,1.0),
                            (1.0,1.0,1.0)),
                    'green': ((0.0,0.0,0.0),
                              (0.4705,0.0,0.0),
                              (1.0,1.0,1.0)),
                    'blue': ((0.0,0.0,0.0),
                             (0.7450,0.0,0.0),
                             (1.0,1.0,1.0))
                }
                
red_temp = mpl.colors.LinearSegmentedColormap('redtemp',
                                                     red_temp_dict,
                                                     N=2048)

jetblack_mayavi = []
for i in np.linspace(0,1.0,255):
    jetblack_mayavi.append(np.array(jetblack(i))*255)
jetblack_mayavi = np.array(jetblack_mayavi)

#Define a colormap
red_temp_mayavi = np.zeros([256,4])
red_temp_mayavi[0:176,0] = np.linspace(0,255,176)
red_temp_mayavi[176:,0] = 255
red_temp_mayavi[121:,1] = np.linspace(0,255,135)
red_temp_mayavi[191:,2] = np.linspace(0,255,65)
red_temp_mayavi[:,3] = 255



def fieldlines(v1, v2, seeds, dS=1):
    """
    Simple 2D Euler fieldline integrator.
    
    Parameters
    ----------
    v1,v2 : numpy.ndarray (2D)
        y and x vector arrays
    seeds : numpy.ndarray
        array of coordinate (array indexes) pairs
    dS : float
        step size (default value is 1.0)
    """
    
    field = []
    #Get extent of domain
    max1 = v1.shape[0]
    max2 = v2.shape[1]
    min2 = 0
    min1 = 0
    #integrate for each fieldline
    for seed in seeds:
        c1,c2 = seed
        out1 = [c1] #first point is seed
        out2 = [c2]
        cnt = 0
        
        while (c1 <= max1 and c1 >= min1) and (c2 <= max2 and c2 >= min2):
            #Interpolate the vector field to the coord
            coords = np.array([[c1],[c2]])
            v1i = ndimage.map_coordinates(v1,coords)[0]
            v2i = ndimage.map_coordinates(v2,coords)[0]
            vi = np.sqrt(v1i**2 + v2i**2)
            
            d1 = ( v1i * dS )/ vi
            d2 = ( v2i * dS )/ vi
            c1 -= d1 
            c2 -= d2
            
            out1.append(c1)
            out2.append(c2)
            
            cnt += 1
            if cnt > 500: # Maximum iteration limit
                print "limit"
                break
        out = np.zeros([len(out1),len(out2)])
        out[0] = out1
        out[1] = out2
        field.append(out)
    
    return np.array(field)