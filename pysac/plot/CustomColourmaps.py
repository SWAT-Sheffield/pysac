# -*- coding: utf-8 -*-
"""
:Created on: Mon Jul 23 14:00:50 2012

:author: Stuart Mumford

This is my custom colourmaps module
"""

import numpy as np
import matplotlib as mpl
from divergingcolourmaps import *

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