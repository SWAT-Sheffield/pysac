# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:36:52 2013

@author: stuart
"""

import pysac.io as sacio

infile = sacio.VACdata("/home/stuart/iceberg_data/VAC/data/3D_data/3D_tube_128_128_128.ini")
print infile

outfile = sacio.VACdata("testout.out",mode='w')


header = {
        'filehead':'3D_mhd33',
        'ndim': 3,
        'it': 0,
        't': 0,
        'nw': 13,
        'nx': [128,128,128],
        'varnames': ['z', 'x', 'y', 'h', 'm1', 'm2', 'm3', 'e', 'b1', 'b2', 'b3',
                     'eb', 'rhob', 'bg1', 'bg2', 'bg3', 'gamma', 'eta', 
                     'grav1', 'grav2', 'grav3'],
        'neqpar': 7,
        'eqpar': [1.66666667, 0.0, 1.0, -274.0, 0.0, 0.0, 0.0]
        }
outfile.header = header
#outfile.x = x
#outfile.w = w
#outfile.write_step()
#outfile.close()




#newinfile = sacio.VACdata("testout.out")
#print newinfile.header
#
#outfile2 = sacio.VACdata("testout.h5",mode='w')
#outfile2.header = infile.header
#outfile2.x = infile.x
#outfile2.w = infile.w
#outfile2.write_step()
#outfile2.close()
#
#newinfile2 = sacio.VACdata("testout.h5")
