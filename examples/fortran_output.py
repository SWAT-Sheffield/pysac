# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:36:52 2013

@author: stuart
"""

import pysac.io as sacio

infile = sacio.VACdata("/home/stuart/iceberg_data/VAC/data/mhdmodes_hpc_hd.out")
print infile

outfile = sacio.VACdata("testout.out",mode='w')
outfile.header = infile.header
outfile.x = infile.x
outfile.w = infile.w
outfile.write_step()
outfile.close()

newinfile = sacio.VACdata("testout.out")
print newinfile.header