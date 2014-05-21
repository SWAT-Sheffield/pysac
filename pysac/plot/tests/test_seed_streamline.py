# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 14:23:07 2014

@author: stuart
"""

import numpy as np

from pysac.plot.mayavi_seed_streamlines import SeedStreamline


def test_creation():
    seeds = np.random.random([10,3])
    sS = SeedStreamline(seed_points = seeds)

def test_change():
    seeds = np.random.random([10,3])
    sS = SeedStreamline(seed_points = seeds)
    seeds = np.random.random([10,3])
    sS.seed_points = seeds
    sS.update_pipeline()