from __future__ import absolute_import

import importlib
import random
import math

import numpy as np

import distmesh as dm

def distmesh2d(shape, boundary_box, fixed_points, l=0.1):
    if shape == "square":
        return square(boundary_box, fixed_points, l)

def square(boundary_box, fixed_points, l=0.1):
    """Square, with size function point and line sources"""
    dist_function = lambda p: dm.drectangle(p, boundary_box[0],boundary_box[2],boundary_box[1],boundary_box[3])
    return dm.distmesh2d(dist_function, dm.huniform, l, boundary_box, fixed_points)
    #fd = lambda p: np.sqrt((p**2).sum(1))-1.0
    #return dm.distmesh2d(fd, dm.huniform, 0.06, (-1,-1,1,1))

