# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import math

import numpy as np

import distmesh as dm

def distmesh2d(shape, bbox, fixed_points, l=0.1):
    if shape == "square":
        return square_mesh(bbox, fixed_points, l)

def distmesh3d(shape, bbox, fixed_points, l=0.1):
    print shape, bbox, fixed_points
    if shape == "cuboid":
        print "hi"
        return cuboid_mesh(bbox, fixed_points, l)
    if shape == "ball":
        return ball(bbox, None, l=l)

def square_mesh(bbox, fixed_points, l=0.1):
    """Square, with size function point and line sources"""
    dist_function = lambda p: dm.drectangle(p, bbox[0], bbox[2], \
                                                bbox[1],bbox[3])
    points, triangles = dm.distmesh2d(dist_function, dm.huniform, l, bbox, fixed_points)
    return points, triangles 

def cuboid_mesh(bbox, fixed_points, l=0.1):
    '''
        Generates 3-simplex, tetrahedral mesh in 3D using DistMesh
    '''
    dist_function = lambda p: dcubiod(p, bbox[0], bbox[1], 
                                        bbox[2], bbox[3],
                                        bbox[4], bbox[5]) 
    points, tetrahedras = dm.distmeshnd(dist_function, dm.huniform, l, bbox) 
    return points, tetrahedras

def dcubiod(p, x1, y1, z1, x2, y2, z2):
    """Signed distance function for cubiod with corners (x1,y1,z1), (x1,y1,z2),
    (x1,y2,z1), (x1,y2,z2), (x2,y1,z1), (x2,y2,z1), (x2,y1,z2), (x2,y2,z2).

    This has an incorrect distance to the four corners. 
    """
    return -np.minimum(np.minimum(np.minimum(np.minimum(np.minimum\
            (-y1+p[:,1],y2-p[:,1]),-x1+p[:,0]),x2-p[:,0]),-z1+p[:,2]), z2-p[:,2])

def ball(bbox, fixed_points, l):                                                                                           
    fd = lambda p: np.sqrt((p**2).sum(1))-float(bbox[5])
    return dm.distmeshnd(fd, dm.huniform, l, bbox, fixed_points) 
