# encoding: utf-8
'''
3D experiment cases
-------------------
'''
from __future__ import absolute_import
import importlib
import os

import numpy as np

from medianshape.simplicial import pointgen3d
from medianshape.simplicial.meshgen import meshgen3d, get_mesh_surface
import medianshape.experiment.inout as inout
    
def equally_spaced_longitudes3ds(): 
    '''
    Hi
    '''
    '''
    # l - initial length of triangle sides. Change it to vary traingle size
    boundary_box = [0,0,0,20,20,20]
    l = 0.6
    mesh = meshgen3d(boundary_box, l, include_corners=False)
    smesh = get_mesh_surface(mesh)
    inout.save_data(mesh, dirname=os.path.abspath("data/mesh"))
    inout.save_data(smesh, dirname=os.path.abspath("data/smesh"))
    print "Saved"
    exit()
    '''
    mesh  = inout.load_mesh3d(dirname='data/mesh')
    smesh = inout.load_mesh3d(dirname='data/smesh')
    '''
    curve1 = pointgen3d.sphere_arc(mesh.bbox, 0, 10)
    curve2 = pointgen3d.sphere_arc(mesh.bbox, 2*np.pi/3, 10)
    curve3 = pointgen3d.sphere_arc(mesh.bbox, 4*np.pi/3, 10)
    '''
    curve2 = pointgen3d.sphere_arc(mesh.bbox, np.pi/18, 10)
    curve1 = pointgen3d.sphere_arc(mesh.bbox, np.pi + np.pi/6, 10)

    shapes = [curve1, curve2]
    input_points  = np.array(shapes)
    lambdas = [0.001]
    mus = [0.00001]
    #lambdas = [0.0000001]
    #mus = [0.000000001]
    is_closed = False
    alpha2 = np.array([0])
    alpha2 = np.append(alpha2, np.linspace(0.4999, 0.5, 3))
    alpha2 = np.append(alpha2, np.linspace(0.500001, 0.50001, 50))
    alpha2 = np.append(alpha2, np.array([1]))
    alpha2 = alpha2.reshape(alpha2.size, 1) 
    alpha1 = (1-alpha2).reshape(alpha2.size, 1)
    alphas = np.hstack((alpha1, alpha2))
    return mesh, smesh, smesh.triangles, smesh.edges, input_points, lambdas, mus, alphas, is_closed
