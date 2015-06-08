# encoding: utf-8
from __future__ import absolute_import
import importlib

import numpy as np

from shapegen import pointgen3d 
from mesh.meshgen import meshgen3d
import utils

def equally_spaced_longitudes3d(): 
    # l - initial length of triangle sides. Change it to vary traingle size
    boundary_box = [0,0,0,20,20,20]
    l = 3.5
    mesh = meshgen3d(boundary_box, l, include_corners=False)
    curve1 = pointgen3d.sphere_arc(mesh.bbox, 0, 10)
    curve2 = pointgen3d.sphere_arc(mesh.bbox, 2*np.pi/3, 10)
    curve3 = pointgen3d.sphere_arc(mesh.bbox, 4*np.pi/3, 10)
    shapes = [curve1, curve2, curve3]
    points  = np.array(shapes)
    lambdas = [0.001]
    mus = [0.00001]
    is_closed = False
    return mesh, mesh.triangles, mesh.edges, points, lambdas, mus, is_closed

def differently_spaced_longitudes3d(): 
    # l - initial length of triangle sides. Change it to vary traingle size
    boundary_box = [0,0,0,20,20,20]
    l = 4
    mesh = meshgen3d(boundary_box, l, include_corners=False)
    curve1 = pointgen3d.sphere_arc(mesh.bbox, 0, 10)
    curve2 = pointgen3d.sphere_arc(mesh.bbox, np.pi/4, 10)
    curve3 = pointgen3d.sphere_arc(mesh.bbox, 9*np.pi/8, 10)
    shapes = [curve1, curve2, curve3]
    points  = np.array(shapes)
    lambdas = [0.001]
    mus = [0.00001]
    is_closed = False
    return mesh, mesh.triangles, mesh.edges, points, lambdas, mus, is_closed

def fine_curves_on_sphere(): 
    mesh = utils.load_mesh3d(dirname="curves_on_sphere")
    input_currents = utils.load_input_currents(len(mesh.edges), 3, dirname='curves_on_sphere')
    t, q, r = utils.load_solutions(len(mesh.triangles), len(mesh.edges), 3, dirname='curves_on_sphere')
    return mesh, input_currents, t, q, r

