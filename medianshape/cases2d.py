# encoding: utf-8
from __future__ import absolute_import
import importlib

import numpy as np

from shapegen import pointgen2d, utils
from mesh.meshgen import meshgen2d
from mesh.mesh import Mesh2D

def triangles():
    mesh = Mesh2D()
    mesh.bbox = [0, 0, 1, 1]
    mesh.set_boundary_points()
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.fixed_points = fixed_points
    mesh.points = np.array([[0,0], [0,0.5], [0,1], [0.5,1], [1,1], \
                            [1,0.5],[1,0],[0.5,0], [0.5, 0.5]])
    mesh.simplices = np.array([[0,8,1],
                                [1,8,2],
                                [2,8,3],
                                [3,8,4],
                                [4,8,5],
                                [5,8,6],
                                [6,8,7],
                                [7,8,0]])

    mesh.edges = np.array([[0,1],
                            [0,7],
                            [0,8],
                            [1,2],
                            [1,8],
                            [2,3],
                            [2,8],
                            [3,4],
                            [3,8],
                            [4,5],
                            [4,8],
                            [5,6],
                            [5,8],
                            [6,7],
                            [6,8],
                            [7,8]])

def ellipses():
    boundary_box = (0,0,1,1)
    l=0.07
    mesh = meshgen2d(boundary_box, l)
    ellipse1 = pointgen2d.sample_ellipse(0.4, 0.2, 10)
    ellipse2 = pointgen2d.sample_ellipse(0.2, 0.4, 10)
    shapes = [ellipse1, ellipse2]
    return mesh, mesh.simplices, mesh.edges, np.array(shapes)

def deform_example():
    boundary_box = (0,0,200,50)
    l=6
    mesh = meshgen2d(boundary_box, l)
    functions= ['curve4', 'curve5']
    points = list()
    for f in functions:
        points.append(pointgen2d.sample_function_mesh(mesh, f))
    return mesh, mesh.simplices, mesh.edges, np.array(points)

