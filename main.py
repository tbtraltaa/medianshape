# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import math

import numpy as np
import matplotlib.pyplot as plt

from mesh.distmesh import distmesh2d
from mesh.mesh import Mesh
from shape_gen.curve_gen import FunctionApprox2d

if __name__ == "__main__":
    mesh = Mesh()
    mesh.points, mesh.simplices = distmesh2d("square", (0,0,1,1),[(0,0), (0,1), (1,0), (1,1)])
    mesh.set_edges()
    mesh.to_string()
    functions = ['myfunc', 'x2', 'x3', 'math.atan', 'math.acos']
    functions = ['myfunc']
    fa = FunctionApprox2d(mesh)
    mesh.plot()
    for f in functions:
        fa.generate_curve(f)
    mesh.orient_simplices_2D()
    #mesh.orient_simplices()
