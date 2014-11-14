# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import math

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix

from mesh.distmesh import distmesh2d
from mesh.mesh import Mesh
from shape_gen.curve_gen import FunctionApprox2d
from mesh.utils import boundary_matrix, simpvol

import mean
from cvxopt import matrix, solvers

if __name__ == "__main__":
    mesh = Mesh()
    mesh.points, mesh.simplices = distmesh2d("square", (0,0,1,1),[(0,0), (0,1), (1,0), (1,1)])
    mesh.set_edges()
    mesh.to_string()
    mesh.orient_simplices_2D()
    functions = ['myfunc', 'x2', 'x5']
    #functions = ['x2']
    fa = FunctionApprox2d(mesh)
    w = simpvol(mesh.points, mesh.edges)
    v = simpvol(mesh.points, mesh.simplices)
    input_currents = list()
    for f in functions:
        input_current = fa.generate_curve(f)
        csr_path = csr_matrix(input_current)
        print "Path vector:\n", csr_path
        mesh.plot()
        fa.plot_curve()
        input_currents.append(input_current)
    lambdas = [0.1, 1, 10, 50]
    for l in lambdas:
        input_currents = np.array(input_currents)
        x, norm = mean.mean(mesh.points, mesh.simplices, mesh.edges, input_currents, l)
        mesh.plot()
        for c in input_currents:
            mesh.plot_curve(c, "red")
        plt = mesh.plot_curve(x, "black")
        plt.show()
        print "Mean", norm
    #mesh.orient_simplices()
