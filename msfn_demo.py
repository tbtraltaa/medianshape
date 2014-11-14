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

import msfn
from cvxopt import matrix, solvers

if __name__ == "__main__":
    mesh = Mesh()
    mesh.points, mesh.simplices = distmesh2d("square", (0,0,1,1),[(0,0), (0,1), (1,0), (1,1)])
    mesh.set_edges()
    mesh.to_string()
    mesh.orient_simplices_2D()
    functions = ['func2']
    #functions = ['x2']
    fa = FunctionApprox2d(mesh)
    boundary_matrix = boundary_matrix(mesh.simplices, mesh.edges)
    w = simpvol(mesh.points, mesh.edges)
    v = simpvol(mesh.points, mesh.simplices)
    m_edges = mesh.edges.shape[0]
    n_simplices = mesh.simplices.shape[0]
    cons = np.hstack((np.identity(m_edges), -np.identity(m_edges)))
    cons = np.hstack((cons, boundary_matrix))
    cons = np.hstack((cons, -boundary_matrix))
    print "Constaint matrix shape", cons.shape
    print "Cons", cons
    lambdas = [0.1, 1, 5, 10, 50]
    for l in lambdas:
        for f in functions:
            input_current = fa.generate_curve(f)
            csr_path = csr_matrix(input_current)
            print "Path vector:\n", csr_path
            mesh.plot()
            fa.plot_curve()
            x, s, norm = msfn.msfn(mesh.points, mesh.simplices, mesh.edges, input_current, l, cons=cons)
            mesh.plot_curve(x)

            print "MSFN", norm
    #mesh.orient_simplices()
