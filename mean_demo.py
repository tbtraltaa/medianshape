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
    # l - initial length of triangle sides 
    # change it to 1 for big traingles
    mesh.points, mesh.simplices = distmesh2d("square", (0,0,1,1),[(0,0), (0,1), (1,0), (1,1)], l=1)
    #mesh.simplices = np.array([[2,1,0],[3,2,0],[4,2,3],[2,4,1]])
    #mesh.edges = np.array([[1,2],[0,1],[1,4],[2,4],[2,3],[0,3]])
    mesh.set_edges()
    mesh.to_string()
    mesh.orient_simplices_2D()
    #functions = ['func2']
    #functions = ['sin1pi']
    functions = ['x2', 'x5']
    fa = FunctionApprox2d(mesh)
    w = simpvol(mesh.points, mesh.edges)
    v = simpvol(mesh.points, mesh.simplices)
    input_currents = list()
    for f in functions:
        input_current = fa.generate_curve(f)
        csr_path = csr_matrix(input_current)
        print "Path vector:\n", csr_path
        mesh.plot()
        plt = fa.plot_curve()
        plt.show()
        input_currents.append(input_current)
    lambdas = [1, 10, 50]
    colors = ['red', 'yellow', 'blue']
    for l in lambdas:
        input_currents = list()
        current1 = np.zeros(shape=(8,1))
        current1[1] = 1
        current1[2] = 1
        current2 = np.zeros(shape=(8,1))
        current2[5] = 1
        current2[6] = 1
        input_currents.append(current1)
        input_currents.append(current2)
        input_currents = np.array(input_currents)
        x,q1,r1,q2,r2, norm = mean.mean(mesh.points, mesh.simplices, mesh.edges, input_currents, l)
        mesh.plot()
        for i, c in enumerate(input_currents):
            if i < 3:
                plt = mesh.plot_curve(c, color=colors[i])
            else:
                plt = mesh.plot_curve(c, color=colors[0])
        title = "lambda=%.01f"%l
        plt = mesh.plot_curve(x, title, "black")
        plt.show()
        print "Mean", norm
        print "x", x
        print "q1", q1
        print "r1", r1
        print "q2", q2
        print "r2", r2
    #mesh.orient_simplices()
