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
    mesh.points, mesh.simplices = distmesh2d("square", (0,0,1,1),[(0,0), (0,1), (1,0), (1,1)])
    mesh.set_edges()

#    mesh.simplices = np.array([[2,1,0],[3,2,0],[4,2,3],[2,4,1]])
#    mesh.edges = np.array([[1,2],[0,1],[1,4],[2,4],[2,3],[0,3]])
#    mesh.points = np.array([[0,0], [0,1],[1,1],[1,0],[0, 0.5],[0.5,1],[1,0.5],[0.5,0], [0.5, 0.5]])
#    mesh.simplices = np.array([[0,8,4],
#                                [4,8,1],
#                                [1,8,5],
#                                [8,2,5],
#                                [8,6,2],
#                                [8,3,6],
#                                [8,7,3],
#                                [8,0,7]])
#    mesh.edges = np.array([[0,4],
#                            [1,4],
#                            [1,5],
#                            [2,5],
#                            [2,6],
#                            [3,6],
#                            [3,7],
#                            [0,7],
#                            [0,8],
#                            [4,8],
#                            [1,8],
#                            [5,8],
#                            [2,8],
#                            [6,8],
#                            [3,8],
#                            [7,8]])
    mesh.to_string()
    mesh.orient_simplices_2D()
    #functions = ['func2']
    #functions = ['sin1pi']
    functions = ['myfunc', 'x2', 'x5']
    fa = FunctionApprox2d(mesh)
    input_currents = list()
    for f in functions:
        input_current = fa.generate_curve(f)
        csr_path = csr_matrix(input_current)
        print "Path vector:\n", csr_path
        mesh.plot()
        plt = fa.plot_curve()
        plt.show()
        input_currents.append(input_current)
    input_currents = np.array(input_currents)
    lambdas = [50, 1, 0.1, 0.01]
    colors = ['red', 'yellow', 'pink']
    for l in lambdas:
#        input_currents = list()
#        current1 = np.zeros(shape=(len(mesh.edges),1))
#        current1[1] = -1
#        current1[2] = -1
#        current2 = np.zeros(shape=(len(mesh.edges),1))
#        current2[5] = 1
#        current2[6] = 1
#        input_currents.append(current1)
#        input_currents.append(current2)
#        input_currents = np.array(input_currents)
        x,q1,r1,q2,r2, norm = mean.mean(mesh.points, mesh.simplices, mesh.edges, input_currents, l)
        plt.figure(facecolor="white", edgecolor=None)
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
#        print "q1", q1
#        print "r1", r1
#        print "q2", q2
#        print "r2", r2
