# encoding: utf-8

from __future__ import absolute_import

import math
import importlib
import random
import itertools


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from scipy.sparse import csr_matrix

from mesh.distmesh import distmesh2d
from mesh.mesh import Mesh
from shape_gen.curve_gen import FunctionApprox2d
from mesh.utils import boundary_matrix, simpvol

import mean
from cvxopt import matrix, solvers

options = {   'default': 1,
            'mass': 2,
            'msfn': 3
        }

if __name__ == "__main__":
    mesh = Mesh()
    # l - initial length of triangle sides 
    # change it to 1 for big traingles
    mesh.points, mesh.simplices = distmesh2d("square", (0,0,1,1),[(0,0), (0,1), (1,0), (1,1)])
    mesh.set_edges()
#    mesh.points = np.array([[0,0], [0,0.5],[0,1],[0.5,1],[1,1],[1,0.5],[1,0],[0.5,0], [0.5, 0.5]])
#    mesh.simplices = np.array([[0,8,1],
#                                [1,8,2],
#                                [2,8,3],
#                                [3,8,4],
#                                [4,8,5],
#                                [5,8,6],
#                                [6,8,7],
#                                [7,8,0]])
#
#    mesh.edges = np.array([[0,1],
#                            [0,7],
#                            [0,8],
#                            [1,2],
#                            [1,8],
#                            [2,3],
#                            [2,8],
#                            [3,4],
#                            [3,8],
#                            [4,5],
#                            [4,8],
#                            [5,6],
#                            [5,8],
#                            [6,7],
#                            [6,8],
#                            [7,8]])
    mesh.to_string()
    mesh.orient_simplices_2D()
    #functions = ['sin1pi', 'sin1pi1']
    #functions = ['x', 'x2', 'x5']
    function_sets = [['x','x2','x5'], ['sin1pi', 'half_sin1pi'], ['sin1pi', 'small_sin1pi']]
    fig_count = 1
    for j, functions in enumerate(function_sets):
        if len(functions) == 2:
            color_set = 'gr'
        else:
            color_set = 'gry'
        colors = itertools.cycle(color_set)
        fa = FunctionApprox2d(mesh)
        input_currents = list()
        fig = plt.figure(figsize=(10,10))
        plt.gca().set_aspect('equal')
        plt.ylim([-0.2, 1.2])
        plt.xlim([-0.2, 1.2])
        mesh.plot()
        for i, f in enumerate(functions):
            input_current = fa.generate_curve(f)
            csr_path = csr_matrix(input_current)
            print "Path vector:\n", csr_path
            fa.plot_curve(color=colors.next())
            input_currents.append(input_current)
        k_currents = len(functions)
        plt.title("Functions")
        plt.savefig("/home/altaa/%d-%s.png"%(fig_count, "-".join(functions)), dpi=fig.dpi)
        fig_count += 1
        input_currents = np.array(input_currents)
        for opt in options:
            lambdas = [0.01, 1, 5, 50]
            for l in lambdas:
#            input_currents = list()
#            current1 = np.zeros(shape=(len(mesh.edges),1))
#            current1[0] = 1
#            current1[3] = 1
#            current1[5] = 1 
#            current1[7] = 1
#            current2 = np.zeros(shape=(len(mesh.edges),1))
#            current2[1] = 1
#            current2[13] = -1
#            current2[11] = -1
#            current2[9] = -1
#            input_currents.append(current1)
#            input_currents.append(current2)
#            input_currents = np.array(input_currents)
#            k_currents = 2
                x, q, r, norm = mean.mean(mesh.points, mesh.simplices, mesh.edges, input_currents, l, options[opt])
                fig = plt.figure(figsize=(10,10))
                plt.gca().set_aspect('equal')
                #plt.ylim([-0.2, 1.2])
                #plt.xlim([-0.2, 1.2])
                #plt.axes([-0.2,-0.2, 1.4,1.4])
                cols = 2
                if opt == 'msfn' and len(functions) == 3:
                    cols = 3
                plt.subplot(cols, 2, 1)
                plt.ylim([-0.2, 1.2])
                plt.xlim([-0.2, 1.2])
                mesh.plot()
                for i, c in enumerate(input_currents):
                    mesh.plot_curve(c, color=colors.next())
                title = "%s, lambda=%.02f"%(opt, l)
                mesh.plot_curve(x, title)
                colors = itertools.cycle(color_set)
                for i in range(r.shape[1]):
                    plt.subplot(cols, 2, 2+i)
                    color = colors.next()
                    #plt.figure(figsize=(10,15))
                    #plt.gca().set_aspect('equal')
                    plt.ylim([-0.2, 1.2])
                    plt.xlim([-0.2, 1.2])
                    mesh.plot()
                    mesh.plot_curve(x)
                    mesh.plot_simplices(r[:,i], color=color)
                    mesh.plot_curve(q[:,i], title=title + ", Q%d&R%d"%(i+1,i+1), color="m", marker='*')
                plt.savefig("/home/altaa/%d-%s-%s-%.02f.png"%(fig_count, "-".join(functions),opt,l), dpi=fig.dpi)
                fig_count += 1
                    #print "q1", q1
                    #print "r1", r1
                    #print "q2", q2
                    #print "r2", r2
