# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import math
import time

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix

from mesh.meshgen import distmesh2d
from mesh.mesh import Mesh2D
from mesh.utils import boundary_matrix, simpvol
from shapegen import pointgen2d, currentgen, utils
import plot2d

import msfn
import shortestpath_msfn

if __name__ == "__main__":
    start = time.time()
    fig = plt.figure(figsize=(19,8))
    mesh = Mesh2D()
    # l - initial length of triangle sides. Change it to vary traingle size
    mesh.bbox = (0,0,1,1)
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.set_boundary_points()
    mesh.points, mesh.simplices = distmesh2d('square', mesh.bbox, mesh.boundary_points, l=0.1)
    mesh.set_edges()
    mesh.orient_simplices_2D()

    points = list()
    points.append(pointgen2d.sample_function_mesh(mesh, 'sin1pi'))
    points = np.array(points)
    vertices, paths, input_currents = currentgen.push_functions_on_mesh_2d(mesh, points, functions=['sin1pi'])
    title = mesh.get_info()
    plot2d.plot_curves_approx(mesh, points, vertices, paths, title)
    plt.show()

    lambdas = [1]
    comb = [1]
    for input_current in input_currents:
        for l in lambdas:
            title = "lambda=%.04f"%l
            x, s, norm = msfn.msfn(mesh.points, mesh.simplices, mesh.edges, input_current, l)
            plot2d.plot_decomposition(mesh, input_currents, comb, x.T, None, s, title, lim=0.2)
            plt.show()
