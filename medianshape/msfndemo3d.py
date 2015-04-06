# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import math
import time

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix

from mesh.distmesh import distmesh3d
from mesh.mesh import Mesh3D
from mesh.utils import boundary_matrix, simpvol
from shape_gen import point_gen, curve_gen, utils
import plotting

import distmesh as dm

import msfn

if __name__ == "__main__":
    start = time.time()
    fig = plt.figure(figsize=(19,8))
    mesh = Mesh3D()
    # l - initial length of triangle sides. Change it to vary traingle size
    mesh.bbox = (0,0,0, 50, 50, 50)
    mesh.set_fixed_points()
    mesh.set_diagonal()
    mesh.points, mesh.simplices= distmesh3d("ball", mesh.bbox, None, l=10)
    mesh.edges = dm.mkt2t(mesh.simplices)
    print mesh.simplices
    print mesh.simplices.shape
    print mesh.get_info()
    print mesh.edges
    #mesh.set_edges()
    #mesh.orient_simplices_2D()
    exit()

    functions = ['curve1']
    points = list()
    points.append(point_gen.sample_function_mesh(mesh, 'curve1'))
    points = np.array(points)
    vertices, paths, input_currents = curve_gen.push_curves_on_mesh(mesh, points)
    title = 'Functions - %s - (%s)' % (mesh.get_info(), ','.join(functions))
    plotting.plot_curves_approx(mesh, points, vertices, paths, title)
    plt.show()

    lambdas = [1]
    comb = [1]
    for input_current in input_currents:
        for l in lambdas:
            title = "lambda=%.04f"%l
            input_current = np.zeros(len(mesh.points))
            input_current[2] = 1
            input_current[20] = 1
            point_indices = np.arange(mesh.points.shape[0], dtype=int).reshape(-1, 1)
            x, s, norm = msfn.msfn(mesh.points, mesh.edges, point_indices, input_current, l)
            plotting.plot_decomposition(mesh, input_currents, comb, None, x, s, title, r_dim=1)
            plt.show()
