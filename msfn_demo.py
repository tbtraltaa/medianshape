# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import math
import time

import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix

from mesh.distmesh import distmesh2d
from mesh.mesh import Mesh
from mesh.utils import boundary_matrix, simpvol
from shape_gen import point_gen, curve_gen, utils
import plotting

import msfn

if __name__ == "__main__":
    start = time.time()
    fig = plt.figure(figsize=(19,8))
    mesh = Mesh()
    # l - initial length of triangle sides. Change it to vary traingle size
    mesh.boundary_box = (0,0,200,50)
    mesh.fixed_points = [(0,0),(200,0),(0,50),(200,50)]
    mesh.points, mesh.simplices = distmesh2d('square', mesh.boundary_box, mesh.fixed_points, l=20)
    mesh.set_edges()
    mesh.orient_simplices_2D()

    functions = ['curve1']
    points, vertices, paths, input_currents = curve_gen.push_curves_on_mesh(mesh, functions)
    title = 'Functions - %s - (%s)' % (mesh.get_info(), ','.join(functions))
    plotting.plot_curves_approx(mesh, points, vertices, paths, title)
    plt.show()

    lambdas = [0.0001]
    comb = [1]
    for input_current in input_currents:
        for l in lambdas:
            title = "lambda=%.04f"%l
            x, s, norm = msfn.msfn(mesh.points, mesh.simplices, mesh.edges, input_current, l)
            plotting.plot_decomposition(mesh, functions, input_currents, comb, None, x, s, title)
            plt.show()
