# encoding: utf-8
'''
MSFN demo 2D
++++++++++++
'''

from __future__ import absolute_import

import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

from medianshape.simplicial.meshgen import distmesh2d
from medianshape.simplicial.mesh import Mesh2D
from medianshape.simplicial import pointgen2d, currentgen
from medianshape.viz import plot2d
from medianshape.core.msfn import msfn

def msfndemo2d():
    '''
    MSFN demo 2D
    '''
    start = time.time()
    fig = plt.figure(figsize=(8,8))
    mesh = Mesh2D()
    # l - initial length of triangle sides. Change it to vary traingle size
    mesh.bbox = (0,0,1,1)
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.set_boundary_points()
    mesh.points, mesh.simplices = distmesh2d('square', mesh.bbox, mesh.boundary_points, l=0.06)
    mesh.set_edges()
    mesh.orient_simplices_2D()

    points = list()
    points.append(pointgen2d.sample_function('sin1pi', mesh.bbox, 20))
    points = np.array(points)
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, mesh.simplices, mesh.edges, points)
    title = mesh.get_info()
    plot2d.plot_curves_approx(mesh, points, vertices, paths, title)
    plt.show()
    fig = plt.figure(figsize=(8,8))

    lambdas = [1, 2, 3, 10, 20 ]
    comb = [1]
    for input_current in input_currents:
        for l in lambdas:
            title = "lambda=%.04f"%l
            x, s, norm = msfn(mesh.points, mesh.simplices, mesh.edges, input_current, l)
            plot2d.plot_decomposition(mesh, input_currents, comb, x.T, None, s, lim = 0.2)
            plt.title(title)
            plt.show()
            fig = plt.figure(figsize=(8,8))

if __name__ == "__main__":
    msfndemo2d()
