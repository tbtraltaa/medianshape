from __future__ import absolute_import 

import numpy as np

from mesh.mesh import Mesh2D
import msfn
import plot2d
from scipy import sparse

import matplotlib.pyplot as plt

def envelope(mesh, input_currents):
    for i, c in enumerate(input_currents):
        comb=[1]
        if i < input_currents.shape[0] -1:
            diff = c - input_currents[i+1]
            x, s, norm = msfn.msfn(mesh.points, mesh.simplices, mesh.edges, diff, 0)
            plot2d.plot_mean(mesh, diff.reshape((1, len(c))), [1], [], lim=0.1)
            plt.show()
            plot2d.plot_decomposition(mesh, input_currents, comb, None, x, s, lim=0.1)
            plt.show()

def adjust_alphas(mesh, input_currents, t, v):
    alphas = list()
    total_sum = 0
    print input_currents.shape
    for i, c in enumerate(input_currents):
        comb=[1]
        diff = c.reshape(-1,1) - t
        x, s, norm = msfn.msfn(mesh.points, mesh.simplices, mesh.edges, diff, 0)
        plotting.plot_mean(mesh, diff.reshape((1, len(c))), [1], [])
        plt.show()
        plotting.plot_decomposition(mesh, input_currents, comb, None, x, s)
        plt.show()
        area =  np.sum(s*v)
        total_sum += area
        alphas.append(area)
    alphas = np.array(alphas).reshape(1,-1)
    alphas = (alphas *1.0/total_sum)/10000
    return alphas

# Saves sparse matrix as text. if the input is not sparse, set is_sparse argument to False.
def sparse_savetxt(fname, matrix, fmt='%d', delimiter=' '):
    if sparse.issparse(matrix):
        if matrix.getformat() !='coo':
            matrix = matrix.asformat('coo')
    else:
        matrix = sparse.coo_matrix(matrix)
    with open(fname, 'w') as f:
        for i in range(len(matrix.row)):
            f.write("%d %d %d\n" % (matrix.row[i], matrix.col[i], matrix.data[i]))

# Loads previously computed mesh, boundary_matrix, input currents, w and v from a directory
def load(dirname='/home/altaa/dumps1'):
    mesh = Mesh()
    mesh.points = np.loadtxt("%s/points.txt"%dirname) 
    mesh.simplices = np.loadtxt("%s/simplices.txt"%dirname)
    mesh.edges = np.loadtxt("%s/edges.txt"%dirname)
    v = np.loadtxt("%s/v.txt"%dirname)
    w = np.loadtxt("%s/w.txt"%dirname)
    b_matrix = sparse.dok_matrix((len(mesh.edges), len(mesh.simplices)), dtype=np.int8)
    with open("%s/b_matrix.txt"%dirname, 'r') as f:
        for line in f.readLines():
            data = line.split()
            b_matrix[int(data[0]), int(data[1])] = np.int8(data[2])
    return mesh, w, v, b_matrix

# Saves mesh, input currents, boundary matrix, w and v.
def save(mesh=None, input_currents=None, b_matrix=None, w=None, v=None, t=None, dirname='/home/altaa/dumps1', **kwargs):
    if mesh is not None:
        np.savetxt('%s/points.txt' % dirname, mesh.points, delimiter=' ')
        np.savetxt('%s/edges.txt' % dirname, mesh.edges, fmt='%d', delimiter=' ')
        np.savetxt('%s/simplices.txt'% dirname, mesh.simplices, fmt='%d', delimiter=' ')
        if mesh.points.shape[1] == 3:
            np.savetxt('%s/triangles.txt'% dirname, mesh.triangles, fmt='%d', delimiter=' ')
    if input_currents is not None:
        for i, c in enumerate(input_currents):
            sparse_savetxt('%s/input_current%d.txt' % (dirname,i), c)
    if b_matrix is not None:
        sparse_savetxt('%s/b_matrix.txt' % dirname, b_matrix)
    if w is not None:
        np.savetxt('%s/w.txt' % dirname, w, delimiter=' ')
    if v is not None:
        np.savetxt('%s/v.txt' % dirname, v, delimiter=' ')
    if t is not None:
        if 'opt' in kwargs and 'lambda_' in kwargs:
            sparse_savetxt("%s/t-%s-lambda-%s.txt"%(dirname, kwargs['opt'], kwargs['lambda_']), t)
        else:
            sparse_savetxt("%s/t.txt"%dirname, t)
