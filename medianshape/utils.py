
from __future__ import absolute_import 

import numpy as np

from scipy import sparse

from mesh.mesh import Mesh2D, Mesh3D


# Loads previously computed mesh, boundary_matrix, input currents, w and v from a directory
def load_mesh2d(dirname='dumps'):
    mesh = Mesh2D()
    mesh.bbox = np.loadtxt("%s/bbox.txt"%dirname)
    mesh.set_boundary_points()
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.points = np.loadtxt("%s/points.txt"%dirname) 
    mesh.edges = np.loadtxt("%s/edges.txt"%dirname)
    mesh.simplices = np.loadtxt("%s/simplices.txt"%dirname)
    return mesh

def load_mesh3d(dirname='dumps'):
    mesh = Mesh3D()
    mesh.bbox = np.loadtxt("%s/bbox.txt"%dirname)
    mesh.set_boundary_points()
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.points = np.loadtxt("%s/points.txt"%dirname) 
    mesh.edges = np.loadtxt("%s/edges.txt"%dirname)
    mesh.triangles = np.loadtxt("%s/triangles.txt"%dirname)
    mesh.simplices = np.loadtxt("%s/simplices.txt"%dirname)
    return mesh

def load_weights_and_boundary(dirname='dumps'):
    v = np.loadtxt("%s/v.txt"%dirname)
    w = np.loadtxt("%s/w.txt"%dirname)
    b_matrix = sparse.dok_matrix((len(mesh.edges), len(mesh.simplices)), dtype=np.int8)
    with open("%s/b_matrix.txt"%dirname, 'r') as f:
        for line in f.readLines():
            data = line.split()
            b_matrix[int(data[0]), int(data[1])] = np.int8(data[2])
    return w, v, b_matrix

# Saves mesh, input currents, boundary matrix, w and v.
def save_data(mesh=None, input_currents=None, b_matrix=None, w=None, v=None, t=None, dirname='dumps', **kwargs):
    if mesh is not None:
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
        if 'mu' in kwargs and 'lambda_' in kwargs:
            sparse_savetxt("%s/t-lambda-%s-mu-%s.txt"%(dirname, kwargs['lambda_'], kwargs['mu']), t)
        else:
            sparse_savetxt("%s/t.txt"%dirname, t)
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
