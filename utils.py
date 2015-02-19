from __future__ import absolute_import 

import numpy as np

from mesh.mesh import Mesh
from scipy.sparse import dok_matrix, coo_matrix

def extract_edges(simplices):
    edges = set()
    for simplex in simplices:
        for i in range(len(simplex)):
            edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
    return np.array(list(edges))

# Saves sparse matrix as text. if the input is not sparse, set is_sparse argument to False.
def sparse_savetxt(fname, matrix, fmt='%d', delimiter=' ', is_sparse=True):
    if not is_sparse:
        matrix = coo_matrix(matrix)
    else:
        matrix = matrix.asformat('coo')
    with open(fname, 'w') as f:
        for i in range(len(matrix.row)):
            f.write("%d %d %d\n" % (matrix.row[i], matrix.col[i], matrix.data[i]))

# Loads previously computed mesh, boundary_matrix, input currents, w and v from a directory
def load(dirname='/home/altaa/dumps1'):
    mesh = Mesh()
    mesh.points = np.loadtxt("%s/points.txt"%dirname) 
    mesh.simplices = np.loadtxt("%s/points.txt"%dirname)
    mesh.edges = np.loadtxt("%s/edges.txt"%dirname)
    v = np.loadtxt("%s/v.txt"%dirname)
    w = np.loadtxt("%s/w.txt"%dirname)
    b_matrix = dok_matrix((len(mesh.edges), len(mesh.simplices)), dtype=np.int8)
    with open("%s/b_matrix.txt"%dirname, 'r') as f:
        for line in f.readLines():
            data = line.split()
            b_matrix[int(data[0]), int(data[1])] = np.int8(data[2])
    
    input_currents = []
    for i in range(1,4):
        curve = np.zero((len(mesh.edges), 1), dtype=np.int8)
        with open("%s/curve%d.txt"%(dirname,i), 'r') as f:
            for line in f.readLines():
                data = line.split()
                curve[int(data[0]), int(data[1])] = np.int8(data[2])
            if input_currents:
                input_currents = np.vstack(input_currents, curve)
            else:
                input_currents = curve
    return mesh, input_currents, b_matrix, w, v

# Saves mesh, input currents, boundary matrix, w and v.
def save(mesh=None, input_currents=None, b_matrix=None, w=None, v=None, t=None, dirname='/home/altaa/dumps1', **kwargs):
    if mesh is not None:
        np.savetxt('%s/points.txt' % dirname, mesh.points, delimiter=' ')
        np.savetxt('%s/edges.txt' % dirname, mesh.edges, fmt='%d', delimiter=' ')
        np.savetxt('%s/simplices.txt'% dirname, mesh.simplices, fmt='%d', delimiter=' ')
    if input_currents is not None:
        print input_currents.shape
        for i, c in enumerate(input_currents):
            sparse_savetxt('%s/input_current%d.txt' % (dirname,i), c, is_sparse=False)
    if b_matrix is not None:
        sparse_savetxt('%s/b_matrix.txt' % dirname, b_matrix)
    if w is not None:
        np.savetxt('%s/w.txt' % dirname, w, delimiter=' ')
    if v is not None:
        np.savetxt('%s/v.txt' % dirname, v, delimiter=' ')
    if t is not None:
        if 'opt' in kwargs and 'lambda_' in kwargs:
            np.savetxt("%s/t-%s-lambda-%s.txt"%(dirname, kwargs['opt'], kwargs['lambda_']), t, \
            fmt="%d", delimiter=" ")
        else:
            np.savetxt("%s/t.txt" % dirname, t, fmt="%d", delimiter=" ")

