# encoding: utf-8

from __future__ import absolute_import

import numpy as np
from scipy import sparse
from cvxopt import matrix, spmatrix, solvers

from mesh.utils import boundary_matrix, simpvol

def msfn(points, simplices, subsimplices, input_current, lambda_, w=[], v=[], cons=[]):
    m_edges = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    if w == []:
        w = simpvol(points, subsimplices)
    if v == []:
        v = simpvol(points, simplices)
    if cons == []:
        b_matrix = boundary_matrix(simplices, subsimplices, format='coo')
        cons = sparse.hstack((b_matrix, -b_matrix))
    c = np.concatenate((lambda_*abs(v), lambda_*abs(v))) 
    c = c.reshape(len(c),1)
    c = matrix(c) 
    cons = spmatrix(cons.data.tolist(), cons.row, cons.col, cons.shape, tc='d')
    input_current = matrix(input_current)
    g = -sparse.identity(2*n_simplices, dtype=np.int8, format='coo')
    h = np.zeros(2*n_simplices)
    G = spmatrix(g.data.tolist(), g.row, g.col, g.shape,  tc='d')
    h = matrix(h)

    sol = solvers.lp(c, G, h, cons, input_current, solver='glpk')
    args = np.array(sol['x'])
    print args
    #args = np.rint(sol['x'])
    norm = sol['primal objective']
    s = (args[0:n_simplices] - args[n_simplices:]).reshape(1, n_simplices).astype(int)
    print b_matrix
    print x
    print s
    return [], s, norm
