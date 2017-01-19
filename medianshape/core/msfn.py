#encoding: utf-8
'''
MSFN - Multiscale flat norm
===========================
'''

from __future__ import absolute_import

import numpy as np
from scipy import sparse
from cvxopt import matrix, spmatrix, solvers

from medianshape import utils

def msfn(points, simplices, subsimplices, input_current, lambda_, w=[], v=[], cons=[]):
    '''
    MSFN - Multiscale flat norm
    '''
    m_subsimplices = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    if w == []:
        w = utils.simpvol(points, subsimplices)
    if v == []:
        v = utils.simpvol(points, simplices)
    if cons == []:
        b_matrix = utils.boundary_matrix(simplices, subsimplices, format='coo')
        m_subsimplices_identity = sparse.identity(m_subsimplices, dtype=np.int8, format='coo')
        cons = sparse.hstack((m_subsimplices_identity, -m_subsimplices_identity, b_matrix, -b_matrix))
    c = np.concatenate((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v))) 
    c = c.reshape(len(c),1)
    c = matrix(c) 
    cons = spmatrix(cons.data.tolist(), cons.row, cons.col, cons.shape, tc='d')
    input_current = matrix(input_current)
    g = -sparse.identity(2*m_subsimplices + 2*n_simplices, dtype=np.int8, format='coo')
    h = np.zeros(2*m_subsimplices + 2*n_simplices)
    G = spmatrix(g.data.tolist(), g.row, g.col, g.shape,  tc='d')
    h = matrix(h)

    sol = solvers.lp(c, G, h, cons, input_current, solver='glpk')
    args = np.rint(sol['x'])
    norm = sol['primal objective']
    x = (args[0:m_subsimplices] - args[m_subsimplices:2*m_subsimplices]).reshape((1,m_subsimplices)).astype(int)
    s = (args[2*m_subsimplices:2*m_subsimplices+n_simplices] - args[2*m_subsimplices+n_simplices:]).reshape(1, n_simplices).astype(int)
    return x, s, norm
