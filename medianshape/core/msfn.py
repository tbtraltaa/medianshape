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

    Accepts simplicial settings, an input current, multiscale factor(:math:`\lambda`).
    Returns flat norm decomposition of the input current and the flat norm.
    Let K be an underlying simplicial complex of dimension q.

    :param float points: points in K.
    :param int simplices: (p+1)-simplices in K, an array of dimension (nx(p+1)) 
        where :math:`p+1 \leq q` and n is the number of p+1-simplices in K. 
    :param int subsimplices: p-simplices in K, an array of dimension (mxp) 
        where :math:`p \leq q` and m is the number of p-simplice in K. 
    :param int input_current: a vector for an input current. 
    :param float lambda_: multiscale factor.  
    :param float w: a vector of subsimplices volumes.
    :param float v: a vector of simplices volumes.
    :param int cons: a constraint matrix A.
    :returns: x, s, norm-- p-chain, (p+1)-chain of flat norm decompostion, flat norm.
    :rtype: int, int, float.
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
