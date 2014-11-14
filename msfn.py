# encoding: utf-8

from __future__ import absolute_import

import numpy as np

from mesh.utils import boundary_matrix, simpvol

import glpk
import pulp

from cvxopt import matrix, solvers

def msfn(points, simplices, subsimplices, input_current, lambda_, v=[], w=[], cons=[]):
    m_edges = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    if w == []:
        w = simpvol(points, subsimplices)
    if v == []:
        v = simpvol(points, simplices)
    if cons == []:
        boundary_matrix = boundary_matrix(simplices, subsimplices)
        cons = np.hstack((np.identity(m_edges), -np.identity(m_edges)))
        cons = np.hstack((cons, boundary_matrix))
        cons = np.hstack((cons, -boundary_matrix))

    c = np.concatenate((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v))) 
    print "C", c
    c = c.reshape(len(c),1)
    print c.shape
    print cons.shape
    print input_current.shape
    c = matrix(c) 
    cons = matrix(cons)
    input_current = matrix(input_current)
    G = matrix(-np.identity(2*m_edges + 2*n_simplices))
    h = matrix(np.zeros(2*m_edges + 2*n_simplices))

    sol = solvers.lp(c, G, h, cons, input_current, solver='glpk')
    args = np.array(sol['x'])
    norm = sol['primal objective']
    x = args[0:m_edges] - args[m_edges:2*m_edges]
    s = args[2*m_edges:2*m_edges+n_simplices] - args[2*m_edges+n_simplices:]
    return x, s, norm
