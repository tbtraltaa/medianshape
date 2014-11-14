# encoding: utf-8

from __future__ import absolute_import

import numpy as np

from mesh.utils import boundary_matrix, simpvol

import glpk
import pulp

from cvxopt import matrix, solvers

def mean(points, simplices, subsimplices, input_currents, lambda_, v=[], w=[], cons=[]):
    m_edges = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    k_currents = len(input_currents)
    input_currents.reshape(k_currents*m_edges,1)
    if w == []:
        w = simpvol(points, subsimplices)
    if v == []:
        v = simpvol(points, simplices)
    if cons == []:
        boundary_matrix = boundary_matrix(simplices, subsimplices)
        sub_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges)))
        sub_cons = np.hstack((sub_cons, -boundary_matrix))
        sub_cons = np.hstack((sub_cons, boundary_matrix))
        sub_cons_col_count = 2*m_edges + 2*n_simplices
        identity_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges)))
        k_identity_cons = np.tile(identity_cons,(k_currents,1))

        c = np.zeros(2*m_edges)
        sub_c = np.concatenate((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v)))
        k_sub_c = np.tile(sub_c, (k_currents,1))
        c = np.concatenate((c, k_sub_c))
        for i in range(0,k_currents):
            cons_row = np.zeros((m_edges, k_currents*(2*m_edges + 2*n_simplices)))
            cons_row[:, 0:2*m_edges] = identity_cons
            sub_cons_start = 2*m_edges + i*sub_cons_col_count
            sub_cons_end = sub_cons_start + sub_cons_col_count
            cons_row[:, sub_cons_start:sub_cons_end] = sub_cons
            if i == 0:
                cons = cons_row
            else:
                cons = np.vstack((cons, cons_row))
        cons = np.hstack((k_identity_cons, cons))

    print "C", c
    print c.shape
    print cons.shape
    print input_currents.shape
    c = matrix(c) 
    cons = matrix(cons)
    input_currents = matrix(input_currents)
    cons_col_count = 2*m_edges + k_currents*(2*m_edges + 2*n_simplices)
    G = matrix(-np.identity(cons_col_count))
    h = matrix(np.zeros(cons_col_count))

    sol = solvers.lp(c, G, h, cons, input_currents, solver='glpk')
    args = np.array(sol['x'])
    norm = sol['primal objective']
    x = args[0:m_edges] - args[m_edges:2*m_edges]
    return x, norm
