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
    input_currents  = input_currents.reshape(k_currents*m_edges,1)
    if w == []:
        w = simpvol(points, subsimplices)
    if v == []:
        v = simpvol(points, simplices)
    if cons == []:
        b_matrix = boundary_matrix(simplices, subsimplices)
        sub_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges), -b_matrix, b_matrix))
        sub_cons_col_count = 2*m_edges + 2*n_simplices
        identity_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges)))
        k_identity_cons = np.tile(identity_cons,(k_currents,1))

        c = np.zeros(2*m_edges)
        sub_c = np.hstack((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v)))
        sub_c = sub_c.reshape(len(sub_c),1)
        k_sub_c = np.tile(sub_c, (k_currents,1))
        c = np.append(c, k_sub_c)
        for i in range(0,k_currents):
            cons_row = np.zeros((m_edges, k_currents*(2*m_edges + 2*n_simplices)))
            sub_cons_start = i*sub_cons_col_count
            sub_cons_end = sub_cons_start + sub_cons_col_count
            cons_row[:, sub_cons_start:sub_cons_end] = sub_cons
            if i == 0:
                cons = cons_row
            else:
                cons = np.vstack((cons, cons_row))
        cons = np.hstack((k_identity_cons, cons))

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
