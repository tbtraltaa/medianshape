# encoding: utf-8

from __future__ import absolute_import

import numpy as np

from mesh.utils import boundary_matrix, simpvol

import glpk
import pulp

from cvxopt import matrix, solvers

def print_cons(sub_cons, cons, c):
    string = ""
    print "Sub cons:"
    for row in sub_cons:
        for col in row:
            string += "%2d" % col
        print string
        string = ""

    string = ""
    print "\n"
    print "Cons"
    for row in cons:
        for col in row:
            string += "%2d" % col
        print string
        string = ""
    print "\n"
    print 'c', c

def mean(points, simplices, subsimplices, input_currents, lambda_, v=[], w=[], cons=[]):
    input_currents = np.array(input_currents)
    m_edges = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    k_currents = len(input_currents)
    input_currents  = input_currents.reshape(k_currents*m_edges,1)
    if w == []:
        w = simpvol(points, subsimplices)
        w[3] = 0.00001
        w[7] = 0.00001
    if v == []:
        v = simpvol(points, simplices)
    if cons == []:
        b_matrix = boundary_matrix(simplices, subsimplices)
        sub_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges), -b_matrix, b_matrix))
        sub_cons_col_count = 2*m_edges + 2*n_simplices
        identity_cons = np.hstack((np.identity(m_edges), -np.identity(m_edges)))
        k_identity_cons = np.tile(identity_cons,(k_currents,1))

        c = np.zeros(2*m_edges)
        sub_c = np.hstack((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v)))
        sub_c = sub_c.reshape(len(sub_c),1)
        k_sub_c = np.tile(sub_c, (k_currents,1))
        k_sub_c = k_sub_c/k_currents
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
    # Uncomment the line below to print sub_cons, cons and c
    #print_cons(sub_cons,cons, c)
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
    q1 = args[2*m_edges:3*m_edges] - args[3*m_edges:4*m_edges]
    r1 = args[4*m_edges:4*m_edges+n_simplices] - args[4*m_edges+n_simplices:4*m_edges+2*n_simplices]

    q2 = args[4*m_edges+2*n_simplices:5*m_edges+2*n_simplices] - args[5*m_edges+2*n_simplices:6*m_edges+2*n_simplices]
    r2 = args[6*m_edges+2*n_simplices:6*m_edges+3*n_simplices] - args[6*m_edges+3*n_simplices:6*m_edges+4*n_simplices]
    
    return x,q1,r1,q2,r2, norm