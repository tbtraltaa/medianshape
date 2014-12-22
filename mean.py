# encoding: utf-8

from __future__ import absolute_import

import numpy as np

from mesh.utils import boundary_matrix, simpvol

import glpk
import pulp

from cvxopt import matrix, solvers
solvers.options['abstol'] = 1e-10
solvers.options['reltol'] = 1e-9
solvers.options['feastol'] = 1e-10

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

def mean(points, simplices, subsimplices, input_currents, lambda_, opt='default', w=[], v=[], cons=[]):
    input_currents = np.array(input_currents)
    m_edges = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    k_currents = len(input_currents)
    sub_cons_count = k_currents 
    input_currents  = input_currents.reshape(k_currents*m_edges,1)
    if opt == 'msfn':
        input_currents = np.vstack((input_currents, np.zeros((m_edges,1))))
        sub_cons_count = k_currents + 1 
    b = input_currents
    print "b", b.shape
    if w == []:
        #w = simpvol(points, subsimplices)
        #w = np.ones(m_edges)
        w = np.ndarray(shape=(m_edges,))
        #w[0:] = 0.09
        w[0:] = 1
        #w[0:] = 0.09050288
        print "w", w[0:10]
    if v == []:
        #v = simpvol(points, simplices)
        v =  np.ndarray(shape=(n_simplices,))
        v[0:] = 0.433
        #v[0:] = 0.003
        #v[0:] = 0.00395007
        print "v", v[0:10]
    if cons == []:
        # Msfn option
        if opt == 'msfn':
            sub_cons_count += 1
        # Mass option
        b_matrix = boundary_matrix(simplices, subsimplices)
        identity_cons = np.hstack((np.identity(m_edges), -np.identity(m_edges)))
        sub_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges), -b_matrix, b_matrix))
        sub_cons_col_count = 2*m_edges + 2*n_simplices
        k_identity_cons = np.tile(identity_cons,(sub_cons_count,1))

        for i in range(0,sub_cons_count):
            cons_row = np.zeros((m_edges, sub_cons_count*(2*m_edges + 2*n_simplices)))
            sub_cons_start = i*sub_cons_col_count
            sub_cons_end = sub_cons_start + sub_cons_col_count
            cons_row[:, sub_cons_start:sub_cons_end] = sub_cons
            if i == 0:
                cons = cons_row
            else:
                cons = np.vstack((cons, cons_row))

        cons = np.hstack((k_identity_cons, cons))
        
    c = np.zeros((2*m_edges,1))
    if opt == 'mass':
        c = np.vstack((abs(w),abs(w)))
    sub_c = np.hstack((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v)))
    sub_c = sub_c.reshape(len(sub_c),1)
    k_sub_c = np.tile(sub_c, (sub_cons_count,1))
    c = np.append(c, k_sub_c)
    c = c/k_currents


    # Uncomment the line below to print sub_cons, cons and c
    #print_cons(sub_cons,cons, c)
    
    #np.savetxt("/home/altaa/dumps1/b-%s.txt"%opt, input_currents, fmt="%d", delimiter=" ")
    #np.savetxt("/home/altaa/dumps1/c-%s.txt"%opt, c, delimiter=" ")
    G = matrix(-np.identity(len(c)))
    h = matrix(np.zeros(len(c)))
    c = matrix(c) 
    cons = matrix(cons)
    b = matrix(b)
    print "w", w[0:10]
    print "v", v[0:10]

    sol = solvers.lp(c, G, h, cons, b, solver='glpk')
    args = np.array(sol['x'])
    norm = sol['primal objective']
    x = args[0:m_edges] - args[m_edges:2*m_edges]
    q = np.zeros((m_edges, sub_cons_count))
    r = np.zeros((n_simplices, sub_cons_count))
    qi_start = 2*m_edges
    for i in range(sub_cons_count):
        qi_end = qi_start + 2*m_edges
        q[:, i] = (args[qi_start: qi_start+m_edges] - args[qi_start+m_edges: qi_end]).reshape(m_edges)
        ri_start = qi_end
        ri_end = ri_start + 2*n_simplices
        r[:, i] = (args[ri_start: ri_start+n_simplices] - args[ri_start+n_simplices: ri_end]).reshape(n_simplices)
        qi_start = ri_end
    return x, q, r, norm

def get_lp_inputs(points, simplices, subsimplices, k_currents, opt='default', w=[], v=[], cons=[]):
    m_edges = subsimplices.shape[0]
    n_simplices = simplices.shape[0]
    if w == []:
        #w = simpvol(points, subsimplices)
        #w = np.ones(m_edges)
        w = np.ndarray(shape=(m_edges,))
        #w[0:] = 0.09
        w[0:] = 1
        #w[0:] = 0.09050288
        print "w", w[0:10]
    if v == []:
        #v = simpvol(points, simplices)
        v =  np.ndarray(shape=(n_simplices,))
        v[0:] = 0.433
        #v[0:] = 0.003
        #v[0:] = 0.00395007
        print "v", v[0:10]
    if cons == []:
        sub_cons_count = k_currents
        # Msfn option
        if opt == 'msfn':
            sub_cons_count += 1
        # Mass option
        b_matrix = boundary_matrix(simplices, subsimplices)
        identity_cons = np.hstack((np.identity(m_edges), -np.identity(m_edges)))
        sub_cons = np.hstack((-np.identity(m_edges), np.identity(m_edges), -b_matrix, b_matrix))
        sub_cons_col_count = 2*m_edges + 2*n_simplices
        k_identity_cons = np.tile(identity_cons,(sub_cons_count,1))


        for i in range(0,sub_cons_count):
            cons_row = np.zeros((m_edges, sub_cons_count*(2*m_edges + 2*n_simplices)))
            sub_cons_start = i*sub_cons_col_count
            sub_cons_end = sub_cons_start + sub_cons_col_count
            cons_row[:, sub_cons_start:sub_cons_end] = sub_cons
            if i == 0:
                cons = cons_row
            else:
                cons = np.vstack((cons, cons_row))

        cons = np.hstack((k_identity_cons, cons))
    return w, v, cons, b_matrix
