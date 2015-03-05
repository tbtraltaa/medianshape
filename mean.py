# encoding: utf-8

from __future__ import absolute_import
import time

import numpy as np
from scipy import sparse

from cvxopt import matrix, spmatrix, solvers
solvers.options['abstol'] = 1e-10
solvers.options['reltol'] = 1e-9
solvers.options['feastol'] = 1e-10
solvers.options['show_progress'] = False

from mesh.utils import boundary_matrix, simpvol

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

def mean(mesh, input_currents, lambda_, opt='default', w=[], v=[], cons=[], mu=0.001, alpha=None, len_cons=False):
    if not isinstance(input_currents, np.ndarray):
        input_currents = np.array(input_currents)
    average_len = np.rint(np.average(np.array([c.nonzero()[0].shape[0] for c in input_currents])))
    print "Average len", average_len
    m_edges = mesh.edges.shape[0]
    n_simplices = mesh.simplices.shape[0]
    k_currents = len(input_currents)
    sub_cons_count = k_currents 
    input_currents  = input_currents.reshape(k_currents*m_edges,1)
    w, v, b_matrix, cons = get_lp_inputs(mesh, k_currents, opt, w, v, [], cons)
    if opt == 'msfn':
        input_currents = np.vstack((input_currents, np.zeros((m_edges,1))))
        sub_cons_count += 1 
    b = input_currents
    print "b", b.shape
    c = np.zeros((2*m_edges,1))
    if opt == 'mass':
        c = np.vstack((abs(w),abs(w)))
        c = c*mu
    sub_c = np.hstack((abs(w), abs(w), lambda_*abs(v), lambda_*abs(v)))
    sub_c = sub_c.reshape(len(sub_c),1)
    k_sub_c = np.tile(sub_c, (sub_cons_count,1))
    if alpha is not None:
        for i in range(k_currents):
            if i < k_currents+1:
                k_sub_c[i*(2*m_edges+2*n_simplices):(i+1)*(2*m_edges+2*n_simplices)] = \
                k_sub_c[i*(2*m_edges+2*n_simplices):(i+1)*(2*m_edges+2*n_simplices)]*alpha[i]
                print "Alpha", i, alpha[i]
    c = np.append(c, k_sub_c)
    if opt != 'mass':
        c = c/k_currents
    #np.savetxt("/home/altaa/dumps1/b-%s.txt"%opt, input_currents, fmt="%d", delimiter=" ")
    #np.savetxt("/home/altaa/dumps1/c-%s.txt"%opt, c, delimiter=" ")
    print "Size of c: ", len(c)
    print "Average len", np.rint(average_len)

    g = -sparse.identity(len(c), dtype=np.int8, format='coo')
    h = np.zeros(len(c))
    if len_cons:
        len_cons_row = np.zeros(g.shape[1], dtype=np.int8)
        len_cons_row[0:m_edges] = -1
        len_cons_row[m_edges:2*m_edges] = -1
        g = sparse.vstack((g, len_cons_row))
        h = np.append(h, np.rint(-average_len))
    G = spmatrix(g.data.tolist(), g.row, g.col, g.shape,  tc='d')
    h = matrix(h)
    c = matrix(c)
    cons = spmatrix(cons.data.tolist(), cons.row, cons.col, cons.shape, tc='d')
    b = matrix(b)

    start = time.time()
    sol = solvers.lp(c, G, h, cons, b, solver='glpk')
    elapsed = time.time() - start
    print 'LP time %f mins.' % (elapsed/60)

    args = np.array(sol['x'])
    args1 = args[np.where(args >=1e-5)]
    args2 = args1[np.where(args1 <= 0.99999)]
    nonint = args2.shape[0]
    args = np.rint(sol['x'])
    norm = sol['primal objective']
    x = args[0:m_edges] - args[m_edges:2*m_edges]
    q = np.zeros((sub_cons_count, m_edges), dtype=int)
    r = np.zeros((sub_cons_count, n_simplices), dtype=int)
    qi_start = 2*m_edges
    for i in range(sub_cons_count):
        qi_end = qi_start + 2*m_edges
        q[i] = (args[qi_start: qi_start+m_edges] - args[qi_start+m_edges: qi_end]).reshape(m_edges,)
        ri_start = qi_end
        ri_end = ri_start + 2*n_simplices
        r[i] = (args[ri_start: ri_start+n_simplices] - args[ri_start+n_simplices: ri_end]).reshape(n_simplices, )
        qi_start = ri_end
    return x, q, r, norm, nonint

def get_lp_inputs(mesh, k_currents, opt='default', w=[], v=[], b_matrix=[], cons=[]):
    m_edges = mesh.edges.shape[0]
    n_simplices = mesh.simplices.shape[0]
    if w == []:
        w = simpvol(mesh.points, mesh.simplices)
    if v == []:
        v = simpvol(mesh.points, mesh.simplices)
    if b_matrix == []:
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges, format='coo')
    if cons == []:
        sub_cons_count = k_currents
        # Msfn option
        if opt == 'msfn':
            sub_cons_count += 1
        m_edges_identity = sparse.identity(m_edges, dtype=np.int8, format='coo')
        identity_cons = sparse.hstack((m_edges_identity, -m_edges_identity))
        sub_cons = sparse.hstack((-m_edges_identity, m_edges_identity, -b_matrix, b_matrix))
        sub_cons_col_count = 2*m_edges + 2*n_simplices
        zero_sub_cons = sparse.coo_matrix((m_edges, sub_cons_col_count), dtype=np.int8)
        for i in range(sub_cons_count):
            if i == 0:
                k_identity_cons = identity_cons
            else:
                k_identity_cons = sparse.vstack((k_identity_cons, identity_cons))
        
        for i in range(0,sub_cons_count):
            for j in range(sub_cons_count):
                if j == i and i == 0:
                    cons_row = sub_cons
                elif j == 0:
                    cons_row = zero_sub_cons 
                elif j == i:
                    cons_row = sparse.hstack((cons_row, sub_cons))
                else:
                    cons_row = sparse.hstack((cons_row, zero_sub_cons))
            if i == 0:
                cons = cons_row
            else:
                cons = sparse.vstack((cons, cons_row))
        cons = sparse.hstack((k_identity_cons, cons))
    # Uncomment the line below to print sub_cons, cons and c
    #print_cons(sub_cons.toarray(), cons.toarray(), c)
    return w, v, b_matrix, cons
