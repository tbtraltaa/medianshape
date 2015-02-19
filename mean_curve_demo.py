# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import time


import numpy as np
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mesh.distmesh import distmesh2d
from mesh.mesh import Mesh
from shape_gen import point_gen, curve_gen, utils
import mean
import plotting

from utils import sparse_savetxt, load, save
from mesh.utils import boundary_matrix, simpvol

from cvxopt import matrix, solvers

#options = ['default', 'mass', 'msfn']
options = ['mass']

def mean_curve_demo(load_data=False, save_data=True):
    #start = time.time()
    lp_times = list()
    start = time.time()
    pdf_file = PdfPages('/home/altaa/figures.pdf')
    fig = plt.figure(figsize=(19,8))
    figcount = 1
    if load_data:
        mesh, input_currents, b_matrix, w, v = load()
    else:
        mesh = Mesh()
        # l - initial length of triangle sides. Change it to vary traingle size
        mesh.boundary_box = (0,0,200,50)
        mesh.fixed_points = [(0,0),(200,0),(0,50),(200,50)]
        mesh.points, mesh.simplices = distmesh2d('square', mesh.boundary_box, mesh.fixed_points, l=7)
        mesh.set_edges()
        mesh.orient_simplices_2D()

#    mesh.points = np.array([[0,0], [0,0.5],[0,1],[0.5,1],[1,1],[1,0.5],[1,0],[0.5,0], [0.5, 0.5]])
#    mesh.simplices = np.array([[0,8,1],
#                                [1,8,2],
#                                [2,8,3],
#                                [3,8,4],
#                                [4,8,5],
#                                [5,8,6],
#                                [6,8,7],
#                                [7,8,0]])
#
#    mesh.edges = np.array([[0,1],
#                            [0,7],
#                            [0,8],
#                            [1,2],
#                            [1,8],
#                            [2,3],
#                            [2,8],
#                            [3,4],
#                            [3,8],
#                            [4,5],
#                            [4,8],
#                            [5,6],
#                            [5,8],
#                            [6,7],
#                            [6,8],
        #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
        #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
        functions= ['curve4', 'curve5']
        combinations = utils.get_combination(len(functions))
        combinations = np.array([[1,1,1]])
        #combinations = np.array([[1,1,1]])
        #w = np.ndarray(shape=(len(mesh.edges),))
        #w[0:] = 1
        #v =  np.ndarray(shape=(len(mesh.simplices),))
        #v[0:] = 0.433
        points, vertices, paths, input_currents = curve_gen.push_curves_on_mesh(mesh, functions)

        figname = '/home/altaa/fig_dump/%d.png'%(figcount)
        title = 'Functions - %s - (%s)' % (mesh.get_info(), ','.join(functions))
        plotting.plot_curves_approx(mesh, points, vertices, paths, title, figname, pdf_file)
        figcount += 1

        w = simpvol(mesh.points, mesh.edges)
        v = simpvol(mesh.points, mesh.simplices)
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)

    print mesh.get_info()
    #mesh.print_detail()
    k_currents = len(input_currents)
    input_currents = np.array(input_currents).reshape(k_currents, mesh.edges.shape[0])
    nonints = list()
    norms = list()
    t_lens = list()
    for opt in options:
        min_norm = 1000
        average_len = np.average(np.array([c.nonzero()[0].shape[0] for c in input_currents]))
        w, v, b_matrix, cons = mean.get_lp_inputs(mesh,  k_currents, opt, w, v, b_matrix)
        #np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
        lambdas = [0.0001]
        mus = [0.0001, 1]
        for l in lambdas:
            #for comb in combinations[:-1,:]:
            #for comb in combinations:
                #input_currents = currents*comb.reshape(comb.size,1) 
            for mu in mus:
                t, q, r, norm, nonint = mean.mean(mesh, input_currents, l, opt, w, v, cons, mu=mu)
                if save_data:
                    save(t=t, opt=opt, lambda_=l)
                nonints.append(nonint)
                norms.append(norm)
                t_len = len(t.nonzero()[0])
                t_lens.append(t_len)
                if norm < min_norm:
                    min_norm = norm
                    min_comb = comb
                    min_t = t
                    min_q = q
                    min_r = r
                    min_currents = input_currents

                title = '%s, lambda=%.04f, mu=%.04f, %s, T_len=%d, T_i_len_ave=%d' % \
                (opt, l, mu, str(comb), t_len, average_len)
                figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.04f'%(figcount, opt, l, mu)
                plotting.plot_mean(mesh, functions, input_currents, comb, t, title, figname, pdf_file)
                figcount += 1

                #figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.04f'%(figcount, opt, l, mu)
                #plotting.plot_curve_and_mean(mesh, functions, input_currents, comb, t, title, \
                #figname, pdf_file)
                #figcount += input_currents.shape[0]

                figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.04f'%(figcount,opt,l, mu)
                plotting.plot_decomposition(mesh, functions, input_currents, comb, t, q, r, title, \
                figname, pdf_file)
                figcount += input_currents.shape[0]
                
                # Plotting the combination with minimum flatnorm difference
            #title = 'Minimum flatnorm difference, %s, lambda=%.04f, %s' % (opt, l, str(comb))
#                figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount, opt, l)
#                plotting.plot_mean(mesh, functions, min_currents, min_comb, min_t, title, \
#                figname, pdf_file)
#                figcount += 1
#
#                figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount,opt,l)
#                plotting.plot_decomposition(mesh, functions, min_currents, comb, min_t, min_q, min_r, \
#                title, figname, pdf_file)
#                figcount += input_currents.shape[0]
    if save_data:
        save(mesh, input_currents, b_matrix, w, v)

    print "non ints", nonints
    print "Norms", norms
    print "t_lens", t_lens
    print "Average len", average_len
    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mean_curve_demo()
