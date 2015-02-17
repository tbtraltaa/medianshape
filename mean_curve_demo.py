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

from utils import sparse_savetxt
from mesh.utils import boundary_matrix, simpvol

from cvxopt import matrix, solvers

options = ['default', 'mass', 'msfn']
#options = ['default']

def load(dirname):
    mesh = Mesh()
    mesh.points = np.loadtxt("%s/points.txt"%dirname) 
    mesh.simplices = np.loadtxt("%s/points.txt"%dirname)
    mesh.edges = np.loadtxt("%s/edges.txt"%dirname)
    v = np.loadtxt("%s/v.txt"%dirname)
    w = np.loadtxt("%s/w.txt"%dirname)
    b_matrix = sparse.dok_matrix((len(mesh.edges), len(mesh.simplices)), dtype=np.int8)
    with open("%s/b_matrix.txt"%dirname, 'r') as f:
        for line in f.readLines():
            data = line.split()
            b_matrix[int(data[0]), int(data[1])] = data[2]
    
    input_currents = []
    for i in range(1,4):
        curve = np.zero((len(mesh.edges), 1), dtype=np.int8)
        with open("%s/curve%d.txt"%(dirname,i), 'r') as f:
            for line in f.readLines():
                data = line.split()
                curve[int(data[0]), int(data[1])] = data[2]
            if input_currents:
                input_currents = np.vstack(input_currents, curve)
            else:
                input_currents = curve
    return mesh, input_currents, b_matrix, w, v
if __name__ == '__main__':
    #start = time.time()
    lp_times = list()
    start = time.time()
    mesh = Mesh()
    #mesh, input_currents, b_matrix, w,v =  load('/home/altaa/dumps1')
    # l - initial length of triangle sides 
    # change it to 1 for big traingles
    mesh.boundary_box = (0,0,200,50)
    mesh.fixed_points = [(0,0),(200,0),(0,50),(200,50)]
    mesh.points, mesh.simplices = distmesh2d('square', mesh.boundary_box, mesh.fixed_points, l=5)
    mesh.set_edges()
    mesh.orient_simplices_2D()
    np.savetxt('/home/altaa/dumps2/points.txt', mesh.points, delimiter=' ')
    np.savetxt('/home/altaa/dumps2/edges.txt', mesh.edges, fmt='%d', delimiter=' ')
    np.savetxt('/home/altaa/dumps2/simplices.txt', mesh.simplices, fmt='%d', delimiter=' ')

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
#                            [7,8]])
    print mesh.get_info()
    #mesh.print_detail()
    pdf_file = PdfPages('/home/altaa/figures2.pdf')
    #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
    #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
    function_sets = [['curve1', 'curve2', 'curve3']]
    figcount = 1
    for j, functions in enumerate(function_sets):
        combinations = utils.get_combination(len(functions))
        combinations = np.array([[1,1,1]])
        fig = plt.figure(figsize=(19,8))
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

        k_currents = len(functions)
        currents = np.array(input_currents).reshape(k_currents, mesh.edges.shape[0])
        w = simpvol(mesh.points, mesh.edges)
        v = simpvol(mesh.points, mesh.simplices)
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)
        np.savetxt('/home/altaa/dumps2/w.txt', w, delimiter=' ')
        np.savetxt('/home/altaa/dumps2/v.txt', v, delimiter=' ')
        sparse_savetxt('/home/altaa/dumps2/b_matrix.txt', b_matrix)
        #np.savetxt('/home/altaa/dumps1/b_matrix.txt', b_matrix, fmt='%d', delimiter=' ')

#            input_currents = list()
#            current1 = np.zeros(shape=(len(mesh.edges),1))
#            current1[0] = 1
#            current1[3] = 1
#            current1[5] = 1 
#            current1[7] = 1
#            current2 = np.zeros(shape=(len(mesh.edges),1))
#            current2[1] = 1
#            current2[13] = -1
#            current2[11] = -1
#            current2[9] = -1
#            input_currents.append(current1)
#            input_currents.append(current2)
#            input_currents = np.array(input_currents)
#            k_currents = 2
        
        nonints = list()
        norms = list()
        t_lens = list()
        for opt in options:
            min_norm = 1000
            average_len = np.average(np.array([c.nonzero()[0].shape[0] for c in input_currents]))
            w, v, b_matrix, cons = mean.get_lp_inputs(mesh,  k_currents, opt, w, v, b_matrix)
            #np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
            lambdas = [0.0001, 0.001, 1]
            for l in lambdas:
                #for comb in combinations[:-1,:]:
                for comb in combinations:
                    input_currents = currents*comb.reshape(comb.size,1) 
                    t, q, r, norm, nonint = mean.mean(mesh, input_currents, l, opt, w, v, cons)
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
                    title = '%s, lambda=%.04f, %s, T_len=%d, T_i_len_ave=%d' % \
                    (opt, l, str(comb), t_len, average_len)
                    figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount, opt, l)
                    plotting.plot_mean(mesh, functions, input_currents, comb, t, title, figname, pdf_file)
                    figcount += 1

                    figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount, opt, l)
                    plotting.plot_curve_and_mean(mesh, functions, input_currents, comb, t, title, \
                    figname, pdf_file)
                    figcount += input_currents.shape[0]

                    figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount,opt,l)
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

    print "non ints", nonints
    print "Norms", norms
    print "t_lens", t_lens
    print "Average len", average_len
    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
