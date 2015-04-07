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
from mesh.mesh import Mesh2D
from shape_gen import point_gen, curve_gen, utils
import mean
import plotting

from utils import sparse_savetxt, load, save, envelope, adjust_alphas
from mesh.utils import boundary_matrix, simpvol, get_subsimplices

from cvxopt import matrix, solvers
import distmesh as dm

#options = ['default', 'mass', 'msfn']
options = ['mass']

def load_mesh(boundary_box=None, fixed_points=None, l=0.02, load_data=False):
    if load_data:
        #TODO add mesh diagonal
        mesh, w, v, b_matrix = load()
    else:
        mesh = Mesh2D()
        # l - initial length of triangle sides. Change it to vary traingle size
        mesh.boundary_box = boundary_box
        mesh.fixed_points = fixed_points
        mesh.diagonal = np.sqrt(mesh.boundary_box[2]**2 + mesh.boundary_box[3]**2)
        mesh.points, mesh.simplices = distmesh2d('square', mesh.boundary_box, mesh.fixed_points, l=l)
        mesh.edges = get_subsimplices(mesh.simplices)
        #mesh.edges, temp = dm.mkt2t(mesh.simplices)
        #print mesh.edges
        #print temp
        mesh.orient_simplices_2D()
#        mesh.points = np.array([[0,0], [0,0.5],[0,1],[0.5,1],[1,1],[1,0.5],[1,0],[0.5,0], [0.5, 0.5]])
#        mesh.simplices = np.array([[0,8,1],
#                                    [1,8,2],
#                                    [2,8,3],
#                                    [3,8,4],
#                                    [4,8,5],
#                                    [5,8,6],
#                                    [6,8,7],
#                                    [7,8,0]])
#
#        mesh.edges = np.array([[0,1],
#                                [0,7],
#                                [0,8],
#                                [1,2],
#                                [1,8],
#                                [2,3],
#                                [2,8],
#                                [3,4],
#                                [3,8],
#                                [4,5],
#                                [4,8],
#                                [5,6],
#                                [5,8],
#                                [6,7],
#                                [6,8],
#                                [7,8]])
        w = simpvol(mesh.points, mesh.edges)
        v = simpvol(mesh.points, mesh.simplices)
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)
    return mesh, w, v, b_matrix

def run_demo(mesh, input_currents, options, lambdas, mus, alphas, w=None, v=None, b_matrix=None, file_doc=None, save=True):
    figcount = 2
    norms = list()
    t_lens = list()
    if w is None:
        w = simpvol(mesh.points, mesh.edges)
    if v is None:
        v = simpvol(mesh.points, mesh.simplices)
    if b_matrix is None:
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)
    k_currents = len(input_currents)
    for opt in options:
        average_len = np.average(np.array([c.nonzero()[0].shape[0] for c in input_currents]))
        w, v, b_matrix, cons = mean.get_lp_inputs(mesh,  k_currents, opt, w, v, b_matrix)
        #np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
        for l in lambdas:
            comb=[1,1,1]
            #for comb in combinations[:-1,:]:
            #for comb in combinations:
                #input_currents = currents*comb.reshape(comb.size,1) 
            for mu in mus:
                for alpha in alphas:
                    t, q, r, norm = mean.mean(mesh, input_currents, l, opt, w, v, cons, mu=mu, alpha=alpha)
                    if save:
                        #save(t=t, opt=opt, lambda_=l)
                        pass
                    norms.append(norm)
                    t_len = len(t.nonzero()[0])
                    t_lens.append(t_len)
                    title = '%s, lambda=%.04f, mu=%.06f, alpha=%s' % \
                    (opt, l, mu, str(alpha))
                    figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.06f'%(figcount, opt, l, mu)
                    if save and file_doc is not None:
                        plotting.plot_mean(mesh, input_currents, comb, t, title, figname, file_doc, save=save)
                        figcount += 1

                        #figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.04f'%(figcount, opt, l, mu)
                        #plotting.plot_curve_and_mean(mesh, input_currents, comb, t, title, \
                        #figname, file_doc, save)
                        #figcount += input_currents.shape[0]

                        #figname = '/home/altaa/fig_dump/%d-%s-%.06f-%.06f'%(figcount,opt,l, mu)
                        #plotting.plot_decomposition(mesh, input_currents, comb, t, q, r, title, \
                        #figname, file_doc, save)
                        #figcount += input_currents.shape[0]
                
                # Plotting the combination with minimum flatnorm difference
            #title = 'Minimum flatnorm difference, %s, lambda=%.04f, %s' % (opt, l, str(comb))
#                figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount, opt, l)
#                plotting.plot_mean(mesh, min_currents, min_comb, min_t, title, \
#                figname, file_doc)
#                figcount += 1
#
#                figname = '/home/altaa/fig_dump/%d-%s-%.04f'%(figcount,opt,l)
#                plotting.plot_decomposition(mesh, min_currents, comb, min_t, min_q, min_r, \
#                title, figname, file_doc)
#                figcount += input_currents.shape[0]
    print "Norms", norms
    print "t_lens", t_lens
    print "Average len", average_len
    return t

def mean_curve_demo(load_data=False, save_data=True):
    lp_times = list()
    start = time.time()
    pdf_file = PdfPages('/home/altaa/figures.pdf')
    fig = plt.figure(figsize=(16,8))
    #fig = plt.figure()
    figcount = 1
    boundary_box = (0,0,200,50)
    fixed_points = [(0,0),(200,0),(0,50),(200,50)]
    l=6
    boundary_box = (0,0,40,40)
    fixed_points = [(0,0),(40,0),(0,40),(40,40)]
    boundary_box = (0,0,1,1)
    fixed_points = [(0,0),(1,0),(0,1),(1,1)]
    l=0.07
    mesh, w, v, b_matrix = load_mesh(boundary_box, fixed_points, l)
    print mesh.get_info()

    #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
    #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
    #functions= ['curve4', 'curve5']
    functions= ['curve1', 'curve2']
    combinations = utils.get_combination(len(functions))
    combinations = np.array([[1,1,1]])
    #combinations = np.array([[1,1,1]])
    ellipse1 = point_gen.sample_ellipse(0.4, 0.2, 10)
    ellipse2 = point_gen.sample_ellipse(0.2, 0.4, 10)
    shapes = [ellipse1, ellipse2]
    #curve1 = point_gen.sample_curve1()
    #curve2 = point_gen.sample_curve2()
    #mesh.plot()
    #plt.scatter(curve1[:,0], curve1[:,1], color='r')
    #plt.scatter(curve2[:,0], curve2[:,1], color='k')
    #plt.show()
    #shapes = [curve1, curve2]
    points = list()
    for f in functions:
        points.append(point_gen.sample_function_mesh(mesh, f))
    points  = np.array(shapes)
    vertices, paths, input_currents = curve_gen.push_curves_on_mesh(mesh, points, is_closed=True)

    figname = '/home/altaa/fig_dump/%d.png'%(figcount)
    title = '%s - (%s)' % (mesh.get_info(), ','.join(functions))
    plotting.plot_curves_approx(mesh, points, vertices, paths, title, figname, pdf_file)
    figcount += 1
    #envelope(mesh, input_currents)
    lambdas = [1]
    mus = [0.0001]
    alpha1 = np.array([0])
    #alpha1 = np.append(alpha1, np.linspace(0.4999, 0.5, 10))
    #alpha1 = np.append(alpha1, np.linspace(0.5, 0.5001, 10))
    alpha1 = np.append(alpha1, np.linspace(0.4, 0.5, 10))
    alpha1 = np.append(alpha1, np.linspace(0.5, 0.6, 10))
    alpha1 = np.append(alpha1, np.array([1]))
    alpha1 = alpha1.reshape(alpha1.size, 1) 
    alpha2 = (1-alpha1).reshape(alpha1.size, 1)
    alphas = np.hstack((alpha1, alpha2))
    alphas = np.ndarray(shape=(1,2), buffer=np.array([0.5, 0.5]))

    t = run_demo(mesh, input_currents, options, lambdas, mus, alphas, w, v, b_matrix, pdf_file)
    alphas = alphas + adjust_alphas(mesh, input_currents, t, v)
    t = run_demo(mesh, input_currents, options, lambdas, mus, alphas, w, v, b_matrix, pdf_file)


    if save_data:
        save(mesh, input_currents, b_matrix, w, v)

    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mean_curve_demo()
