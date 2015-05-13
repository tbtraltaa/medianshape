# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import time


import numpy as np
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mesh.meshgen import distmesh2d
from mesh.mesh import Mesh2D
from shapegen import pointgen2d, currentgen, utils
import median
import plot2d

from utils import sparse_savetxt, load, save
from mesh.utils import boundary_matrix, boundary_matrix_2d, simpvol, get_subsimplices

from cvxopt import matrix, solvers
import distmesh as dm

options = ['MRSMS']

def triangles():
    mesh = Mesh2D()
    mesh.bbox = [0, 0, 1, 1]
    mesh.set_boundary_points()
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.fixed_points = fixed_points
    mesh.points = np.array([[0,0], [0,0.5], [0,1], [0.5,1], [1,1], \
                            [1,0.5],[1,0],[0.5,0], [0.5, 0.5]])
    mesh.simplices = np.array([[0,8,1],
                                [1,8,2],
                                [2,8,3],
                                [3,8,4],
                                [4,8,5],
                                [5,8,6],
                                [6,8,7],
                                [7,8,0]])

    mesh.edges = np.array([[0,1],
                            [0,7],
                            [0,8],
                            [1,2],
                            [1,8],
                            [2,3],
                            [2,8],
                            [3,4],
                            [3,8],
                            [4,5],
                            [4,8],
                            [5,6],
                            [5,8],
                            [6,7],
                            [6,8],
                            [7,8]])
def ellipses():
    boundary_box = (0,0,1,1)
    l=0.07
    mesh, w, v, b_matrix = load_mesh(boundary_box, l)
    ellipse1 = pointgen2d.sample_ellipse(0.4, 0.2, 10)
    ellipse2 = pointgen2d.sample_ellipse(0.2, 0.4, 10)
    shapes = [ellipse1, ellipse2]
    return mesh, mesh.simplices, mesh.edges, w, v, b_matrix, np.array(shapes),l
    

def load_mesh(boundary_box=None, l=0.02, fixed_points=None, include_corners=True, load_data=False):
    if load_data:
        #TODO add mesh diagonal
        mesh, w, v, b_matrix = load()
    else:
        mesh = Mesh2D()
        #l - initial length of triangle sides. Change it to vary traingle size
        mesh.bbox = boundary_box
        mesh.set_boundary_points()
        mesh.set_diagonal()
        mesh.set_boundary_values()
        mesh.fixed_points = fixed_points
        if include_corners:
            if mesh.fixed_points is not None:
                mesh.fixed_points = np.vstack((mesh.fixed_points, mesh.boundary_points))
            else:
                mesh.fixed_points = mesh.boundary_points
        mesh.points, mesh.simplices = distmesh2d('square', mesh.bbox, mesh.fixed_points, l=l)
        mesh.edges = get_subsimplices(mesh.simplices)
        mesh.orient_simplices_2D()
        w = simpvol(mesh.points, mesh.edges)
        v = simpvol(mesh.points, mesh.simplices)
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)
    return mesh, w, v, b_matrix

def run_demo(mesh, input_currents, options, lambdas, mus, w=None, v=None, b_matrix=None, file_doc=None, save_data=True):
    figcount = 2
    if w is None:
        w = simpvol(mesh.points, mesh.edges)
    if v is None:
        v = simpvol(mesh.points, mesh.simplices)
    if b_matrix is None:
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)
    k_currents = len(input_currents)
    for opt in options:
        w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, mesh.simplices, mesh.edges,  k_currents, opt, w, v, b_matrix)
        for l in lambdas:
            comb=[1,1,1]
            #for comb in combinations[:-1,:]:
            #for comb in combinations:
                #input_currents = currents*comb.reshape(comb.size,1) 
            for mu in mus:
                t, q, r, norm = median.median(mesh.points, mesh.simplices, mesh.edges, \
                input_currents, l, opt, w, v, cons, mu=mu)
                if save_data:
                    save(t=t, opt=opt, lambda_=l)
                title = '%s, lambda=%.04f, mu=%.04f'  % \
                (opt, l, mu)
                figname = 'output/figures/%d-%s-%.04f-%.06f'%(figcount, opt, l, mu)
                if save and file_doc is not None:
                    plot2d.plot_median(mesh, input_currents, comb, t, title, figname, file_doc, save=save)
                    plt.tight_layout()
                    plt.show()
                    fig = plt.figure(figsize=(8,8))
                    figcount += 1

                    #figname = 'output/figures/%d-%s-%.04f-%.04f'%(figcount, opt, l, mu)
                    #plot2d.plot_curve_and_median(mesh, input_currents, comb, t, title, \
                    #figname, file_doc, save)
                    #figcount += input_currents.shape[0]

                    figname = 'output/figures/%d-%s-%.06f-%.06f'%(figcount,opt,l, mu)
                    plot2d.plot_decomposition(mesh, input_currents, comb, t, q, r, title, \
                    figname, file_doc, save)
                    figcount += input_currents.shape[0]
                
    return t

def mediandemo2d(load_data=False, save_data=True):
    lp_times = list()
    start = time.time()
    pdf_file = PdfPages('output/figures.pdf')
    fig = plt.figure(figsize=(8,8))
    #fig = plt.figure()
    figcount = 1
    boundary_box = (0,0,200,50)
    l=6
    #boundary_box = (0,0,40,40)
    #fixed_points = [(0,0),(40,0),(0,40),(40,40)]

    #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
    #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
    #functions= ['curve4', 'curve5']
    functions= ['curve1', 'curve2']
    combinations = utils.get_combination(len(functions))
    combinations = np.array([[1,1,1]])
    #combinations = np.array([[1,1,1]])
    #curve1 = pointgen2d.sample_curve1()
    #curve2 = pointgen2d.sample_curve2()
    #mesh.plot()
    #plt.scatter(curve1[:,0], curve1[:,1], color='r')
    #plt.scatter(curve2[:,0], curve2[:,1], color='k')
    #plt.show()
    #shapes = [curve1, curve2]
    #plot2d.plot_curve(mesh, c)
    #points = list()
    #for f in functions:
     #   points.append(pointgen2d.sample_function_mesh(mesh, f))
    mesh, simplices, subsimplices, w, v, b_matrix, points, l  = ellipses() 
    print mesh.get_info()
    #vertices, paths, input_currents = currentgen.push_functions_on_mesh_2d(mesh, points, is_closed=False, functions=functions)
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, simplices, subsimplices, points, True)
    figname = 'output/figures/%d.png'%(figcount)
    title = mesh.get_info()
    plot2d.plot_curves_approx(mesh, points, vertices, paths, title, figname, pdf_file)
    plt.tight_layout()
    plt.show()
    fig = plt.figure(figsize=(8,8))
    figcount += 1
    lambdas = [0.001]
    lambdas = [1]
    mus = [0.0001]

    t = run_demo(mesh, input_currents, options, lambdas, mus, w, v, b_matrix, pdf_file)
    if save_data:
        save(mesh, input_currents, b_matrix, w, v)

    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mediandemo2d()
