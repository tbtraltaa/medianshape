# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import time


import numpy as np
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from mesh.meshgen import distmesh3d, scipy_mesh3d
from mesh.mesh import Mesh3D
from shapegen import pointgen3d, currentgen, utils
import mean

from utils import sparse_savetxt, load, save, envelope, adjust_alphas
from mesh.utils import boundary_matrix, simpvol, get_subsimplices

from cvxopt import matrix, solvers
import distmesh as dm
import plot3d

#options = ['default', 'mass', 'msfn']
options = ['mass']

def load_mesh(boundary_box=None, l=0.2, fixed_points=None, include_corners=True, load_data=False):
    if load_data:
        #TODO add mesh diagonal
        mesh, w, v, b_matrix = load()
    else:
        mesh = Mesh3D()
        # l - initial length of triangle sides. Change it to vary traingle size
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
        #mesh.fixed_points = np.array([[0.5, 0.5, 0.5]])
        mesh.points, mesh.simplices= distmesh3d("sphere", mesh.bbox, mesh.fixed_points, l=l)
        #mesh.points, mesh.simplices = scipy_mesh3d(mesh.bbox, mesh.fixed_points, l)
        #mesh.points, mesh.simplices= meshpy_cuboid3d(mesh.bbox, mesh.fixed_points.tolist(), max_volume=(1**3)*1.0/(6*np.sqrt(2)))
        #plot3d.plotmesh3d(mesh)
        mesh.triangles = get_subsimplices(mesh.simplices)
        mesh.edges = get_subsimplices(mesh.triangles)
        w = simpvol(mesh.points, mesh.edges)
        v = simpvol(mesh.points, mesh.triangles)
        b_matrix = boundary_matrix(mesh.triangles, mesh.edges)
    return mesh, w, v, b_matrix

def run_demo(mesh, input_currents, options, lambdas, mus, alphas, w=None, v=None, b_matrix=None, file_doc=None, save=True):
    figcount = 2
    norms = list()
    t_lens = list()
    if w is None:
        w = simpvol(mesh.points, mesh.edges)
    if v is None:
        v = simpvol(mesh.points, mesh.triangles)
    if b_matrix is None:
        b_matrix = boundary_matrix(mesh.triangles, mesh.edges)
    k_currents = len(input_currents)
    for opt in options:
        average_len = np.average(np.array([c.nonzero()[0].shape[0] for c in input_currents]))
        w, v, b_matrix, cons = mean.get_lp_inputs(mesh.points, mesh.triangles, mesh.edges,  k_currents, opt, w, v, b_matrix)
        #np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
        for l in lambdas:
            comb=[1,1,1]
            #for comb in combinations[:-1,:]:
            #for comb in combinations:
                #input_currents = currents*comb.reshape(comb.size,1) 
            for mu in mus:
                for alpha in alphas:
                    t, q, r, norm = mean.mean(mesh.points, mesh.triangles, mesh.edges, input_currents, l, opt, w, v, cons, mu=mu, alpha=alpha)
                    if save:
                        save(t=t, opt=opt, lambda_=l)
                    norms.append(norm)
                    t_len = len(t.nonzero()[0])
                    t_lens.append(t_len)
                    title = '%s, lambda=%.04f, mu=%.06f'  % \
                    (opt, l, mu)
                    figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.06f'%(figcount, opt, l, mu)
                    if save and file_doc is not None:
                        plot3d.plot_mean(mesh, input_currents, comb, t, title, figname, file_doc, save=save)
                        #plt.show()
                        figcount += 1

                        #figname = '/home/altaa/fig_dump/%d-%s-%.04f-%.04f'%(figcount, opt, l, mu)
                        #plotting.plot_curve_and_mean(mesh, input_currents, comb, t, title, \
                        #figname, file_doc, save)
                        #figcount += input_currents.shape[0]

                        figname = '/home/altaa/fig_dump/%d-%s-%.06f-%.06f'%(figcount,opt,l, mu)
                        plot3d.plot_decomposition(mesh, input_currents, comb, t, q, r, title, \
                        figname, file_doc, save)
                        figcount += input_currents.shape[0]
                        #plt.show()
                
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
    return t

def mean_curve_demo(load_data=False, save_data=True):
    lp_times = list()
    start = time.time()
    pdf_file = PdfPages('/home/altaa/figures.pdf')

    fig = plt.figure(figsize=(12,8))
    figcount = 1
    boundary_box = (0,0,200,50)
    fixed_points = [(0,0),(200,0),(0,50),(200,50)]
    l=6
    boundary_box = (0,0,40,40)
    fixed_points = [(0,0),(40,0),(0,40),(40,40)]
    boundary_box = (0,0,1,1)
    fixed_points = [(0,0),(1,0),(0,1),(1,1)]
    l=0.07
    boundary_box = [0,0,0,20,20,20]
    l=1
    mesh, w, v, b_matrix = load_mesh(boundary_box, l, include_corners=True)
    print mesh.get_info()

    #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
    #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
    #functions= ['curve4', 'curve5']
    combinations = np.array([[1,1,1]])
    #combinations = np.array([[1,1,1]])
    #curve1 = pointgen3d.curve1(mesh.bbox)
    #curve2 = pointgen3d.curve2(mesh.bbox)
    #shapes = [curve1, curve2]
    curve1 = pointgen3d.sphere_arc(mesh.bbox, 0, 10)
    curve2 = pointgen3d.sphere_arc(mesh.bbox, 2*np.pi/3, 10)
    curve3 = pointgen3d.sphere_arc(mesh.bbox, 4*np.pi/3, 10)
    shapes = [curve1, curve2, curve3]
    points  = np.array(shapes)
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, points, is_closed=False)

    figname = '/home/altaa/fig_dump/%d.png'%(figcount)
    title = '%s-l=%0.2f' % (mesh.get_info(),l)
    plot3d.plot_curves_approx(mesh, points, vertices, paths, title, figname, pdf_file)
    figcount += 1
    #plt.show()
    #envelope(mesh, input_currents)
    if save_data:
        save(mesh, input_currents, b_matrix, w, v)
    lambdas = [0.001]
    mus = [0.00001]
    #alphas = np.ndarray(shape=(1,2), buffer=np.array([0.5, 0.5]))
    alphas = np.ndarray(shape=(1,3), buffer=np.array([1.0/3, 1.0/3, 1.0/3]))
    t = run_demo(mesh, input_currents, options, lambdas, mus, alphas, w, v, b_matrix, pdf_file)
    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mean_curve_demo()
