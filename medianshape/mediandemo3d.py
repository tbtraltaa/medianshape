#encoding: utf-8

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
import median

from utils import sparse_savetxt, load, save
from mesh.utils import boundary_matrix, simpvol, get_subsimplices

from cvxopt import matrix, solvers
import distmesh as dm
import plot3d

def equators():
    # l - initial length of triangle sides. Change it to vary traingle size
    boundary_box = [0,0,0,20,20,20]
    l = 4
    mesh, w, v, b_matrix = load_mesh(boundary_box, l, include_corners=False)
    points = pointgen3d.sphere_equator(mesh.bbox, 4)
    simplices = mesh.edges
    subsimplices = np.arange(mesh.points.shape[0]).reshape(-1,1)
    w = simpvol(mesh.points, subsimplices)
    v = simpvol(mesh.points, simplices)
    b_matrix = boundary_matrix(simplices, subsimplices)
    dim = 0
    return mesh, simplices, subsimplices, w, v, b_matrix, points.reshape(-1, 1, points.shape[1]), l, dim

def equally_spaced_longitudes(): 
    # l - initial length of triangle sides. Change it to vary traingle size
    boundary_box = [0,0,0,20,20,20]
    l = 4
    mesh, w, v, b_matrix = load_mesh(boundary_box, l, include_corners=False)
    curve1 = pointgen3d.sphere_arc(mesh.bbox, 0, 10)
    curve2 = pointgen3d.sphere_arc(mesh.bbox, 2*np.pi/3, 10)
    curve3 = pointgen3d.sphere_arc(mesh.bbox, 4*np.pi/3, 10)
    shapes = [curve1, curve2, curve3]
    points  = np.array(shapes)
    return mesh, mesh.triangles, mesh.edges, w, v, b_matrix, points, l, 1

def differently_spaced_longitudes(): 
    # l - initial length of triangle sides. Change it to vary traingle size
    boundary_box = [0,0,0,20,20,20]
    l = 4
    mesh, w, v, b_matrix = load_mesh(boundary_box, l, include_corners=False)
    curve1 = pointgen3d.sphere_arc(mesh.bbox, 0, 10)
    curve2 = pointgen3d.sphere_arc(mesh.bbox, np.pi/4, 10)
    curve3 = pointgen3d.sphere_arc(mesh.bbox, 9*np.pi/8, 10)
    shapes = [curve1, curve2, curve3]
    points  = np.array(shapes)
    return mesh, mesh.triangles, mesh.edges, w, v, b_matrix, points, l, 1

def load_mesh(boundary_box=None, l=0.2, fixed_points=None, include_corners=True, load_data=False):
    if load_data:
        #TODO add mesh diagonal
        mesh, w, v, b_matrix = load()
    else:
        mesh = Mesh3D()
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

def run_demo(mesh, simplices, subsimplices, input_currents, lambdas, mus, w=None, v=None, b_matrix=None, file_doc=None, save_data=True, dim=1):
    figcount = 2
    norms = list()
    t_lens = list()
    if w is None:
        w = simpvol(mesh.points, subsimplices)
    if v is None:
        v = simpvol(mesh.points, simplices)
    if b_matrix is None:
        b_matrix = boundary_matrix(simplices, subsimplices)
    k_currents = len(input_currents)
    w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, simplices, subsimplices,  k_currents, w, v, b_matrix)
    #np.savetxt('output/dumps/cons-%s.txt'%cons, fmt='%d', delimiter=' ')
    for l in lambdas:
        for mu in mus:
            t, q, r, norm = median.median(mesh.points, simplices, subsimplices, input_currents, l, w, v, cons, mu=mu)
            if save_data:
                save(t=t, lambda_=l, mu=mu)
            norms.append(norm)
            t_len = len(t.nonzero()[0])
            t_lens.append(t_len)
            title = 'MRSMS, lambda=%.04f, mu=%.06f'%(l, mu)
            figname = 'output/figures/%d-%.04f-%.06f'%(figcount, l, mu)
            if save_data and file_doc is not None:
                fig = plt.figure(figsize=(8,8))
                plot3d.plot_median(mesh, input_currents, t, title, figname, file_doc, save=save, dim=dim)
                plt.tight_layout()
                plt.show()
                fig = plt.figure(figsize=(8,8))
                figcount += 1

                #figname = 'output/figures/%d-%s-%.04f-%.04f'%(figcount, l, mu)
                #plotting.plot_curve_and_median(mesh, input_currents, comb, t, title, \
                #figname, file_doc, save)
                #figcount += input_currents.shape[0]

                figname = 'output/figures/%d-%.06f-%.06f'%(figcount, l, mu)
                plot3d.plot_decomposition(mesh, input_currents, t, q, r, title, \
                figname, file_doc, save, dim=dim)
                figcount += input_currents.shape[0]
    return t

def mediandemo3d(load_data=False, save_data=True):
    lp_times = list()
    start = time.time()
    pdf_file = PdfPages('output/figures-l=2.5.pdf')

    fig = plt.figure(figsize=(8,8))
    figcount = 1
    mesh, simplices, subsimplices, w, v, b_matrix, points, l, dim = equally_spaced_longitudes()
    print mesh.get_info()
    plot3d.plotmesh3d(mesh, mesh.get_info())
    pdf_file.savefig(fig, pad_inches=-1, box_inches='tight')
    plt.savefig("output/figures/mesh.png")
    fig.tight_layout()
    plt.show()
    fig = plt.figure(figsize=(8,8))
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, simplices, subsimplices, points, is_closed=False)

    figname = 'output/figures/%d.png'%(figcount)
    title = '%s-l=%0.2f' % (mesh.get_info(),l)
    title = 'Curve approximation'
    plot3d.plot_curves_approx(mesh, points, vertices, paths, title, figname, pdf_file)
    figcount += 1
    fig.tight_layout()
    plt.show()
    if save_data:
        save(mesh, input_currents, b_matrix, w, v)
    lambdas = [0.001]
    mus = [0.00001]
    t = run_demo(mesh, simplices, subsimplices, input_currents, lambdas, mus, w, v, b_matrix, pdf_file, dim=dim)
    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mediandemo3d()
