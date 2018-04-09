#encoding: utf-8

'''
Median shape demo 3D
--------------------
'''

from __future__ import absolute_import

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from matplotlib.backends.backend_pdf import PdfPages

from medianshape.simplicial import currentgen 
from medianshape.simplicial.mesh_to_curve import smoothen
from medianshape.utils import simplicial_curve_to_points
from medianshape.experiment.median import runmedians as run, cases3d

from medianshape.viz import plot3d 

def show_median3d():
    '''
    Shows the results of a saved experiment case.
    '''
    mesh, input_currents, t, q, r = cases3d.get_saved_case(dirname='data/curves_on_sphere')
    fig = plt.figure(figsize=(8,8))
    plot3d.plotmesh3d(mesh, mesh.get_info())
    fig.tight_layout()
    plt.show()
    plt.figure(figsize=(8,8))
    title = r"$\lambda=0.0010$, $\mu=0.000010$"
    plot3d.plot_median3d(mesh, input_currents, t, title=title) 
    plt.show()
    plot3d.plot_decomposition3d(mesh, input_currents, t, q, r, title=title)


def mediandemo3d(outdir='data', save=True):
    '''
    Median shape demo in 3D. The experiment case is chosen from 'medianshape.cases3d'.
    Given the experiment case, it gets input currents from the input curves by pushing
    the underlying simplicial complex. Then input the simplicial setting to Median LP which
    solves median current. The experiment result is saved in a given output directory(data).
    '''
    lp_times = list()
    start = time.time()
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/mediandemo3d.pdf'%outdir)

    fig = plt.figure(figsize=(8,8))
    figcount = 1
    mesh, simplices, subsimplices, points, lambdas, mus, is_closed \
    = cases3d.equally_spaced_longitudes3d()
    #= cases3d.tunnel_loops_on_torus_surface3d()
    #= cases3d.handle_loops_on_torus_surface3d()
    #= cases3d.torus_surface3d()
    #= cases3d.equally_spaced_longitudes3d()
    print mesh.get_info()
    if mesh.surf_points is None:
        title = mesh.get_info() 
        title = None
        plot3d.plotmesh3d(mesh, title, '%s/figures/%d'%(outdir, figcount), pdf_file, save)
    else:
        title = None
        plot3d.plot_simplices3d(mesh, np.ones(len(mesh.triangles)), title, '%s/figures/%d'%(outdir, figcount), pdf_file, save)
    figcount += 1
    fig.tight_layout()
    fig = plt.figure(figsize=(8,8))
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh.points, mesh.edges, points, is_closed, mesh.surf_points)

    figname = '%s/figures/%d'%(outdir, figcount)
    title = r'$Curve$ $approximation$'
    plot3d.plot_curves_approx3d(mesh, points, vertices, paths, figname, pdf_file, save, title=None)
    figcount += 1
    fig.tight_layout()
    t_list, q_list, r_list = run.runmedians3d(mesh, simplices, subsimplices, input_currents, lambdas, mus, file_doc=pdf_file, save=save, figcount=figcount)

    # To visualize the medianshape decompsition
    title = None
    for i, t in enumerate(t_list):
        t_points = simplicial_curve_to_points(mesh.points, mesh.edges, t)
        smooth_t = smoothen(t_points)
        smooth_t_tmp = smooth_t.copy()
        start_point = mesh.points[vertices[0][0],:].reshape(-1,3)
        end_point = mesh.points[vertices[0][-1],:].reshape(-1,3)
        dist1 = cdist(t_points[0,:].reshape(-1,3), start_point )[0][0]
        dist2 = cdist(t_points[-1,:].reshape(-1,3), start_point )[0][0]
        if dist1 < dist2:
            smooth_t = np.append(smooth_t.reshape(-1,3), end_point, axis=0)
            smooth_t = np.insert(smooth_t.reshape(-1,3), 0, start_point, axis=0)
        else:
            smooth_t = np.insert(smooth_t, 0, end_point, axis=0)
            smooth_t = np.append(smooth_t, start_point, axis=0)

        figname = '%s/figures/%d'%(outdir, figcount)
        fig = plt.figure(figsize=(8,8))
        plot3d.plot_median3d(mesh, input_currents, t, title, figname, pdf_file, save=save,linewidth=3)
        plt.tight_layout()
        
        figcount += 1
        figname = '%s/figures/%d'%(outdir, figcount)
        fig = plt.figure(figsize=(8,8))
        plot3d.plot_median3d(mesh, input_currents, t, title, figname, pdf_file, save=save,linewidth=3)
        ax = plt.gca(projection='3d')
        ax.plot(smooth_t[:,0],smooth_t[:,1], smooth_t[:,2], 'k', linewidth = 3.0)
        plt.savefig('%s.png'%figname, dpi=300)
        
        figcount += 1
        figname = '%s/figures/%d'%(outdir, figcount)
        fig = plt.figure(figsize=(8,8))
        plot3d.plot_median3d(mesh, input_currents, [], title, figname, pdf_file, save=save,linewidth=3)
        ax = plt.gca(projection='3d')
        ax.plot(smooth_t[:,0],smooth_t[:,1], smooth_t[:,2], 'k', linewidth = 3.0)
        plt.savefig('%s.png'%figname, dpi=300)
        
        figname = '%s/figures/%d'%(outdir, figcount)
        plot3d.plot_decomposition3d(mesh, input_currents, t, q_list[i], r_list[i], title, \
        figname, pdf_file, save, linewidth=3)
        figcount += input_currents.shape[0]
       

        mesh, simplices, subsimplices, input_points, lambdas, mus, is_closed \
        = cases3d.equally_spaced_longitudes3d(n=smooth_t.shape[0])
        
        start_point = input_points[0][0,:].reshape(-1,3)
        end_point = input_points[0][-1,:].reshape(-1,3)
        if dist1 < dist2:
            smooth_t_tmp = np.append(smooth_t_tmp.reshape(-1,3), end_point, axis=0)
            smooth_t_tmp = np.insert(smooth_t_tmp.reshape(-1,3), 0, start_point, axis=0)
        else:
            smooth_t_tmp = np.insert(smooth_t_tmp, 0, end_point, axis=0)
            smooth_t_tmp = np.append(smooth_t_tmp, start_point, axis=0)
        
        figname = '%s/figures/%d'%(outdir, figcount)
        fig = plt.figure(figsize=(8,8))
        plt.gca().set_aspect('equal')
        ax = plt.gca(projection='3d')
        colors = plot3d.get_colors(len(input_currents))
        for i, input_point in enumerate(input_points):
            ax.plot(input_point[:,0],input_point[:,1], input_point[:,2], color=colors[i], linewidth = 3.0)
        ax.plot(smooth_t_tmp[:,0],smooth_t_tmp[:,1], smooth_t_tmp[:,2], 'k', linewidth = 3.0)
        plt.savefig('%s.png'%figname, dpi=300)
        pdf_file.savefig(fig)
    if save:
        pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    plt.show()

if __name__ == '__main__':
    mediandemo3d(save=True)
    #show_median3d()
