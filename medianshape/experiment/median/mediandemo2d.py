# encoding: utf-8
'''
Median shape demo 2D
--------------------

'''
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from matplotlib.backends.backend_pdf import PdfPages

from medianshape.simplicial import currentgen
from medianshape.simplicial.mesh_to_curve import smoothen
from medianshape.experiment.median import runmedians as run, cases2d

from medianshape.viz import plot2d 
from medianshape.utils import simplicial_curve_to_points

def mediandemo2d(outdir='data', save=True):
    '''
    Median shape demo in 2D. The experiment case is chosen from 'medianshape.cases2d'.
    Given the experiment case, it gets input currents from the input curves by pushing
    the underlying simplicial complex. Then input the simplicial setting to Median LP which
    solves median current. The experiment result is saved in an output directory(data)'.

    '''
    lp_times = list()
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/mediandemo2d.pdf'%outdir)
    fig = plt.figure(figsize=(8,8))
    figcount = 1

    mesh, simplices, subsimplices, points, lambdas, mus, is_closed = cases2d.multicurves2d() 
    
    start = time.time()
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh.points, mesh.edges, points, is_closed=is_closed)
    elapsed = time.time() - start
    print 'Elapsed time %d' % elapsed
    print mesh.points.shape[0]
    figname = '%s/figures/%d'%(outdir, figcount)
    title = mesh.get_info() + r' - $Curve$ $approximation$'
    plot2d.plot_curves_approx2d(mesh, points, vertices, paths, title, figname, pdf_file, save)
    #plt.show()
    figcount += 1
    
    t_list, q_list, r_list = run.runmedians2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, save=save)
    # To visualize the medianshape decompsition
    title = ''
    for i, t in enumerate(t_list):
        t_points = simplicial_curve_to_points(mesh.points, mesh.edges, t)
        smooth_t = smoothen(t_points)
        '''
        start_point = points[0][0,:].reshape(-1,2)
        end_point = points[0][-1,:].reshape(-1,2)
        dist1 = cdist(t_points[0,:].reshape(-1,2), start_point )[0][0]
        dist2 = cdist(t_points[-1,:].reshape(-1,2), start_point )[0][0]
        if dist1 < dist2:
            smooth_t = np.append(smooth_t.reshape(-1,2), end_point, axis=0)
            smooth_t = np.insert(smooth_t.reshape(-1,2), 0, start_point, axis=0)
        else:
            smooth_t = np.insert(smooth_t, 0, end_point, axis=0)
            smooth_t = np.append(smooth_t, start_point, axis=0)
        '''
        figname = '%s/figures/%d'%(outdir, figcount)
        plot2d.plot_median2d(mesh, input_currents, t, title, figname, pdf_file, save=save,linewidth=3)
        
        figcount += 1
        figname = '%s/figures/%d'%(outdir, figcount)
        plot2d.plot_median2d(mesh, input_currents, t, title, figname, pdf_file, save=save,linewidth=3)
        plt.plot(smooth_t[:,0],smooth_t[:,1], 'k', linewidth = 3.0)
        plt.savefig('%s.png'%figname, dpi=300)
        
        figcount += 1
        figname = '%s/figures/%d'%(outdir, figcount)
        plot2d.plot_median2d(mesh, input_currents, [], title, figname, pdf_file, save=save,linewidth=3)
        plt.plot(smooth_t[:,0],smooth_t[:,1], 'k', linewidth = 3.0)
        plt.savefig('%s.png'%figname, dpi=300)
        
        figname = '%s/figures/%d'%(outdir, figcount)
        plot2d.plot_decomposition2d(mesh, input_currents, t, q_list[i], r_list[i], title, \
        figname, pdf_file, save, linewidth=3)
        figcount += input_currents.shape[0]
        
       
        mesh, simplices, subsimplices, input_points, lambdas, mus, is_closed = cases2d.multicurves2d(n=smooth_t.shape[0]) 
        figname = '%s/figures/%d'%(outdir, figcount)
        plt.clf()
        fig = plt.gca().figure
        plt.gca().set_aspect('equal')
        colors = plot2d.get_colors(len(input_currents))
        for i, input_point in enumerate(input_points):
            plt.plot(input_point[:,0],input_point[:,1], color=colors[i], linewidth = 3.0)
        plt.plot(smooth_t[:,0],smooth_t[:,1], 'k', linewidth = 3.0)
        plt.savefig('%s.png'%figname, dpi=300)
        pdf_file.savefig(fig)

        plt.show()



    if pdf_file is not None:
        pdf_file.close()
    
if __name__ == '__main__':
    mediandemo2d()
