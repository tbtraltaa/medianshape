'''
++++++++++++++++++++
Running median shape
++++++++++++++++++++

'''
from __future__ import absolute_import 

import numpy as np
import matplotlib.pyplot as plt

from medianshape.core import median
import medianshape.utils as utils 
from medianshape.viz import plot2d, plot3d
import medianshape.experiment.inout as inout

def rundeform2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, alphas, w=None, v=None, b_matrix=None, file_doc=None, save=True, outdir='data'):
    '''
    Hi
    '''
    figcount = 2
    if w is None:
        w = utils.simpvol(mesh.points, subsimplices)
    if v is None:
        v = utils.simpvol(mesh.points, simplices)
    if b_matrix is None:
        b_matrix = utils.boundary_matrix(simplices, subsimplices)   
    k_currents = len(input_currents)
    w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, simplices, subsimplices,  k_currents, w, v, b_matrix)
    #np.savetxt('%s/dumps/cons.txt', cons, fmt='%d', delimiter=' ')
    for l in lambdas:
        for mu in mus:
            for alpha in alphas:
                t, q, r, norm = median.median(mesh.points, simplices, subsimplices, \
                input_currents, l, w, v, cons, mu=mu, alphas=alpha)
                if save:
                    inout.save_data(t=t, lambda_=l, mu=mu)
                title = r'$MRSMS$, $\lambda=%.06f$, $\mu=%.06f$, $\alpha=[%.06f$ $%.06f]$'%(l, mu, alpha[0], alpha[1])
                figname = '%s/figures/%d'%(outdir, figcount)
                plot2d.plot_median2d(mesh, input_currents, t, title, figname, file_doc, save=save)
                figcount += 1

def rundeform3d(mesh, simplices, subsimplices, input_currents, lambdas, mus, alphas, w=None, v=None, b_matrix=None, file_doc=None, save=True, figcount=1, outdir='data'):
    '''
    Hi
    '''
    figcount = figcount
    if w is None:
        w = utils.simpvol(mesh.points, mesh.edges)
    if v is None:
        v = utils.simpvol(mesh.points, mesh.triangles)
    if b_matrix is None:
        b_matrix = utils.boundary_matrix(mesh.triangles, mesh.edges)
    k_currents = len(input_currents)
    print k_currents
    w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, simplices, subsimplices,  k_currents, w, v, b_matrix)
    #np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
    for l in lambdas:
        comb=[1,1,1]
        #for comb in combinations[:-1,:]:
        #for comb in combinations:
            #input_currents = currents*comb.reshape(comb.size,1) 
        for mu in mus:
            for alpha in alphas:
                t, q, r, norm = median.median(mesh.points, mesh.triangles, mesh.edges, input_currents, l, w, v, cons, mu=mu, alphas=alpha)
                if save:
                    inout.save_data(t=t, lambda_=l, mu=mu)
                title = r'$MRSMS$, $\lambda=%.06f$, $\mu=%.06f$, $\alpha=[%.06f$ $%.06f]$'%(l, mu, alpha[0], alpha[1])
                figname = '%s/figures/%d'%(outdir, figcount)
                if save and file_doc is not None:
                    fig = plt.figure(figsize=(8,8))
                    plot3d.plot_median3d(mesh, input_currents, t, title, figname, file_doc, save=save)
                    figcount += 1
