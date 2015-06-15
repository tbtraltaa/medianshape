
import numpy as np
import matplotlib.pyplot as plt

import medianshape.core.median as median
from medianshape import utils
from medianshape.viz import plot2d, plot3d
import medianshape.experiment.inout as inout

def runmedians2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, w=None, v=None, b_matrix=None, file_doc=None, save=True, outdir='../data/output'):
    figcount = 2
    if w is None:
        w = utils.simpvol(mesh.points, subsimplices)
    if v is None:
        v = utils.simpvol(mesh.points, simplices)
    if b_matrix is None:
        b_matrix = utils.boundary_matrix(simplices, subsimplices)
    k_currents = len(input_currents)
    w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, simplices, subsimplices,  k_currents, w, v, b_matrix)
    for l in lambdas:
        for mu in mus:
            t, q, r, norm = median.median(mesh.points, mesh.simplices, mesh.edges, \
            input_currents, l, w, v, cons, mu=mu)
            if save:
                inout.save_data(t=t, lambda_=l, mu=mu)
            title = 'MRSMS, lambda=%.06f, mu=%.06f'  % (l, mu)
            figname = '%s/figures/%d-%.06f-%.06f'%(outdir, figcount, l, mu)
            plot2d.plot_median2d(mesh, input_currents, t, title, figname, file_doc, save=save)
            plt.show()
            fig = plt.figure(figsize=(8,8))
            figcount += 1

            #figname = '%s/figures/%d-%.04f-%.04f'%(outdir, figcount, l, mu)
            #plot2d.plot_curve_and_median2d(mesh, input_currents, t, title, \
            #figname, file_doc, save=save)
            #figcount += input_currents.shape[0]

            figname = '%s/figures/%d-%.06f-%.06f'%(outdir, figcount, l, mu)
            plot2d.plot_decomposition2d(mesh, input_currents, t, q, r, title, \
            figname, file_doc, save)
            figcount += input_currents.shape[0]
    if save:
        inout.save_data(mesh, input_currents, b_matrix, w, v)

def runmedians3d(mesh, simplices, subsimplices, input_currents, lambdas, mus, w=None, v=None, b_matrix=None, file_doc=None, save=True, outdir='../data/output'):
    figcount = 2
    if w is None:
        w = utils.simpvol(mesh.points, subsimplices)
    if v is None:
        v = utils.simpvol(mesh.points, simplices)
    if b_matrix is None:
        b_matrix = utils.boundary_matrix(simplices, subsimplices)
    k_currents = len(input_currents)
    w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, simplices, subsimplices,  k_currents, w, v, b_matrix)
    #np.savetxt('%s/dumps/cons-%s.txt'%cons, fmt='%d', delimiter=' ')
    for l in lambdas:
        for mu in mus:
            t, q, r, norm = median.median(mesh.points, simplices, subsimplices, input_currents, l, w, v, cons, mu=mu)
            if save:
                inout.save_data(t=t, lambda_=l, mu=mu)
            title = 'MRSMS, lambda=%.06f, mu=%.06f'%(l, mu)
            figname = '%s/figures/%d-%.06f-%.06f'%(outdir, figcount, l, mu)
            fig = plt.figure(figsize=(8,8))
            plot3d.plot_median3d(mesh, input_currents, t, title, figname, file_doc, save)
            plt.tight_layout()
            plt.show()
            fig = plt.figure(figsize=(8,8))
            figcount += 1

            #figname = '%s/figures/%d-%s-%.04f-%.04f'%(outdir, figcount, l, mu)
            #plotting.plot_curve_and_median3d(mesh, input_currents, comb, t, title, \
            #figname, file_doc, save)
            #figcount += input_currents.shape[0]

            figname = '%s/figures/%d-%.06f-%.06f'%(outdir, figcount, l, mu)
            plot3d.plot_decomposition3d(mesh, input_currents, t, q, r, title, \
            figname, file_doc, save)
            figcount += input_currents.shape[0]
    if save:
        inout.save_data(mesh, input_currents, b_matrix, w, v)
