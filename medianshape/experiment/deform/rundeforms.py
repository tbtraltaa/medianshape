
from __future__ import absolute_import 

import numpy as np
import matplotlib.pyplot as plt

from medianshape.core import median
import medianshape.utils as utils 
from medianshape.viz import plot2d, plot3d
import medianshape.experiment.inout as inout

def rundeform2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, alphas, w=None, v=None, b_matrix=None, file_doc=None, save=True, outdir='data/output'):
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

                #figname = '%s/figures/%d-%.04f-%.04f'%(outdir, figcount, l, mu)
                #plot2d.plot_curve_and_median(mesh, input_currents, t, title, \
                #figname, file_doc, save)
                #figcount += input_currents.shape[0]

                #figname = '%s/figures/%d-%.06f-%.06f'%(outdir, figcount, l, mu)
                #plot2d.plot_decomposition(mesh, input_currents, t, q, r, title, \
                #figname, file_doc, save)
                #figcount += input_currents.shape[0]

def rundeforms3d(mesh, input_currents, options, lambdas, mus, alphas, w=None, v=None, b_matrix=None, file_doc=None, save_data=True):
    '''
    Hi
    '''
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
        w, v, b_matrix, cons = median.get_lp_inputs(mesh.points, mesh.triangles, mesh.edges,  k_currents, opt, w, v, b_matrix)
        #np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
        for l in lambdas:
            comb=[1,1,1]
            #for comb in combinations[:-1,:]:
            #for comb in combinations:
                #input_currents = currents*comb.reshape(comb.size,1) 
            for mu in mus:
                for alpha in alphas:
                    t, q, r, norm = median.median(mesh.points, mesh.triangles, mesh.edges, input_currents, l, opt, w, v, cons, mu=mu, alphas=alpha)
                    if save_data:
                        save(t=t, opt=opt, lambda_=l)
                    norms.append(norm)
                    t_len = len(t.nonzero()[0])
                    t_lens.append(t_len)
                    title = '%s, lambda=%.04f, mu=%.06f'  % \
                    (opt, l, mu)
                    figname = 'output/figures/%d-%s-%.04f-%.06f'%(figcount, opt, l, mu)
                    if save and file_doc is not None:
                        plot3d.plot_median(mesh, input_currents, comb, t, title, figname, file_doc, save=save)
                        #plt.show()
                        figcount += 1

                        #figname = 'output/figures/%d-%s-%.04f-%.04f'%(figcount, opt, l, mu)
                        #plotting.plot_curve_and_median(mesh, input_currents, comb, t, title, \
                        #figname, file_doc, save)
                        #figcount += input_currents.shape[0]

                        figname = 'output/figures/%d-%s-%.06f-%.06f'%(figcount,opt,l, mu)
                        plot3d.plot_decomposition(mesh, input_currents, comb, t, q, r, title, \
                        figname, file_doc, save)
                        figcount += input_currents.shape[0]
                        #plt.show()
                
                # Plotting the combination with minimum flatnorm difference
            #title = 'Minimum flatnorm difference, %s, lambda=%.04f, %s' % (opt, l, str(comb))
#                figname = 'output/figures/%d-%s-%.04f'%(figcount, opt, l)
#                plotting.plot_median(mesh, min_currents, min_comb, min_t, title, \
#                figname, file_doc)
#                figcount += 1
#
#                figname = 'output/figures/%d-%s-%.04f'%(figcount,opt,l)
#                plotting.plot_decomposition(mesh, min_currents, comb, min_t, min_q, min_r, \
#                title, figname, file_doc)
#                figcount += input_currents.shape[0]
