# encoding: utf-8

from __future__ import absolute_import

import importlib
import time


import numpy as np
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from shapegen import currentgen
import plot2d

from mesh.utils import boundary_matrix, simpvol, get_subsimplices
from runmedians import rundeform2d

import cases2d

def deform2d(outdir='output', save=True):
    lp_times = list()
    start = time.time()
    figcount = 1
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/deform2d.pdf'%outdir)
    fig = plt.figure(figsize=(14,4))
    mesh, simplices, subsimplices, points, lambdas, mus, alphas, is_closed \
    = cases2d.ellipsesdeform2d() 
    print mesh.get_info()
    print alphas

    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, simplices, subsimplices, points, is_closed=is_closed)
    figname = '%s/figures/%d.png'%(outdir, figcount)
    title = mesh.get_info()
    plot2d.plot_curves_approx2d(mesh, points, vertices, paths, title, figname, pdf_file, save=save)
    plt.tight_layout()
    plt.show()
    fig = plt.figure(figsize=(14,4))
    figcount += 1

    rundeform2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, alphas, file_doc=pdf_file, save=save)
    if save:
        pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    deform2d()
