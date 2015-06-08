#encoding: utf-8

from __future__ import absolute_import

import importlib
import time

import numpy as np

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

from mesh.utils import boundary_matrix, simpvol, get_subsimplices
from shapegen import currentgen, utils
from runmedians import runmedians3d

import plot3d
import cases3d

def show_median3d():
    mesh, input_currents, t, q, r = cases3d.fine_curves_on_sphere()
    fig = plt.figure(figsize=(8,8))
    plot3d.plotmesh3d(mesh, mesh.get_info())
    fig.tight_layout()
    plt.show()
    plt.figure(figsize=(8,8))
    plot3d.plot_median3d(mesh, input_currents, t, title="Median")
    plt.show()
    plot3d.plot_decomposition3d(mesh, input_currents, t, q, r, title='Flatnorn decomposition')


def mediandemo3d(outdir='output', save=True):
    lp_times = list()
    start = time.time()
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/mediandemo3d.pdf'%outdir)

    fig = plt.figure(figsize=(8,8))
    figcount = 1
    mesh, simplices, subsimplices, points, lambdas, mus, is_closed \
    = cases3d.equally_spaced_longitudes3d()
    print mesh.get_info()
    plot3d.plotmesh3d(mesh, mesh.get_info(), '%s/figures/mesh'%outdir, pdf_file, save)
    fig.tight_layout()
    plt.show()
    fig = plt.figure(figsize=(8,8))
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, simplices, subsimplices, points, is_closed=is_closed)

    figname = '%s/figures/%d.png'%(outdir, figcount)
    title = mesh.get_info()
    title = 'Curve approximation'
    plot3d.plot_curves_approx3d(mesh, points, vertices, paths, title, figname, pdf_file, save)
    figcount += 1
    fig.tight_layout()
    plt.show()
    runmedians3d(mesh, simplices, subsimplices, input_currents, lambdas, mus, file_doc=pdf_file, save=save)
    if save:
        pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    #mediandemo3d(save=False)
    show_median3d()
