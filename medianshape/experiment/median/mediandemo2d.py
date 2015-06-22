# encoding: utf-8
import os
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from medianshape.simplicial import currentgen
from medianshape.experiment.median import runmedians as run, cases2d

from medianshape.viz import plot2d 

def mediandemo2d(outdir='data/output', save=True):
    lp_times = list()
    start = time.time()
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/mediandemo2d.pdf'%outdir)
    fig = plt.figure(figsize=(8,8))
    figcount = 1

    mesh, simplices, subsimplices, points, lambdas, mus, is_closed \
    = cases2d.multicurves2d() 
    print mesh.get_info()

    #vertices, paths, input_currents = currentgen.push_functions_on_mesh_2d(mesh, points, is_closed=False, functions=functions)
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, simplices, subsimplices, points, is_closed=is_closed)

    figname = '%s/figures/%d.png'%(outdir, figcount)
    title = mesh.get_info()
    plot2d.plot_curves_approx2d(mesh, points, vertices, paths, title, figname, pdf_file, save)
    plt.show()
    fig = plt.figure(figsize=(8,8))
    figcount += 1

    run.runmedians2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, file_doc=pdf_file, save=save)

    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mediandemo2d()
