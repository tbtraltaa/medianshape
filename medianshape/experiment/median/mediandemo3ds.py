#encoding: utf-8
'''
Median shape demo 3D
--------------------
'''

from __future__ import absolute_import

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from medianshape.simplicial import currentgen 
from medianshape.experiment.median import runmedians as run, cases3d

from medianshape.viz import plot3d 

def mediandemo3ds(outdir='data/output', save=True):
    '''
    Hi
    '''
    lp_times = list()
    start = time.time()
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/mediandemo3ds.pdf'%outdir)

    fig = plt.figure(figsize=(8,8))
    figcount = 1
    mesh, smesh, simplices, subsimplices, points, lambdas, mus, is_closed \
    = cases3d.equally_spaced_longitudes3ds()
    print smesh.get_info()
    title = r"$%d$ $points$, $%d$ $triangles$ $and$ $%d$ $edges$"% (len(smesh.surf_points), \
    len(smesh.triangles), len(smesh.edges))
    plot3d.plotmesh3d(mesh, title, '%s/figures/%d'%(outdir, figcount), pdf_file, save)
    figcount += 1
    fig.tight_layout()
    #plt.show()
    fig = plt.figure(figsize=(8,8))
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(smesh, simplices, subsimplices, points, is_closed, smesh.surf_points)

    figname = '%s/figures/%d'%(outdir, figcount)
    plot3d.plot_curves_approx3d(smesh, points, vertices, paths, figname, pdf_file, save)
    figcount += 1
    fig.tight_layout()
   # plt.show()
    run.runmedians3d(smesh, simplices, subsimplices, input_currents, lambdas, mus, file_doc=pdf_file, save=save, figcount=figcount)
    if save:
        pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)

if __name__ == '__main__':
    mediandemo3ds(save=True)
    #show_median3d()
