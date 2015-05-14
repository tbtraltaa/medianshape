# encoding: utf-8

from __future__ import absolute_import

import importlib
import random
import time

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from shapegen import currentgen, utils
import mesh.utils as mutils
from runmedians import runmedians2d

import plot2d
import cases2d

def mediandemo2d(outdir='output', save=True):
    lp_times = list()
    start = time.time()
    pdf_file = None
    if save:
        pdf_file = PdfPages('%s/mediandemo2d.pdf'%outdir)
    fig = plt.figure(figsize=(8,8))
    figcount = 1
    boundary_box = (0,0,200,50)
    l=6
    #boundary_box = (0,0,40,40)
    #fixed_points = [(0,0),(40,0),(0,40),(40,40)]

    #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
    #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
    #functions= ['curve4', 'curve5']
    functions= ['curve1', 'curve2']
    #combinations = np.array([[1,1,1]])
    #curve1 = pointgen2d.sample_curve1()
    #curve2 = pointgen2d.sample_curve2()
    #mesh.plot()
    #plt.scatter(curve1[:,0], curve1[:,1], color='r')
    #plt.scatter(curve2[:,0], curve2[:,1], color='k')
    #plt.show()
    #shapes = [curve1, curve2]
    #plot2d.plot_curve(mesh, c)
    #points = list()
    #for f in functions:
     #   points.append(pointgen2d.sample_function_mesh(mesh, f))
    mesh, simplices, subsimplices, points = cases2d.ellipses() 
    print mesh.get_info()
    #vertices, paths, input_currents = currentgen.push_functions_on_mesh_2d(mesh, points, is_closed=False, functions=functions)
    vertices, paths, input_currents = currentgen.push_curves_on_mesh(mesh, simplices, subsimplices, points, True)
    figname = '%s/figures/%d.png'%(outdir, figcount)
    title = mesh.get_info()
    plot2d.plot_curves_approx2d(mesh, points, vertices, paths, title, figname, pdf_file, save)
    plt.show()
    fig = plt.figure(figsize=(8,8))
    figcount += 1
    lambdas = [0.001]
    lambdas = [1]
    mus = [0.0001]

    runmedians2d(mesh, simplices, subsimplices, input_currents, lambdas, mus, file_doc=pdf_file, save=save)

    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    
if __name__ == '__main__':
    mediandemo2d()
