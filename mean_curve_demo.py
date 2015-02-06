# encoding: utf-8

from __future__ import absolute_import

import math
import importlib
import random
import itertools
import time


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from scipy.sparse import csr_matrix

from mesh.distmesh import distmesh2d
from mesh.mesh import Mesh
from shape_gen import point_gen, curve_gen, utils
from mesh.utils import boundary_matrix, simpvol

import mean
from cvxopt import matrix, solvers
from matplotlib.backends.backend_pdf import PdfPages

options = ['default', 'mass', 'msfn']
#options = ['default']

if __name__ == '__main__':
    start = time.time()
    mesh = Mesh()
    # l - initial length of triangle sides 
    # change it to 1 for big traingles
    mesh.boundary_box = (0,0,200,50)
    mesh.fixed_points = [(0,0),(200,0),(0,50),(200,50)]
    mesh.points, mesh.simplices = distmesh2d('square', mesh.boundary_box, mesh.fixed_points, l=6)
    mesh.set_edges()
    np.savetxt('/home/altaa/dumps1/points.txt', mesh.points, delimiter=' ')
    np.savetxt('/home/altaa/dumps1/edges.txt', mesh.edges, fmt='%d', delimiter=' ')
    np.savetxt('/home/altaa/dumps1/simplices.txt', mesh.simplices, fmt='%d', delimiter=' ')

#    mesh.points = np.array([[0,0], [0,0.5],[0,1],[0.5,1],[1,1],[1,0.5],[1,0],[0.5,0], [0.5, 0.5]])
#    mesh.simplices = np.array([[0,8,1],
#                                [1,8,2],
#                                [2,8,3],
#                                [3,8,4],
#                                [4,8,5],
#                                [5,8,6],
#                                [6,8,7],
#                                [7,8,0]])
#
#    mesh.edges = np.array([[0,1],
#                            [0,7],
#                            [0,8],
#                            [1,2],
#                            [1,8],
#                            [2,3],
#                            [2,8],
#                            [3,4],
#                            [3,8],
#                            [4,5],
#                            [4,8],
#                            [5,6],
#                            [5,8],
#                            [6,7],
#                            [6,8],
#                            [7,8]])

    mesh.orient_simplices_2D()
    mesh.print_detail()

    pdf_file = PdfPages('/home/altaa/figures1.pdf')
    #function_sets = [['sin1pi','half_sin1pi'], ['x', 'x2', 'x5']]
    #function_sets = [['curve1', 'curve2', 'curve3', 'curve4', 'curve5']]
    function_sets = [['curve1', 'curve2', 'curve3']]

    figcount = 1
    for j, functions in enumerate(function_sets):
        combinations = utils.get_combination(len(functions))
        if len(functions) == 2:
            color_set = 'gr'
        elif len(functions) == 3:
            color_set = 'gry'
        elif len(functions) == 5:
            color_set = 'grcym' 
        colors = itertools.cycle(color_set)
        input_currents = list()
        fig = plt.figure(figsize=(19,8))
        plt.gca().set_aspect('equal')
        plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+5])
        plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
        mesh.plot()
        #w = np.ndarray(shape=(len(mesh.edges),))
        #w[0:] = 1

        #v =  np.ndarray(shape=(len(mesh.simplices),))
        #v[0:] = 0.433

        for i, f in enumerate(functions):
            points = point_gen.sample_function_mesh(f, mesh)
            input_current, path, closest_vertices = curve_gen.generate_curve_on_mesh(points, mesh, func_str=f) 
            np.savetxt('/home/altaa/dumps1/%s.txt'%f, input_current.reshape(len(input_current),1), fmt='%d', delimiter=' ')
            csr_path = csr_matrix(input_current)
            print 'Path vector:\n', csr_path
            curve_gen.plot_curve(mesh, points, closest_vertices, path, color=colors.next())
            input_currents.append(input_current)
        k_currents = len(functions)
        currents = np.array(input_currents).reshape(k_currents, mesh.edges.shape[0])
        w = simpvol(mesh.points, mesh.edges)
        v = simpvol(mesh.points, mesh.simplices)
        b_matrix = boundary_matrix(mesh.simplices, mesh.edges)
        np.savetxt('/home/altaa/dumps1/w.txt', w, delimiter=' ')
        np.savetxt('/home/altaa/dumps1/v.txt', v, delimiter=' ')
        np.savetxt('/home/altaa/dumps1/b_matrix.txt', b_matrix, fmt='%d', delimiter=' ')
        
        plt.title('Functions - %s - (%s)' % (mesh.get_info(), ','.join(functions)), fontsize=20)
        figname = '/home/altaa/fig_dump/%d.png'%(figcount)
        plt.savefig(figname, dpi=fig.dpi)
        figcount += 1
        pdf_file.savefig(fig)
        for opt in options:
            w, v, b_matrix, cons = mean.get_lp_inputs(mesh.points, mesh.simplices, mesh.edges, k_currents,
            opt, w, v, b_matrix)
            np.savetxt('/home/altaa/dumps1/cons-%s.txt'%opt, cons, fmt='%d', delimiter=' ')
            lambdas = [0.0001, 0.001]
            for l in lambdas:
#            input_currents = list()
#            current1 = np.zeros(shape=(len(mesh.edges),1))
#            current1[0] = 1
#            current1[3] = 1
#            current1[5] = 1 
#            current1[7] = 1
#            current2 = np.zeros(shape=(len(mesh.edges),1))
#            current2[1] = 1
#            current2[13] = -1
#            current2[11] = -1
#            current2[9] = -1
#            input_currents.append(current1)
#            input_currents.append(current2)
#            input_currents = np.array(input_currents)
#            k_currents = 2
                for comb in combinations[:-1,:]:
                    input_currents = currents*comb.reshape(comb.size,1) 
                    x, q, r, norm = mean.mean(mesh.points, mesh.simplices, mesh.edges, input_currents, l, opt, w, v, cons)
                    fig.clf()
                    plt.gca().set_aspect('equal')
                    plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+20])
                    plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
                    rows = 2
                    if opt == 'msfn' and len(functions) == 3:
                        rows = 3
                    elif len(functions) == 6:
                        rows = 4
                    #plt.subplot(rows, 2, 1)
                    mesh.plot()
                    for i, c in enumerate(input_currents):
                        mesh.plot_curve(c, color=colors.next(), label='%s, %d'%(functions[i], comb[i]), linewidth=5)
                    title = '%s, lambda=%.04f, %s' % (opt, l, str(comb))
                    mesh.plot_curve(x, title)
                    plt.legend(loc='upper right')
                    figname = '/home/altaa/fig_dump/%d-%s-%.04f.png'%(figcount, opt, l)
                    plt.savefig(figname, dpi=fig.dpi)
                    pdf_file.savefig(fig)
                    figcount += 1
#                    colors = itertools.cycle(color_set)
#                    for i, c in enumerate(input_currents):
#                        fig.clf()                    
#                        plt.gca().set_aspect('equal')
#                        plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+15])
#                        plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
#                        mesh.plot()
#                        mesh.plot_curve(c, color=colors.next(), linewidth=5, \
#                        label='%s, %d'%(functions[i], comb[i]))
#                        #title = '%s, lambda=%.04f, (%s,%d)' % (opt, l, functions[i], comb[i])
#                        mesh.plot_curve(x, title, label='Mean')
#                        plt.legend(loc='upper right')
#                        figname = '/home/altaa/fig_dump/%d-%s-%.04f.png'%(figcount, opt, l)
#                        plt.savefig(figname, dpi=fig.dpi)
#                        pdf_file.savefig(fig)
#                        figcount += 1
#                    for i in range(r.shape[1]):
#                        color = colors.next()
#                        fig.clf()
#                        plt.gca().set_aspect('equal')
#                        plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+20])
#                        plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
#                        mesh.plot()
#                        mesh.plot_simplices(r[:,i], color=color)
#                        mesh.plot_curve(q[:,i], title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
#                        mesh.plot_curve(x, linewidth=4, label="Mean")
#                        if opt =='msfn' and i== r.shape[1]-1:
#                            pass
#                        else:
#                            mesh.plot_curve(input_currents[i], color='r', ls='--', \
#                            label='%s, %d'%(functions[i], comb[i]))
#                        plt.legend(loc='upper right')
#                        figname = '/home/altaa/fig_dump/%d-%s-%.04f.png'%(figcount,opt,l)
#                        plt.savefig(figname, dpi=fig.dpi)
#                        pdf_file.savefig(fig)
#                        figcount += 1

                    #figname = '/home/altaa/fig_dump/%d-%s-%s-%.04f.png'%(figcount, '-'.join(functions),opt,l)
                    #plt.savefig(figname, dpi=fig.dpi)
                    #figcount += 1
                        #print 'q1', q1
                        #print 'r1', r1
                        #print 'q2', q2
                        #print 'r2', r2
    pdf_file.close()
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
