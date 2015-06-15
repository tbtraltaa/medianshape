# encoding: utf-8

from __future__ import absolute_import

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm

def get_colors(n):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct 
        RGB color.'''
    if n < 4:
        if n == 1:
            return 'r'
        elif n == 2:
            return 'gr'
        elif n == 3:
            return 'gry'
        elif n == 4:
            return 'grcy' 
    else:
        return cm.rainbow(np.linspace(0, 1, n))

def plotmesh2d(mesh):
    plt.triplot(mesh.points[:,0], mesh.points[:,1], mesh.simplices.copy())
    plt.tight_layout()
    #plt.scatter(mesh.points[:,0], mesh.points[:,1])

def plot_curve2d(mesh, func_path, title=None, color="black", marker=None, linewidth=3, ls='-', label=""):
    if type(func_path) == list:
        func_path = np.array(func_path)
    if func_path.dtype != int:
        func_path = func_path.astype(int)
    nonzero_edges = func_path.nonzero()[0]
    for i, edge_idx in enumerate(nonzero_edges):
        edge = mesh.edges[edge_idx]
        points = mesh.points[edge]
        plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls, label=label if i==0 else "")

    if title:
        plt.title(title)

# Plot simplices
def plot_simplices2d(mesh, simplices, title=None, color="y", label=""):
    ax = plt.gca()
    #ccw_symbol = u'\u2941'
    #cw_symbol = u'\u21BB'
    for i, idx in enumerate(simplices.nonzero()[0]):
        hatch = ''
        if simplices[idx] == -1:
            hatch = '.' 
        simplex = plt.Polygon(mesh.points[mesh.simplices[idx]], closed=True, fill=True, fc=color, hatch=hatch, label=label if i==0 else "")
        ax.add_patch(simplex)

def plot_curves_approx2d(mesh, points, vertices, paths, title="", figname=None, file_doc=None, save=True, lim=5):
    colors = get_colors(len(points))
    fig = plt.gca().figure
    plt.clf()
    plt.gca().set_aspect('equal')
    lim = float(mesh.ymax)/10
    plt.ylim([mesh.ymin-lim, mesh.ymax+lim])
    plt.xlim([mesh.xmin-lim, mesh.xmax+lim])
    plotmesh2d(mesh)
    for i, path in enumerate(paths):
        plot_curve_approx2d(mesh, points[i], vertices[i], path, color=colors[i], label="T%d"%(i+1))
    plt.title(title)
    plt.tight_layout()
    if save and figname:
        plt.savefig('%s.png'%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_approx2d(mesh, input_points, closest_vertices, path, title=None, color="red", linewidth=3, label=""):
    plt.plot(input_points[:,0], input_points[:,1], c=color, ls="--", label='Input points')
    plt.title(title)
    for i, edge in enumerate(path):
        points = mesh.points[edge]
        plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, label=label if i==0 else "")
    plt.scatter(mesh.points[closest_vertices][:,0], mesh.points[closest_vertices][:,1], s=100, label='Closest vertices')
    plt.scatter(input_points[:,0], input_points[:,1], c=color)
    plt.legend(loc='lower right')

def plot_median2d(mesh, input_currents, t, title='', figname="", file_doc=None, save=True, lim=5):
    colors= get_colors(len(input_currents))
    plt.clf()
    fig = plt.gca().figure
    plt.gca().set_aspect('equal')
    lim = float(mesh.ymax)/10
    plt.ylim([mesh.ymin-lim, mesh.ymax+lim])
    plt.xlim([mesh.xmin-lim, mesh.xmax+lim])
    plotmesh2d(mesh)
    for i, c in enumerate(input_currents):
        plot_curve2d(mesh,  c, color=colors[i], label='T%d'%(i+1), linewidth=5)
    plot_curve2d(mesh, t, title, label='Median')
    plt.legend(loc='lower right')
    plt.tight_layout()
    if save and figname:
        plt.savefig("%s.png"%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_and_median2d(mesh, input_currents, t, title=None, figname=None, file_doc=None, save=True, lim=5):
    colors= get_colors(len(input_currents))
    fig = plt.gca().figure
    for i, c in enumerate(input_currents):
        fig.clf()                    
        plt.gca().set_aspect('equal')
        lim = mesh.ymax*1.0/10
        plt.ylim([mesh.ymin-lim, mesh.ymax+lim])
        plt.xlim([mesh.xmin-lim, mesh.xmax+lim])
        plotmesh2d(mesh)
        plot_curve2d(mesh, c, color=colors[i], linewidth=5, \
        label='T%d'%(i+1))
        plot_curve2d(mesh, t, title, label='median')
        plt.legend(loc='lower right')
        plt.tight_layout()
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)

def plot_decomposition2d(mesh, input_currents, t, q, r, title='', figname=None, file_doc=None, save=True, lim=None):
    colors= get_colors(len(input_currents))
    fig = plt.gca().figure
    for i, r_i in enumerate(r):
        fig.clf()
        plt.gca().set_aspect('equal')
        if lim is None:
            lim = mesh.ymax*1.0/10
        print lim
        plt.ylim([mesh.ymin-lim, mesh.ymax+lim])
        plt.xlim([mesh.xmin-lim, mesh.xmax+lim])
        plotmesh2d(mesh)
        plot_simplices2d(mesh, r_i, color=colors[i], label="R%d"%(i+1))
        if q is not None:
            plot_curve2d(mesh, q[i], title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
        plot_curve2d(mesh, input_currents[i], color='r', ls='--', label='T%d'%(i+1))
        if t is not None:
            plot_curve2d(mesh, t, linewidth=4, label="Median")
        plt.legend(loc='lower right')
        plt.tight_layout()
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)
        plt.show()
        fig = plt.figure(figsize=(8,8))

