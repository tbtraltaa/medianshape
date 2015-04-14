# encoding: utf-8

from __future__ import absolute_import

import itertools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotmesh3d(mesh):
    dim = mesh.points.shape[1]
    if dim == 2:
        from distmesh.plotting import SimplexCollection
        ax = plt.gca()
        c = SimplexCollection()
        ax.add_collection(c)
        ax.set_xlim(mesh.bbox[2])
        ax.set_ylim(mesh.bbox[3])
        ax.set_aspect('equal')
        ax.set_axis_off()
        c.set_simplices((mesh.points, mesh.simplices))
    elif dim == 3:
        import mpl_toolkits.mplot3d
        from distmesh.plotting import axes_simpplot3d
        ax = plt.gca(projection='3d')
        ax.set_xlim(mesh.bbox[3])
        ax.set_ylim(mesh.bbox[4])
        ax.set_zlim(mesh.bbox[5])
        ax.cla()
        axes_simpplot3d(ax, mesh.points, mesh.simplices, mesh.points[:,1] > 0)
        ax.set_title('Retriangulation')
    else:
        print "Plotting only supported in dimensions 2 and 3." 

def plot_curve(mesh, func_path, title=None, color="black", marker=None, linewidth=3, ls='-', label=""):
    ax = plt.gca(projection='3d')
    if type(func_path) == list:
        func_path = np.array(func_path)
    if func_path.dtype != int:
        func_path = func_path.astype(int)
    nonzero_edges = func_path.nonzero()[0]
    for i, edge_idx in enumerate(nonzero_edges):
        edge = mesh.edges[edge_idx]
        points = mesh.points[edge]
        if i == len(nonzero_edges)-1:
            ax.plot(points[:,0], points[:,1], points[:,2], color, linewidth=linewidth, marker=marker, ls=ls, label=label)
        else:
            ax.plot(points[:,0], points[:,1], points[:,2], color, linewidth=linewidth, marker=marker, ls=ls)
    if title:
        plt.title(title)

# Plot simplices
def plot_simplices(mesh, simplices, title=None, color="y"):
    ax = plt.gca(projection='3d')
    #ccw_symbol = u'\u2941'
    #cw_symbol = u'\u21BB'
    for i in simplices.nonzero()[0]:
        hatch = ''
        if simplices[i] == -1:
            hatch = '.' 
        simplex = plt.Polygon(mesh.points[self.simplices[i]], closed=True, fill=True, fc=color, hatch=hatch)
        ax.add_patch(simplex)

def plot_curves_approx(mesh, points, vertices, paths, title="", figname=None, file_doc=None, save=True, lim=5):
    color_set = "r"
    if len(paths) == 2:
        color_set = 'gr'
    elif len(paths) == 3:
        color_set = 'gry'
    elif len(paths) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    fig = plt.gca(projection='3d').figure
    plt.clf()
    plt.gca().set_aspect('equal')
    lim = float(mesh.bbox[3])/10
    plt.ylim([mesh.bbox[1]-lim, mesh.bbox[4]+lim])
    plt.xlim([mesh.bbox[0]-lim, mesh.bbox[3]+lim])
    plotmesh3d(mesh)
    for i, path in enumerate(paths):
        plot_curve_approx(mesh, points[i], vertices[i], path, color=colors.next())
    plt.title(title)
    if save and figname:
        plt.savefig('%s.png'%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_approx(mesh, input_points, closest_vertices, path, title=None, color="red", linewidth=3, label=""):
    ax = plt.gca(projection='3d')
    ax.plot(input_points[:,0], input_points[:,1], input_points[:,2], c=color, ls="--")
    plt.title(title)
    for i, edge in enumerate(path):
        points = mesh.points[edge]
        ax.plot(points[:,0], points[:,1], points[:,2],color, linewidth=linewidth, label=label)
    ax.scatter(mesh.points[closest_vertices][:,0], mesh.points[closest_vertices][:,1], mesh.points[closest_vertices][:,2], s=100)
    ax.scatter(input_points[:,0], input_points[:,1], input_points[:,2], c=color)


def plot_mean(mesh, input_currents, comb, t, title='', figname="", file_doc=None, save=True, lim=5):
    color_set = "r"
    if len(input_currents) == 2:
        color_set = 'gr'
    elif len(input_currents) == 3:
        color_set = 'gry'
    elif len(input_currents) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    plt.clf()
    fig = plt.gca(projection='3d').figure
    plt.gca().set_aspect('equal')
    lim = float(mesh.bbox[3])/10
    plt.ylim([mesh.bbox[1]-lim, mesh.bbox[4]+lim])
    plt.xlim([mesh.bbox[0]-lim, mesh.bbox[3]+lim])
    #mesh.plot()
    for i, c in enumerate(input_currents):
        plot_curve(mesh, c, color=colors.next(), label='current%d, %d'%(i+1, comb[i]), linewidth=5)
    plot_curve(mesh, t, title)
    plt.legend(loc='upper right')
    if save and figname:
        plt.savefig("%s.png"%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_and_mean(mesh, input_currents, comb, t, title=None, figname=None, file_doc=None, save=True, lim=5):
    color_set = "r"
    if len(functions) == 2:
        color_set = 'gr'
    elif len(functions) == 3:
        color_set = 'gry'
    elif len(functions) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    ax = plt.gca(projection='3d')
    fig = ax.figure
    for i, c in enumerate(input_currents):
        fig.clf()                    
        plt.gca().set_aspect('equal')
        lim = mesh.bbox[3]/10
        plt.ylim([mesh.bbox[1]-lim, mesh.bbox[4]+lim])
        plt.xlim([mesh.bbox[0]-lim, mesh.bbox[3]+lim])
        plot_curve(mesh, c, color=colors.next(), linewidth=5, \
        label='current%d, %d'%(i+1, comb[i]))
        plot_curve(mesh, t, title, label='Mean')
        plt.legend(loc='upper right')
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)

def plot_decomposition(mesh, input_currents, comb, t, q, r, title='', figname=None, file_doc=None, save=True, lim=5):
    color_set = "r"
    if len(input_currents) == 2:
        color_set = 'gr'
    elif len(input_currents) == 3:
        color_set = 'gry'
    elif len(input_currents) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    ax = plt.gca(projection='3d')
    fig = ax.figure
    for i, r_i in enumerate(r):
        color = colors.next()
        fig.clf()
        plt.gca().set_aspect('equal')
        lim = mesh.bbox[3]/10
        plt.ylim([mesh.bbox[1]-lim, mesh.bbox[4]+lim])
        plt.xlim([mesh.bbox[0]-lim, mesh.bbox[3]+lim])
        plot_simplices(mesh, r_i, color=color)
        plot_curve(mesh, q[i], title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
        if t is not None:
            plot_curve(mesh, t, linewidth=4, label="Mean")
        if i < input_currents.shape[0]:
            plot_curve(mesh, input_currents[i], color='r', ls='--', \
            label='current%d, %d'%(i+1, comb[i]))
        plt.legend(loc='upper right')
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)
