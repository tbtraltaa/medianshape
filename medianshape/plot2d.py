# encoding: utf-8

from __future__ import absolute_import

import itertools

import matplotlib.pyplot as plt

def plot(mesh):
    plt.triplot(mesh.points[:,0], mesh.points[:,1], mesh.simplices.copy())
    #plt.scatter(mesh.points[:,0], mesh.points[:,1])

def plot_curve(mesh, func_path, title=None, color="black", marker=None, linewidth=3, ls='-', label=""):
    if type(func_path) == list:
        func_path = np.array(func_path)
    if func_path.dtype != int:
        func_path = func_path.astype(int)
    nonzero_edges = func_path.nonzero()[0]
    for i, edge_idx in enumerate(nonzero_edges):
        edge = mesh.edges[edge_idx]
        points = mesh.points[edge]
        if i == len(nonzero_edges)-1:
            plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls, label=label)
        else:
            plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls)

    if title:
        plt.title(title)

# Plot simplices
def plot_simplices(mesh, simplices, title=None, color="y"):
    ax = plt.gca()
    #ccw_symbol = u'\u2941'
    #cw_symbol = u'\u21BB'
    for i in simplices.nonzero()[0]:
        hatch = ''
        if simplices[i] == -1:
            hatch = '.' 
        simplex = plt.Polygon(mesh.points[mesh.simplices[i]], closed=True, fill=True, fc=color, hatch=hatch)
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
    fig = plt.gca().figure
    plt.clf()
    plt.gca().set_aspect('equal')
    lim = float(mesh.bbox[3])/10
    plt.ylim([mesh.bbox[1]-lim, mesh.bbox[3]+lim])
    plt.xlim([mesh.bbox[0]-lim, mesh.bbox[2]+lim])
    plot(mesh)
    for i, path in enumerate(paths):
        plot_curve_approx(mesh, points[i], vertices[i], path, color=colors.next())
    plt.title(title)
    if save and figname:
        plt.savefig('%s.png'%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_approx(mesh, input_points, closest_vertices, path, title=None, color="red", linewidth=3, label=""):
    plt.plot(input_points[:,0], input_points[:,1], c=color, ls="--")
    plt.title(title)
    for i, edge in enumerate(path):
        points = mesh.points[edge]
        plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, label=label)
    plt.scatter(mesh.points[closest_vertices][:,0], mesh.points[closest_vertices][:,1], s=100)
    plt.scatter(input_points[:,0], input_points[:,1], c=color)

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
    fig = plt.gca().figure
    plt.gca().set_aspect('equal')
    lim = float(mesh.bbox[3])/10
    plt.ylim([mesh.bbox[1]-lim, mesh.bbox[3]+lim])
    plt.xlim([mesh.bbox[0]-lim, mesh.bbox[2]+lim])
    plot(mesh)
    for i, c in enumerate(input_currents):
        plot_curve(mesh,  c, color=colors.next(), label='current%d, %d'%(i+1, comb[i]), linewidth=5)
    plot_curve(mesh, t, title, label='Median curve')
    plt.legend(loc='upper right')
    if save and figname:
        plt.savefig("%s.png"%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_and_mean(mesh, input_currents, comb, t, title=None, figname=None, file_doc=None, save=True, lim=5):
    color_set = "r"
    if len(input_currents) == 2:
        color_set = 'gr'
    elif len(input_currents) == 3:
        color_set = 'gry'
    elif len(input_currents) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    fig = plt.gca().figure
    for i, c in enumerate(input_currents):
        fig.clf()                    
        plt.gca().set_aspect('equal')
        lim = mesh.bbox[3]/10
        plt.ylim([mesh.bbox[1]-lim, mesh.bbox[3]+lim])
        plt.xlim([mesh.bbox[0]-lim, mesh.bbox[2]+lim])
        plot(mesh)
        plot_curve(mesh, c, color=colors.next(), linewidth=5, \
        label='current%d, %d'%(i+1, comb[i]))
        plot_curve(mesh, t, title, label='Mean')
        plt.legend(loc='upper right')
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)

def plot_decomposition(mesh, input_currents, comb, t, q, r, title='', figname=None, file_doc=None, save=True, lim=5, r_dim=2):
    color_set = "r"
    if len(input_currents) == 2:
        color_set = 'gr'
    elif len(input_currents) == 3:
        color_set = 'gry'
    elif len(input_currents) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    fig = plt.gca().figure
    for i, r_i in enumerate(r):
        color = colors.next()
        fig.clf()
        plt.gca().set_aspect('equal')
        lim = mesh.bbox[3]/10
        plt.ylim([mesh.bbox[1]-lim, mesh.bbox[3]+lim])
        plt.xlim([mesh.bbox[0]-lim, mesh.bbox[2]+lim])
        plot(mesh)
        if r_dim == 1:
            plot_curve(mesh, r_i, title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
            plt.scatter(mesh.points[q[i]][:, 0], mesh.points[q[i]][:,1], color='r')
        elif r_dim ==2:
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

