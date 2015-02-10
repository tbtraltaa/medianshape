# encoding: utf-8

from __future__ import absolute_import

import itertools

import matplotlib.pyplot as plt

def plot_curves_approx(mesh, points, vertices, paths, title="", figname=None, file_doc=None, save=True):
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
    plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+5])
    plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
    mesh.plot()
    for i, path in enumerate(paths):
        plot_curve_approx(mesh, points[i], vertices[i], path, color=colors.next())
    plt.title(title, fontsize=20)
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

def plot_mean(mesh, functions, input_currents, comb, t, title, figname="", file_doc=None, save=True):
    if len(functions) == 2:
        color_set = 'gr'
    elif len(functions) == 3:
        color_set = 'gry'
    elif len(functions) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    plt.clf()
    fig = plt.gca().figure
    plt.gca().set_aspect('equal')
    plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+20])
    plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
    mesh.plot()
    for i, c in enumerate(input_currents):
        mesh.plot_curve(c, color=colors.next(), label='%s, %d'%(functions[i], comb[i]), linewidth=5)
    mesh.plot_curve(t, title)
    plt.legend(loc='upper right')
    if save and figname:
        plt.savefig("%s.png"%figname, dpi=fig.dpi)
    if save and file_doc:
        file_doc.savefig(fig)

def plot_curve_and_mean(mesh, functions, input_currents, comb, t, title=None, save=True, figname=None, file_doc=None):
    if len(functions) == 2:
        color_set = 'gr'
    elif len(functions) == 3:
        color_set = 'gry'
    elif len(functions) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    fig = plt.gca().figure
    for i, c in enumerate(input_currents):
        plt.clf()                    
        plt.gca().set_aspect('equal')
        plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+15])
        plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
        mesh.plot()
        mesh.plot_curve(c, color=colors.next(), linewidth=5, \
        label='%s, %d'%(functions[i], comb[i]))
        mesh.plot_curve(t, title, label='Mean')
        plt.legend(loc='upper right')
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)

def plot_decomposition(mesh, functions, input_currents, comb, t, q, r, title=None, figname=None, file_doc=None, save=True):
    if len(functions) == 2:
        color_set = 'gr'
    elif len(functions) == 3:
        color_set = 'gry'
    elif len(functions) == 5:
        color_set = 'grcym' 
    colors = itertools.cycle(color_set)
    fig = plt.gca().figure
    for i, r_i in enumerate(r):
        color = colors.next()
        fig.clf()
        plt.gca().set_aspect('equal')
        plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+20])
        plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
        mesh.plot()
        mesh.plot_simplices(r_i, color=color)
        mesh.plot_curve(q[i], title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
        mesh.plot_curve(t, linewidth=4, label="Mean")
        if i < input_currents.shape[0]:
            mesh.plot_curve(input_currents[i], color='r', ls='--', \
            label='%s, %d'%(functions[i], comb[i]))
        plt.legend(loc='upper right')
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)

