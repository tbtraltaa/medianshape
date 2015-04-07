# encoding: utf-8

from __future__ import absolute_import

import itertools

import matplotlib.pyplot as plt

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
        plt.show()
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
        plt.show()
    else:
        print "Plotting only supported in dimensions 2 and 3." 

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
    lim = float(mesh.boundary_box[3])/10
    plt.ylim([mesh.boundary_box[1]-lim, mesh.boundary_box[3]+lim])
    plt.xlim([mesh.boundary_box[0]-lim, mesh.boundary_box[2]+lim])
    mesh.plot()
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
    lim = float(mesh.boundary_box[3])/10
    plt.ylim([mesh.boundary_box[1]-lim, mesh.boundary_box[3]+lim])
    plt.xlim([mesh.boundary_box[0]-lim, mesh.boundary_box[2]+lim])
    mesh.plot()
    for i, c in enumerate(input_currents):
        mesh.plot_curve(c, color=colors.next(), label='current%d, %d'%(i+1, comb[i]), linewidth=5)
    mesh.plot_curve(t, title)
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
    fig = plt.gca().figure
    for i, c in enumerate(input_currents):
        fig.clf()                    
        plt.gca().set_aspect('equal')
        lim = mesh.boundary_box[3]/10
        plt.ylim([mesh.boundary_box[1]-lim, mesh.boundary_box[3]+lim])
        plt.xlim([mesh.boundary_box[0]-lim, mesh.boundary_box[2]+lim])
        mesh.plot()
        mesh.plot_curve(c, color=colors.next(), linewidth=5, \
        label='current%d, %d'%(i+1, comb[i]))
        mesh.plot_curve(t, title, label='Mean')
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
    print r
    for i, r_i in enumerate(r):
        color = colors.next()
        fig.clf()
        plt.gca().set_aspect('equal')
        lim = mesh.boundary_box[3]/10
        plt.ylim([mesh.boundary_box[1]-lim, mesh.boundary_box[3]+lim])
        plt.xlim([mesh.boundary_box[0]-lim, mesh.boundary_box[2]+lim])
        mesh.plot()
        if r_dim == 1:
            mesh.plot_curve(r_i, title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
            plt.scatter(mesh.points[q[i]][:, 0], mesh.points[q[i]][:,1], color='r')
        elif r_dim ==2:
            mesh.plot_simplices(r_i, color=color)
            mesh.plot_curve(q[i], title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
            if t is not None:
                mesh.plot_curve(t, linewidth=4, label="Mean")
            if i < input_currents.shape[0]:
                mesh.plot_curve(input_currents[i], color='r', ls='--', \
                label='current%d, %d'%(i+1, comb[i]))
        plt.legend(loc='upper right')
        if save and figname:
            plt.savefig("%s-%d.png" % (figname, i), dpi=fig.dpi)
        if save and file_doc:
            file_doc.savefig(fig)

