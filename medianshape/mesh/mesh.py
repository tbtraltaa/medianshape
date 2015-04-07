# encoding: utf-8
from __future__ import absolute_import

import numpy as np
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

from mesh.utils import boundary

class Mesh2D():
    '''

    '''
    def __init__(self, *args, **kwargs):
        self.points = kwargs.get('points', None)
        self.simplices = kwargs.get('simplices', None)
        self.edges = kwargs.get('edges', None)
        self.bbox = kwargs.get('bbox', None)
        self.fixed_points = kwargs.get('fixed_points', None)
        self.diagonal = kwargs.get('diagonal', None)
        if self.bbox is not None:
            if self.fixed_points is None: 
                self.set_fixed_points()
            if self.diagonal is None:
                self.set_diagonal()

    def set_edges(self):
        edges = set()
        for simplex in self.simplices:
            for i in range(len(simplex)):
                edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
        self.edges = np.array(list(edges), dtype=int)

    def set_fixed_points(self):
        self.fixed_points = list()
        if self.bbox is not None:
            if len(self.bbox) == 4:
                self.fixed_points = [[self.bbox[0], self.bbox[1]],\
                                        [self.bbox[0], self.bbox[3]],\
                                        [self.bbox[2], self.bbox[1]],\
                                        [self.bbox[2], self.bbox[3]]]
    def set_diagonal(self):
        self.diagonal = np.sqrt(self.bbox[2]**2 + mesh.bbox[3]**2)

    # get simplices based on simplices vector
    def get_simplices(self, simplices_vector, opts=None):
        if simplices_vector.dtype != int:
            simplices_vector = simplices_vector.astype(int)
        # All simplices despite their orientations
        if not opts:
            simplices = self.simplices[simplices_vector.nonzero()[0]]
        # Counter clockwise simplices
        elif opts == 'ccw':
            simplices = self.simplices[np.where(simplices_vector == 1)[0]] 
        # Clockwise simplices
        elif opts == 'cw':
            simplices = self.simplices[np.where(simplices_vector == -1)[0]]
        return simplices 

    def orient_simplices_2D(self):
        for i, simplex in enumerate(self.simplices):
            if self.right_hand_rule(simplex) < 0: 
                self.simplices[i] = self.simplices[i,::-1]

    def right_hand_rule(self, simplex):
        edges = list()
        edges.append(boundary(simplex, 0))
        edges.append(boundary(simplex, 1))
        edge_points = self.points[np.array(edges)]
        v1 = edge_points[0][1] - edge_points[0][0]
        v2 = edge_points[1][1] - edge_points[1][0]
        direction = np.cross(v1,v2)
        return direction

    def get_info(self):
        return "Mesh info: %d points, %d triangles and %d edges"% \
        (len(self.points), len(self.simplices), len(self.edges))

    def print_detail(self):
        print self.get_info()
        print "Mesh Points:"
        for i, p in enumerate(self.points):
                print i, p
        print "Point numbers in triangle:"
        for i, t in enumerate(self.simplices):
                print i, t
        print "Edges in mesh:"
        for i, edge in enumerate(self.edges):
            print i, edge

    def plot(self):
        plt.triplot(self.points[:,0], self.points[:,1], self.simplices.copy())
        #plt.scatter(self.points[:,0], self.points[:,1])

    def plot_curve(self, func_path, title=None, color="black", marker=None, linewidth=3, ls='-', label=""):
        if type(func_path) == list:
            func_path = np.array(func_path)
        if func_path.dtype != int:
            func_path = func_path.astype(int)
        nonzero_edges = func_path.nonzero()[0]
        for i, edge_idx in enumerate(nonzero_edges):
            edge = self.edges[edge_idx]
            points = self.points[edge]
            if i == len(nonzero_edges)-1:
                plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls, label=label)
            else:
                plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls)

        if title:
            plt.title(title)

    # Plot simplices
    def plot_simplices(self, simplices, title=None, color="y"):
        ax = plt.gca()
        #ccw_symbol = u'\u2941'
        #cw_symbol = u'\u21BB'
        for i in simplices.nonzero()[0]:
            hatch = ''
            if simplices[i] == -1:
                hatch = '.' 
            simplex = plt.Polygon(self.points[self.simplices[i]], closed=True, fill=True, fc=color, hatch=hatch)
            ax.add_patch(simplex)

class Mesh3D():
    '''

    '''
    def __init__(self, *args, **kwargs):
        self.points = kwargs.get('points', [])
        self.simplices = kwargs.get('simplices', [])
        self.edges = kwargs.get('edges', [])
        self.triangles = kwargs.get('triangles', [])
        self.bbox = kwargs.get('bbox', [])
        self.fixed_points = kwargs.get('fixed_points', [])
        self.diagonal = kwargs.get('diagonal', 1e10)
        if len(self.bbox) != 0:
            self.set_diagonal()

    def set_fixed_points(self):
        self.fixed_points = list()
        if len(self.bbox) != 0 and len(self.bbox) == 6:
            self.fixed_points = [[self.bbox[0], self.bbox[1], self.bbox[2]],\
                                    [self.bbox[0], self.bbox[1], self.bbox[5]],\
                                    [self.bbox[0], self.bbox[4], self.bbox[2]],\
                                    [self.bbox[0], self.bbox[4], self.bbox[5]],\
                                    [self.bbox[3], self.bbox[1], self.bbox[2]],\
                                    [self.bbox[3], self.bbox[4], self.bbox[2]],\
                                    [self.bbox[3], self.bbox[1], self.bbox[5]],\
                                    [self.bbox[3], self.bbox[4], self.bbox[5]]]
    def set_diagonal(self):
        self.diagonal = pdist([[self.bbox[0],self.bbox[1], self.bbox[2]],\
                                [self.bbox[3], self.bbox[4], self.bbox[5]]])

    # get simplices based on simplices vector
    def get_simplices(self, simplices_vector, opts=None):
        if simplices_vector.dtype != int:
            simplices_vector = simplices_vector.astype(int)
        # All simplices despite their orientations
        if not opts:
            simplices = self.simplices[simplices_vector.nonzero()[0]]
        # Counter clockwise simplices
        elif opts == 'ccw':
            simplices = self.simplices[np.where(simplices_vector == 1)[0]] 
        # Clockwise simplices
        elif opts == 'cw':
            simplices = self.simplices[np.where(simplices_vector == -1)[0]]
        return simplices 

    def get_info(self):
        return "Mesh info: %d points, %d tetrahedras, %d triangles and %d edges"% \
        (len(self.points), len(self.simplices), len(self.triangles), len(self.edges))

    def print_detail(self):
        print self.get_info()
        print "Mesh Points:"
        for i, p in enumerate(self.points):
                print i, p
        print "Point numbers in tetrahedras:"
        for i, t in enumerate(self.simplices):
                print i, t
        print "Point numbers in triangles:"
        for i, t in enumerate(self.triangles):
                print i, t
        print "Edges in mesh:"
        for i, edge in enumerate(self.edges):
            print i, edge

    def plot_curve(self, func_path, title=None, color="black", marker=None, linewidth=3, ls='-', label=""):
        if type(func_path) == list:
            func_path = np.array(func_path)
        if func_path.dtype != int:
            func_path = func_path.astype(int)
        nonzero_edges = func_path.nonzero()[0]
        for i, edge_idx in enumerate(nonzero_edges):
            edge = self.edges[edge_idx]
            points = self.points[edge]
            if i == len(nonzero_edges)-1:
                plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls, label=label)
            else:
                plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, marker=marker, ls=ls)

        if title:
            plt.title(title)

    # Plot simplices
    def plot_simplices(self, simplices, title=None, color="y"):
        ax = plt.gca()
        #ccw_symbol = u'\u2941'
        #cw_symbol = u'\u21BB'
        for i in simplices.nonzero()[0]:
            hatch = ''
            if simplices[i] == -1:
                hatch = '.' 
            simplex = plt.Polygon(self.points[self.simplices[i]], closed=True, fill=True, fc=color, hatch=hatch)
            ax.add_patch(simplex)
