# encoding: utf-8
from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt

from mesh.utils import boundary

class Mesh():
    def __init(self):
        self.points = None
        self.simplices = None
        self.edges = None

    def set_edges(self):
        edges = set()
        for simplex in self.simplices:
            for i in range(len(simplex)):
                edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
        self.edges = np.array(list(edges))

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

    def to_string(self):
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
        return plt

    def plot_curve(self, func_path, title=None, color="black"):
        for i, orient in enumerate(func_path):
            if orient != 0:
                edge = self.edges[i]
                points = self.points[edge]
                plt.plot(points[:,0], points[:,1], color)
        if title:
            plt.title(title)
        return plt
