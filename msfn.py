# encoding: utf-8

from __future__ import absolute_import

import numpy as np

from mesh.utils import boundary_matrix, simpvol

def msfn(points, simplices, subsimplices, input_current, lambda, v=None, w=None, cons=None):
    if not w:
        w = simpvol(points, subsimplices)
    if not v:
        v = simpvol(points, simplices)
    if not cons:
        m_edges = mesh.edges.shape[0]
        n_simplices = mesh.simplices.shape[0]
        boundary_matrix = boundary_matrix(simplices, subsimplices)
        cons = np.hstack((np.identity(m_edges), -np.identity(m_edges)))
        cons = np.hstack((cons, boundary_matrix))
        cons = np.hstack((cons, -boundary_matrix))

    
    
