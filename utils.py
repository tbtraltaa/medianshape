from __future__ import absolute_import 

import numpy as np
from scipy.sparse import dok_matrix

def extract_edges(simplices):
    edges = set()
    for simplex in simplices:
        for i in range(len(simplex)):
            edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
    return np.array(list(edges))

def sparse_savetxt(fname, matrix, fmt='%d', delimiter=' ', is_sparse=True):
    if not is_sparse:
        matrix = coo_matrix(matrix)
    else:
        matrix = matrix.asformat('coo')
    with open(fname, 'w') as f:
        for i in range(len(matrix.row)):
            f.write("%d %d %d\n" % (matrix.row[i], matrix.col[i], matrix.data[i]))
