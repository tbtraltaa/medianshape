# encoding: utf-8
'''
----------------------
Median surface demo 3D
−−−−−−−−−−−−−−−−--−−−−
'''
from __future__ import absolute_import
import importlib
import time
from shutil import copyfile

import numpy as np
import matplotlib.pyplot as plt
from medianshape.core import median

def load_tetgen_mesh(fname):
    '''
    Loads .node, .face, .ele files. Usually their indexing starts from 1.
    '''
    points = np.loadtxt("%s.node"%fname, dtype=float, skiprows=1)[:,1:4]
    triangles = np.loadtxt("%s.face"%fname, dtype=int, skiprows=1)[:,1:5]
    tetras = np.loadtxt("%s.ele"%fname, dtype=int, skiprows=1)[:,1:5]
    return points, triangles, tetras

def save_tetgen_mesh(t, q, r, fname="/home/altaa/tet/tet.1"):
    copyfile("%s.node"%fname, "%s.1.node"%fname)
    triangles = np.loadtxt("%s.face"%fname, dtype=int, skiprows=1)
    tetras = np.loadtxt("%s.ele"%fname, dtype=int, skiprows=1)
    print triangles.shape
    if len(np.argwhere(t!=0))!=0:
        triangles[np.argwhere(t!=0)][4] = -3
    '''
    flag = -3
    for tmp in q:
        if len(np.argwhere(tmp!=0))!=0:
            triangles[np.argwhere(tmp!=0)][4] = flag -1
    for tmp in r:
        if len(np.argwhere(tmp!=0))!=0:
            r[np.argwhere(tmp!=0)] = flag-1
    print tetras.shape
    print r.shape
    tetras = np.hstack((tetras,r))
    np.savetxt("%s.1.ele"%fname, header="%d %d %d"%(len(tetras), 4, 1))
    '''
    np.savetxt("%s.1.face"%fname, triangles, fmt="%d", header="%d %d"%(len(triangles), 1))
    print t.shape
    print np.unique(t)
    print np.unique(q)
    print np.unique(r)

def surfaces3d(fname="/home/altaa/tet/tet.1"):
    points, faces, tetras = load_tetgen_mesh(fname)
    inputsurf1 = (faces[:,3] == -1).astype(int)
    inputsurf2 = (faces[:,3] == -2).astype(int)
    inputcurrents = np.vstack((inputsurf1, inputsurf2))
    print faces
    triangles = faces[:,:-1] - 1
    print triangles
    tetras = tetras - 1
    lambda_ = 0.0000001
    mu = 0.000001
    return points, tetras, triangles, inputcurrents, lambda_, mu

def mediansurfdemo3d(outdir='data/output', save=True):
    '''
    Median shape demo for median surface
    '''
    start = time.time()
    points, simplices, subsimplices, inputcurrents, _lambda, mu = surfaces3d()
    w, v, b_matrix, cons = median.get_lp_inputs(points, simplices, subsimplices, len(inputcurrents))
    t, q, r, norm = median.median(points, simplices, subsimplices, \
    inputcurrents, _lambda, w, v, cons, mu=mu)
    elapsed = time.time() - start
    print 'Elapsed time %f mins.' % (elapsed/60)
    save_tetgen_mesh(t, q, r)

if __name__ == '__main__':
    mediansurfdemo3d(save=True)
