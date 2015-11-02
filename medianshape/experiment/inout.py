
import os.path

import numpy as np
from scipy import sparse

from medianshape.simplicial.mesh import Mesh2D, Mesh3D
import medianshape.utils as utils 


# Loads previously computed mesh, boundary_matrix, input currents, w and v from a directory
def load_mesh2d(dirname='data/dumps'):
    '''
    Hi
    '''
    mesh = Mesh2D()
    if os.path.exists(dirname):
        if os.path.exists("%s/points.txt"%dirname):
            mesh.points = np.loadtxt("%s/points.txt"%dirname, dtype=np.float) 
            mesh.bbox = (np.amin(mesh.points[:,0]),\
                            np.amin(mesh.points[:,1]),\
                            np.amax(mesh.points[:,0]),\
                            np.amax(mesh.points[:,1]))

            mesh.set_boundary_points()
            mesh.set_diagonal()
            mesh.set_boundary_values()
        else:
            print "Can't load points. <points.txt> file doesn't exist"
        if os.path.exists("%s/edges.txt"%dirname):
            mesh.edges= np.loadtxt("%s/edges.txt"%dirname, dtype=np.int) 
        else:
            print "Can't load edges. <edges.txt> file doesn't exist"
        if os.path.exists("%s/simplices.txt"%dirname):
            mesh.simplices = np.loadtxt("%s/simplices.txt"%dirname, dtype=np.int) 
        else:
            print "Can't load simplices. <simplices.txt> file doesn't exist"
    else:
        print "%s directory doesn't exist"%dirname
    return mesh

def load_mesh3d(dirname='data/dumps'):
    '''
    Hi
    '''
    mesh = Mesh3D()
    if os.path.exists(dirname):
        if os.path.exists("%s/points.txt"%dirname):
            mesh.points = np.loadtxt("%s/points.txt"%dirname, dtype=np.float) 
            mesh.bbox = (np.amin(mesh.points[:,0]),\
                            np.amin(mesh.points[:,1]),\
                            np.amin(mesh.points[:,2]),\
                            np.amax(mesh.points[:,0]),\
                            np.amax(mesh.points[:,1]),\
                            np.amax(mesh.points[:,2]))

            mesh.set_boundary_points()
            mesh.set_diagonal()
            mesh.set_boundary_values()
        else:
            print "Can't load points. <points.txt> file doesn't exist"
        if os.path.exists("%s/edges.txt"%dirname):
            mesh.edges= np.loadtxt("%s/edges.txt"%dirname, dtype=np.int) 
        else:
            print "Can't load edges. <edges.txt> file doesn't exist"
        if os.path.exists("%s/simplices.txt"%dirname):
            mesh.simplices = np.loadtxt("%s/simplices.txt"%dirname, dtype=np.int) 
        else:
            print "Can't load simplices. <simplices.txt> file doesn't exist"
        if os.path.exists("%s/triangles.txt"%dirname):
            mesh.triangles = np.loadtxt("%s/triangles.txt"%dirname, dtype=np.int)
        else:
            mesh.triangles = utils.get_subsimplices(mesh.simplices)
            print "Can't load triangles. <triangles.txt> file doesn't exist"
    else:
        print "%s directory doesn't exist"%dirname
    return mesh

def load_weights_and_boundary(n_simplices, m_subsimplices, dirname='data/dumps'):
    '''
    Hi
    '''
    w = np.zeros(shape=(m_subsimplices, 1))
    v = np.zeros(shape=(n_simplices, 1))
    b_matrix = sparse.dok_matrix((m_subsimplices, n_simplices), dtype=np.int8)
    if os.path.exists("%s/w.txt"%dirname):
        w = np.loadtxt("%s/w.txt"%dirname, dtype=np.float) 
    else:
        print "Can't load the weight vector, w. <%s/w.txt> file doesn't exist"%dirname
    if os.path.exists("%s/v.txt"%dirname):
        v = np.loadtxt("%s/v.txt"%dirname, dtype=np.float) 
    else:
        print "Can't load the weight vector, v. <%s/v.txt> file doesn't exist"%dirname
    if os.path.exists("%s/b_matrix.txt"%dirname):
        with open("%s/b_matrix.txt"%dirname, 'r') as f:
            for line in f.readlines():
                data = line.split()
                b_matrix[int(data[0]), int(data[1])] = np.int8(data[2])
    else:
        print "Can't load boundary matrix. <%s/b_matrix.txt> file doesn't exist"%dirname
    return w, v, b_matrix

def load_input_currents(m_subsimplices, k, dirname='data/dumps'):
    '''
    Hi
    '''
    input_currents = np.zeros(shape=(k, m_subsimplices), dtype=np.int) 
    for i in range(k):
        if os.path.isfile('%s/input_current%d.txt'%(dirname, i)):
            with open("%s/input_current%d.txt"%(dirname, i), 'r') as f:
                for line in f.readlines():
                    data = line.split()
                    input_currents[i, int(data[1])] = np.int(data[2])
        else:
            print "Can't load input current %d. <input_current%d.txt> doesn't exist"%(i, i)
    return input_currents

def load_solutions(n_simplices, m_subsimplices, k, dirname='data/dumps'):
    '''
    Hi
    '''
    x = np.zeros(shape=(2*m_subsimplices + 2*k*m_subsimplices + 2*k*n_simplices,1), dtype=np.int)
    t = np.zeros((m_subsimplices, 1), dtype=np.int) 
    q = np.zeros((k, m_subsimplices), dtype=int)
    r = np.zeros((k, n_simplices), dtype=int)
    if os.path.exists('%s/x.txt'%dirname):
        with open("%s/x.txt"%dirname, 'r') as f:
            for line in f.readlines():
                data = line.split()
                x[int(data[0])] = np.int(data[1])
        t = x[0:m_subsimplices] - x[m_subsimplices:2*m_subsimplices]
        qi_start = 2*m_subsimplices
        for i in range(k):
            qi_end = qi_start + 2*m_subsimplices
            q[i] = (x[qi_start: qi_start+m_subsimplices] - x[qi_start+m_subsimplices: qi_end]).reshape(m_subsimplices,)
            ri_start = qi_end
            ri_end = ri_start + 2*n_simplices
            r[i] = (x[ri_start: ri_start+n_simplices] - x[ri_start+n_simplices: ri_end]).reshape(n_simplices, )
            qi_start = ri_end
    else:
        print "Can't load the solution. <x.txt> file doesn't exist"
    return t, q, r
    
# Saves mesh, input currents, boundary matrix, w and v.
def save_data(mesh=None, input_currents=None, b_matrix=None, w=None, v=None, t=None, dirname=os.path.abspath('data/dumps'), **kwargs):
    '''
    Hi
    '''
    if mesh is not None:
        np.savetxt('%s/edges.txt' % dirname, mesh.edges, fmt='%d', delimiter=' ')
        np.savetxt('%s/simplices.txt'% dirname, mesh.simplices, fmt='%d', delimiter=' ')
        if mesh.points.shape[1] == 3:
            np.savetxt('%s/triangles.txt'% dirname, mesh.triangles, fmt='%d', delimiter=' ')
    if input_currents is not None:
        for i, c in enumerate(input_currents):
            sparse_savetxt('%s/input_current%d.txt' % (dirname,i), c)
    if b_matrix is not None:
        sparse_savetxt('%s/b_matrix.txt' % dirname, b_matrix)
    if w is not None:
        np.savetxt('%s/w.txt' % dirname, w, delimiter=' ')
    if v is not None:
        np.savetxt('%s/v.txt' % dirname, v, delimiter=' ')
    if t is not None:
        if 'mu' in kwargs and 'lambda_' in kwargs:
            sparse_savetxt("%s/t-lambda-%s-mu-%s.txt"%(dirname, kwargs['lambda_'], kwargs['mu']), t)
        else:
            sparse_savetxt("%s/t.txt"%dirname, t)
# Saves sparse matrix as text. if the input is not sparse, set is_sparse argument to False.
def sparse_savetxt(fname, matrix, fmt='%d', delimiter=' '):
    '''
    Hi
    '''
    if sparse.issparse(matrix):
        if matrix.getformat() !='coo':
            matrix = matrix.asformat('coo')
    else:
        matrix = sparse.coo_matrix(matrix)
    with open(fname, 'w') as f:
        for i in range(len(matrix.row)):
            f.write("%d %d %d\n" % (matrix.row[i], matrix.col[i], matrix.data[i]))
