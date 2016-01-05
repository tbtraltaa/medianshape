**Getting started**
===================

It is a python package called Mass Regularized Simplicial Median Shape (MRSMS). The purpose of the package is to compute the median shape of input shapes given. The inputs can be any dimensional, closed or open such as curves, 2D shapes or 3D objects and they is represented by **currents** which is a mathematical concepts adapted from Geometric Measure Theory(GMT). We used **flat norm** which is also a concept from GMT to measure distances between the median current and input currents. The input currents are approximated in an underlying simplicial complex and the *simplicial flat norm* used to compare the simplicial currents with the median. Finding the median which locates within minimum distance from each input currents is an optimization problem and we constructed the median shape problem as a Linear Programming (LP) problem and solved it using LP solvers such as *cplex* and *cvxopt*. There are several experiments included in the package which works on 1-currents in 2D (the codimension one case) and 3D (the codimension two case). MRSMS itself is dimension-free and the same MRSMS runs for both 2D and 3D experiments. That is the generality adapted from Geometric Measure Theory. 

Install
-------
``pip install medianshape``

Requirements
------------

*    Babel==1.3
*    Cython==0.21.1
*    Jinja2==2.7.3
*    MarkupSafe==0.23
*    MeshPy==2014.1
*    PuLP==1.5.6
*    PyDistMesh==1.2
*    PyVTK==0.4.85
*    Pygments==2.0.2
*    Sphinx==1.3.1
*    alabaster==0.7.3
*    argparse==1.2.1
*    cplex==12.5.1.0
*    cvxopt==1.1.7
*    decorator==3.4.0
*    docutils==0.12
*    matplotlib==1.4.0
*    mock==1.0.1
*    nose==1.3.4
*    numpy==1.9.0
*    ply==3.4
*    py==1.4.24
*    pyparsing==1.5.7
*    pytest==2.6.2
*    python-dateutil==2.2
*    pytools==2014.3
*    pytz==2014.10
*    scipy==0.14.0
*    six==1.8.0
*    snowballstemmer==1.2.0
*    sphinx-rtd-theme==0.1.7
*    wsgiref==0.1.2

Please refer to their own documentation for different version compatability.
You can use **Anaconda** tool to easily install numpy, scipy and scikit-learn. 
Anaconda is a collection python packages for scientific computation and it provides
a package manager and an environment manager additionally.
