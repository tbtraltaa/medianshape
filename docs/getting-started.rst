**Getting started**
===================

It is a python package called **Mass Regularized Simplicial Median Shape (MRSMS)**. The purpose of the package is to compute the median shape of input shapes given. The inputs can be any dimensional, closed or open such as curves, 2D shapes or 3D objects and they is represented by **currents** which is a mathematical concepts adapted from Geometric Measure Theory(GMT). We used **flat norm** which is also a concept from GMT to measure distances between the median current and input currents. The input currents are approximated in an underlying simplicial complex and the *simplicial flat norm* used to compare the simplicial currents with the median. Finding the median which locates within minimum distance from each input currents is an optimization problem and we constructed the median shape problem as a **Linear Program(LP) or Integer Optimization problem** and solved it using LP solvers such as **cplex** and **cvxopt**. There are several experiments included in the package which works on 1-currents in 2D (the codimension one case) and 3D (the codimension two case). MRSMS itself is dimension-free and the same MRSMS runs for both 2D and 3D experiments. That is the generality adapted from Geometric Measure Theory. 

Install
-------
``pip install medianshape``

Run source code using conda package management
----------------------------------------------

Linux(Ubuntu16.04)
------------

* Install conda
* conda create -name medianshape python=2
* source ~/.bashrc
* source activate medianshape
* conda install scipy=0.17.1
* conda install numpy=1.11.1
* conda install matplotlib=1.5.1
* conda install cvxopt=1.1.8

<!--Linux 32-bit
------------

sudo apt-get install lib32ncurses5
sudo apt-get install lib32z1
-->

Installing PyDistMesh1.2
------------------------
* conda install cython=0.24.1 #needed for pydistmesh
* sudo apt-get install libblas-dev liblapack-dev
* pip install pydistmesh

Export the `medianshape` library location to PYTHONPATH
-------------------------------------------------------
To use Medianshape library, put the following lines in the files, ~/.profile or ~/.bashrc
* export PYTHONPATH=$PYTHONPATH:$HOME/medianshape/

Install cplex
-------------
To install the CPLEX-Python modules on your system, use the script setuy.py located in yourCplexhome/python/PLATFORM. If you want to install the CPLEX-Python modules in a nondefault location, use the option --home to identify the installation directory. For example, to install the CPLEX-Python modules in the default location, use the following command from the command line:

* sudo chmode a+x cplex_studio1251.linux-x86-32.bin
* sudo ./cplex_studio1251.linux-x86-32.bin
* sudo chmod -R 777 <cplex_HOME_folder> - cplex lacks permission to a create folder when you install it from conda environment.
* source activate medianshape
* python setup.py install

cplex throws error when row and col args are not explicitly typecasted to int.
prob.linear_constraints.set_coefficients(zip(cons.row.astype(int), cons.col.astype(int), cons.data.astype(float)))

Run tetview
-----------
* Download tetview-linux
* gunzip tetview-linux.gz
* chmod +x tetview-linux
* sudo apt-get install libglu1-mesa:i386
* Tetview requires 32-bit version of the package
* Tetview doesn't run on Ubuntu 16.05 64-bit
* Install Wine
* Download tetview-win.exe
* Run it using Wine

Install tetgen
--------------

* sudo apt-get install tetgen

Documentation
-------------
* conda install sphinx
* source activate medianshape
* cd medinashape/docs
* make html
* google-chrome medianshape/docs/_build/html/index.html

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
