**Getting started**
===================

MedianShape_ is a python package called **Mass Regularized Simplicial Median Shape (MRSMS)**. The purpose of the package is to compute the median of shapes. Input shapes can be any dimensional, closed or open such as curves, 2D shapes or 3D objects represented by **currents** which is a mathematical concepts adapted from Geometric Measure Theory(GMT). We used **flat norm** which is also a concept from GMT to measure distances between a median current and input currents. Input currents are approximated in an underlying simplicial complex and *simplicial flat norm* is used to compare simplicial currents with median. Finding a median which locates within minimum distance from each input current is an optimization problem and we constructed the median shape problem as a **Linear Program(LP)** by relaxing **Integer Optimization problem** and solved it using LP solvers such as **cplex** and **cvxopt**. There are several experiments included in the package which works on 1-currents in 2D (codimension one cases) and 3D (codimension two cases). MRSMS itself is dimension-free and the same MRSMS runs for both 2D and 3D experiments. That is the generality adapted from Geometric Measure Theory. 

.. _MedianShape: https://github.com/tbtraltaa/medianshape

Run source code using conda package management
----------------------------------------------

Linux
-----

We used Ubuntu16.04.

* Install conda
* conda create -name medianshape python=2
* source ~/.bashrc
* source activate medianshape
* conda install scipy=0.17.1
* conda install numpy=1.11.1
* conda install matplotlib=1.5.1
* conda install cvxopt=1.1.8

.. Linux 32-bit

.. sudo apt-get install lib32ncurses5
.. sudo apt-get install lib32z1

Installing PyDistMesh
---------------------
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

* NumPy_
* SciPy_
* PyDistMesh_
* cplex_ or cvxopt_
* matplotlib_ (optional)

Please refer to their own documentation for different version compatability.
You can use Anaconda_ tool to easily install numpy, scipy. 
Anaconda is a collection python packages for scientific computation and it provides conda_
a package manager and an environment manager additionally.

.. _NumPy: http://numpy.org/
.. _SciPy: https://scipy.org/
.. _PyDistMesh: https://pypi.python.org/pypi/PyDistMesh/1.2
.. _cvxopt: http://cvxopt.org
.. _cplex: https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
.. _matplotlib: http://matplotlib.org
.. _Anaconda: https://www.continuum.io
.. _conda: https://conda.io
