from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

setup(name="medianshape",
    version='1.0.0',
    description='Median Shape',
    long_description=long_description,
    keywords='Simplicial Median shape',
    author='WSU',
    author_email='',
    url = 'https://github.com',
    license = 'LICENSE.txt',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'cvxopt',
        'numpy',
        'scipy',
        'matplotlib',
        'PyDistMesh',
    ]

)
