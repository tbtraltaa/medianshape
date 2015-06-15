from setuptools import setup, find_packages

with open('README.rst') as f:
    long_description = f.read()

setup(name="MedianShape"
    version = '0.1.1',
    description = 'Median Shape',
    long_description = readme(),
    classifiers=[
        'Development Status :: 0.1',
        'License ::  MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Median Shapes :: Shape Statistics',
    ],
    keywords='',
    author='WSU',
    author_email='',
    url = 'https://github.com',
    license = ''
    test_suite='',
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
