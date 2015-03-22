from setuptools import setup, find_packages

with open('README.rst') as f:
    long_description = f.read()

setup(name="MedianShape"
    version = '0.1.1',
    description = 'Median Shape',
    keywords='',
    author='Team',
    auther_email='',
    url = 'https://github.com',
    test_suite='',
    include_package_data=True,
    packages=find_packages(),
)
