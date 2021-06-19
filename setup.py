#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='ndpolator',
      version='0.1.0',
      description='Multi-dimensional linear interpolation, extrapolation and imputation',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Andrej Prsa',
      author_email='aprsa@villanova.edu',
      url='https://www.github.com/aprsa/ndpolator',
      download_url = 'https://github.com/aprsa/ndpolator/tarball/0.1.0',
      packages=['ndpolator'],
      install_requires=[],
      classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development :: Libraries :: Python Modules",
       ],
     )
