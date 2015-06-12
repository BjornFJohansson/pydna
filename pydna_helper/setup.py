#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

__version__ = "Undefined"
for line in open('pydna_helper/__init__.py'):
    if line.startswith('__'):
        exec(line.strip())

setup( name='pydna_helper',
       description='pydna utilities',
       version         =__version__,
       author          =__author__,
       author_email    =__email__,
       packages=['pydna_helper'],
       zip_safe = False,
       install_requires =[ "pydna>=0.6.1", "percache>=0.2.1", "appdirs>=1.2.0"],
       keywords = "bioinformatics",
       classifiers = ['Development Status :: 3 - Alpha',
                      'Environment :: Console',
                      'Intended Audience :: Education',
                      'Intended Audience :: Science/Research',
                      'License :: OSI Approved :: BSD License',
                      'Programming Language :: Python :: 2.7',
                      'Topic :: Education',
                      'Topic :: Scientific/Engineering :: Bio-Informatics',])
