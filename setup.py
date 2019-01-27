#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('pydna/__init__.py', encoding="utf-8"):
    if line.startswith('__') and not line.startswith(('__version', '__long')):
        exec(line.strip())

from setuptools import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(  name            = 'pydna',
        version=versioneer.get_version().split("+", 1)[0],
        product_version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        author          = __author__,
        author_email    = __email__,
        zip_safe = False,
        packages=['pydna'],
        url='http://pypi.python.org/pypi/pydna/',
        license='LICENSE.txt',
        description='''Contains classes and code for representing double
                     stranded DNA and functions for simulating homologous
                     recombination between DNA molecules.''',
        long_description=long_description,
        long_description_content_type='text/markdown',
        install_requires = ["appdirs", "biopython", "networkx", "prettytable", "pyparsing", "requests"],
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Operating System :: OS Independent',
                       'Programming Language :: Python :: 3.5',
                       'Programming Language :: Python :: 3.6',
                       'Programming Language :: Python :: 3.7',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
