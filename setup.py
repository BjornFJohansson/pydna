#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('pydna/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

from setuptools import setup

setup(  name            = 'pydna',
        version=versioneer.get_version()[:5],
        cmdclass=versioneer.get_cmdclass(),
        author          = __author__,
        author_email    = __email__,
        packages=['pydna',
                  'pydna.py_rstr_max',],
        url='http://pypi.python.org/pypi/pydna/',
        license='LICENSE.txt',
        description='''Contains classes and code for representing double
                     stranded DNA and functions for simulating homologous
                     recombination between DNA molecules.''',
        long_description=open('README.rst').read(),
        setup_requires=['pytest-runner', "biopython", "networkx", "appdirs", "prettytable", "ordered-set", "pyparsing"],
        tests_require=['pytest', "biopython", "networkx", "appdirs", "prettytable", "ordered-set", "pyparsing"],
        install_requires = [ "biopython", "networkx", "appdirs", "prettytable", "ordered-set", "pyparsing"],
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Operating System :: OS Independent',
                       'Programming Language :: Python :: 3.5',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
