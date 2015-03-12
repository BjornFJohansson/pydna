#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'pydna/_version.py'
versioneer.versionfile_build = 'pydna/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = '' # dirname like 'myproject-1.2.0'

# Read author etc..
for line in open('pydna/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

# Change line ending to windows for all text files
import os
for root, dirs, files in os.walk(os.path.abspath(os.path.dirname(__file__))):
    for name in files:
        if not name.lower().endswith(".txt"):
            continue
        filename = os.path.join(root, name)
        with open(filename, "rb") as f:
            data = f.read()
        temp = data.replace('\r\n', '\n')
        temp = temp.replace('\r',   '\n')
        temp = temp.replace('\n',   '\r\n')
        if not data == temp:
            with open(filename, "wb") as f:
                f.write(temp)
                print("changed", filename)

from setuptools import setup #, find_packages

import textwrap, sys, nose

setup(  name='pydna',
        author          =__author__,
        author_email    =__email__,
        version=versioneer.get_version()[:5],
        cmdclass=versioneer.get_cmdclass(),        
        #cmdclass = {'build_ext': optional_build_ext},
        #packages = find_packages(),
        packages=['pydna',
                  'pydna.py_rstr_max',],

        #ext_modules=[Extension('pydna.find_sub_strings', ['pydna/findsubstrings_numpy_arrays_cython.c'])],
        #scripts=[],
        url='http://pypi.python.org/pypi/pydna/',
        license='LICENSE.txt',
        description='''Contains classes and code for representing double
                     stranded DNA and functions for simulating homologous
                     recombination between DNA molecules.''',
        long_description=open('README.rst').read(),
        install_requires =[
        "biopython",
        "networkx",
        "appdirs",
        "prettytable",],
        #test_suite = 'nose.collector',
        test_suite="run_tests.load_my_tests",
        #include_dirs = [numpy.get_include()],
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 2.7',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',]
        )

