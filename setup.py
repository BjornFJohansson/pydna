#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer
#versioneer.VCS = 'git'
#versioneer.versionfile_source = 'pydna/_version.py'
#versioneer.versionfile_build = 'pydna/_version.py'
#versioneer.tag_prefix = '' # tags are like 1.2.0
#versioneer.parentdir_prefix = '' # dirname like 'myproject-1.2.0'

# Read author etc..
for line in open('pydna/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

from setuptools import setup
import textwrap, sys

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
        install_requires =[
        "biopython",
        "networkx",
        "appdirs",
        "prettytable"],
        test_suite="run_tests.load_my_tests",
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
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])

