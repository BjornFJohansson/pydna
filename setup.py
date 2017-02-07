#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('pydna/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

from setuptools import setup
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()

setup(  name            = 'pydna',
        version=versioneer.get_version()[:5],
        cmdclass=versioneer.get_cmdclass(),
        author          = __author__,
        author_email    = __email__,
        zip_safe = False,
        packages=['pydna',
                  'pydna._py_rstr_max',],
        url='http://pypi.python.org/pypi/pydna/',
        license='LICENSE.txt',
        description='''Contains classes and code for representing double
                     stranded DNA and functions for simulating homologous
                     recombination between DNA molecules.''',
        long_description=read_md('README.md'),
        setup_requires=['pytest-runner', "appdirs", "biopython", "prettytable",  "networkx", "ordered-set", "pyparsing", "requests"],
        tests_require=['pytest',         "appdirs", "biopython", "prettytable",  "networkx", "ordered-set", "pyparsing", "requests"],
        install_requires = [             "appdirs", "biopython", "prettytable",  "networkx", "ordered-set", "pyparsing", "requests"],
        keywords = "bioinformatics",
        platforms=['CPython 3.5'],
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
