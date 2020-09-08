#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Read __author__, __email__. from __init__.py
__author__ = "__author__"
__email__  = "__email__"
for line in open("pydna/__init__.py"):
    if line.startswith("__author__") or line.startswith("__email__"):
        exec(line.strip())

from setuptools import setup

from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pydna",
    author=__author__,
    author_email=__email__,
    zip_safe=False,
    packages=["pydna"],
    url="https://github.com/BjornFJohansson/pydna",
    license="LICENSE.txt",
    description="""Contains classes and code for representing double
                     stranded DNA and functions for simulating homologous
                     recombination between DNA molecules.""",
    long_description=long_description,
    long_description_content_type="text/markdown",
    setup_requires =["pytest-runner", "setuptools_scm"],
    tests_require=["pytest"],
    use_scm_version={"write_to": "pydna/_version.py"},
    install_requires=[
        "appdirs",
        "biopython",
        "networkx",
        "prettytable",
        "pyparsing",
        "requests",
    ],
    keywords="bioinformatics",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Education",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
