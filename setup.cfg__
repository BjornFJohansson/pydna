#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

from setuptools import setup
from setuptools import find_packages
from os import path
import re

# Read __author__, __email__. from __init__.py
__author__ = "__author__"
__email__ = "__email__"
for line in open("src/pydna/__init__.py"):
    if line.startswith("__author__") or line.startswith("__email__"):
        exec(line.strip())

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


install_requires = []
with open("requirements.txt", encoding="utf-8") as f:
    for line in f.readlines():
        install_requires.append(re.split(r"(<|>|=)=", line)[0])

setup(name="pydna",
      author=__author__,
      author_email=__email__,
      zip_safe=False,
      packages=find_packages("src"),
      package_dir={"": "src"},
      url="https://github.com/BjornFJohansson/pydna",
      license="LICENSE.txt",
      description=("Representing double stranded DNA and functions for "
                   "simulating cloning and homologous recombination"
                   " between DNA molecules."),
      long_description=long_description,
      long_description_content_type="text/markdown",
      setup_requires=["setuptools_scm"],
      use_scm_version={"write_to": "src/pydna/_version.py"},
      install_requires=install_requires,
      keywords="bioinformatics",
      classifiers=["Development Status :: 4 - Beta",
                   "Environment :: Console",
                   "Intended Audience :: Education",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: BSD License",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python :: 3.7",
                   "Programming Language :: Python :: 3.8",
                   "Programming Language :: Python :: 3.9",
                   "Programming Language :: Python :: 3.10",
                   "Programming Language :: Python :: 3.11",
                   "Topic :: Education",
                   "Topic :: Scientific/Engineering :: Bio-Informatics", ],)
