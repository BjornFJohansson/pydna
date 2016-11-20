|icon| pydna
============

|GitHub stars|\ |GitHub issues|

Planning genetic constructs with many parts and assembly steps, such as
recombinant metabolic pathways are often difficult to properly document.

The Pydna python package provide a human-readable formal description of
cloning and assembly strategies in Python which also allows for
simulation and verification of cloning strategies. Pydna can be though
of as executable documentation for molecular biology.

Pydna provides simulation of:

-  Restriction digestion
-  Ligation
-  PCR
-  Primer design
-  Gibson assembly
-  Homologous recombination
-  Gel electrophoresis of DNA (**NEW Feature**)

Any sub-cloning experiment can be described in pydna, and execution
yield the sequence of the of the resulting DNA molecule and intermediate
steps. A cloning strategy expressed in pydna is *complete*,
*unambiguous* and *stable*. Pydna has been designed to be understandable
for biologists with some basic understanding of Python.

Pydna can formalize planning and sharing of cloning strategies and is
especially useful for complex or combinatorial DNA molecule
constructions.

Look at some assembly strategies made in the Jupyter notebook format
`here <http://nbviewer.ipython.org/github/BjornFJohansson/ypk-xylose-pathways/blob/master/index.ipynb>`__.

There is an open access paper in BMC Bioinformatics describing pydna:

|abstr|

Please reference the above paper when using pydna.

Usage
-----

Most pydna functionality is implemented as methods for the double
stranded DNA sequence record classes Dseq and Dseqrecord, which are
subclasses of the `Biopython <http://biopython.org/wiki/Main_Page>`__
`Seq <http://biopython.org/wiki/Seq>`__ and
`SeqRecord <http://biopython.org/wiki/SeqRecord>`__ classes.

These classes make cut and paste cloning and PCR very simple:

::

    >>> import pydna
    >>> seq = pydna.Dseq("GGATCCAAA","TTTGGATCC",ovhg=0)
    >>> seq
    Dseq(-9)
    GGATCCAAA
    CCTAGGTTT
    >>> from Bio.Restriction import BamHI
    >>> a,b = seq.cut(BamHI)
    >>> a
    Dseq(-5)
    G
    CCTAG
    >>> b
    Dseq(-8)
    GATCCAAA
        GTTT
    >>> a+b
    Dseq(-9)
    GGATCCAAA
    CCTAGGTTT
    >>> b+a
    Dseq(-13)
    GATCCAAAG
        GTTTCCTAG
    >>> b+a+b
    Dseq(-17)
    GATCCAAAGGATCCAAA
        GTTTCCTAGGTTT
    >>> b+a+a
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/dist-packages/pydna/dsdna.py", line 217, in __add__
        raise TypeError("sticky ends not compatible!")
    TypeError: sticky ends not compatible!
    >>>

As the example above shows, pydna keeps track of sticky ends.

Notably, homologous recombination and Gibson assembly between linear DNA
fragments can be easily simulated without any additional information
besides the primary sequence of the fragments.

Gel electrophoresis of DNA fragments can be simulated using the included
gel module by `Bruno Silva <https://github.com/bruno2git>`__:

.. figure:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/gel.png
   :alt: simulated agarose gel

   alt text

Look at an example notebook with a gel simulation
`here <http://nbviewer.jupyter.org/github/BjornFJohansson/pydna/blob/master/scripts/gel_inline_ex.ipynb>`__.

Pydna can be very compact. The nine lines of Python below, simulates the
construction of a recombinant plasmid. DNA sequences are downloaded from
Genbank by accession numbers that are guaranteed to be stable over time.

::

    import pydna
    gb = pydna.Genbank("myself@email.com") # Tell Genbank who you are!
    gene = gb.nucleotide("X06997") # Kluyveromyces lactis LAC12 gene for lactose permease.
    primer_f,primer_r = pydna.parse(''' >760_KlLAC12_rv (20-mer)
                                        ttaaacagattctgcctctg

                                        >759_KlLAC12_fw (19-mer)
                                        aaatggcagatcattcgag
                                        ''', ds=False)
    pcr_prod = pydna.pcr(primer_f,primer_r, gene)
    vector = gb.nucleotide("AJ001614") # pCAPs cloning vector
    from Bio.Restriction import EcoRV
    lin_vector = vector.linearize(EcoRV)
    rec_vec =  ( lin_vector + pcr_prod ).looped()

Pydna can automate the simulation of `sub
cloning <http://en.wikipedia.org/wiki/Subcloning>`__ experiments using
python. This is helpful to generate examples for teaching purposes.

Read the `documentation <http://pydna.readthedocs.io/index.html>`__ or
the
`cookbook <https://www.dropbox.com/sh/4re9a0wk03m95z4/AABpu4zwq4IuKUvK0Iy9Io0Fa?dl=0>`__
with example files for further information.

Please post a message in the `google
group <https://groups.google.com/d/forum/pydna>`__ for pydna if you have
problems, questions or comments.

Feedback is very welcome!

Documentation
-------------

Documentation is built using `Sphinx <http://www.sphinx-doc.org/>`__
from `docstrings <https://www.python.org/dev/peps/pep-0257/>`__ in the
code and displayed at readthedocs |Documentation Status|

The `numpy <www.numpy.org>`__ `docstring
format <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`__
is used.

Installation using conda on Anaconda
------------------------------------

The absolutely best way of installing and using pydna is to use the free
`Anaconda <https://store.continuum.io/cshop/anaconda>`__ or
`Miniconda <http://conda.pydata.org/miniconda.html>`__ python
distributions.

Anaconda is a large download (about 400 Mb) while Miniconda is about
40-50 Mb.

Once Anaconda (or Miniconda) is installed, the conda package manager can
be used to install pydna. Pydna and its dependencies are available from
the `BjornFJohansson <https://anaconda.org/bjornfjohansson>`__ package
channel ast `Anaconda.org <https://anaconda.org>`__.

The first step is to add the channel by typing the command below
followed by return:

::

    conda config --append channels BjornFJohansson

Then pydna can be installed by typing the command below followed by
return:

::

    conda install pydna

This works on Windows, MacOSX and Linux, and installs all necessary and
optional dependencies automatically (see below).

Installation using pip
----------------------

The second best way of installing pydna is with pip, the officially
`recommended <http://python-packaging-user-guide.readthedocs.org/en/latest>`__
tool.

Pip is included in recent Python versions.

Pip installs the minimal installation requirements automatically, but
not the optional requirements (see below). This will probably not work
directly on windows, as biopython is not directly installable.

Linux:
~~~~~~

::

    bjorn@bjorn-UL30A:~/pydna$ sudo pip install pydna

Windows:
~~~~~~~~

Installing biopython on Windows can be tricky. The biopython site has
`executable installers <http://biopython.org/wiki/Download>`__. Read
`here <http://biopython.org/DIST/docs/install/Installation.html>`__ on
how to install biopython requirements such as Numpy. Christoph Gohlke at
University of California, Irvine has compiled many `binary
installers <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`__ for Windows
wich include most requirements.

When the requrements are installed you can pip install pydna from the
Windows terminal:

::

    C:\> pip install pydna

Installation from Source
------------------------

If you install from source, you need to install all dependencies
separately (listed above). Download one of the source installers from
the pypi site or from Github and extract the file. Open the pydna source
code directory (containing the setup.py file) in terminal and type:

::

    python setup.py install

Source Code
-----------

Pydna is developed on
`Github <https://github.com/BjornFJohansson/pydna>`__.

Minimal installation requirements
---------------------------------

Pydna is currently developed on and for Python 3.5. Pydna versions
before 1.0.0 were compatible with python 2.7 only. The list below is the
minimal requirements for installing pydna. Biopython has c-extensions,
but the other modules are pure python.

-  `Python 3.5 <http://www.python.org>`__
-  `biopython >= 1.65 <http://pypi.python.org/pypi/biopython>`__
-  `networkx >= 1.8.1 <http://pypi.python.org/pypi/networkx>`__
-  `pyparsing >= 2.1.10 <https://pypi.python.org/pypi/pyparsing>`__
-  `appdirs >=1.3.0 <https://pypi.python.org/pypi/appdirs>`__
-  `prettytable>=0.7.2 <https://pypi.python.org/pypi/PrettyTable>`__
-  `ordered\_set>=2.0.1 <https://pypi.python.org/pypi/ordered-set>`__

Optional Requirements
---------------------

Pydna has been designed to be used from the Jupyter notebook. If
`IPython <https://ipython.org/>`__ and `Jupyter <http://jupyter.org/>`__
are installed, importing ipython notebooks as modules among are
supported among other things.

If the modules listed below are installed, gel simulation functionality
is available.

-  `numpy <http://www.numpy.org>`__
-  `scipy <https://www.scipy.org>`__
-  `matplotlib <http://matplotlib.org>`__
-  `mpldatacursor <https://pypi.python.org/pypi/mpldatacursor>`__
-  `pint >= 0.7.2 <https://pypi.python.org/pypi/pint>`__

The pydna conda package installs the optional requirements listed above
as well as:

-  `ipython <https://pypi.python.org/pypi/ipython>`__
-  `jupyter <https://pypi.python.org/pypi/jupyter>`__

Requirements for running tests
------------------------------

-  `pytest>=3.0.3 <https://pypi.python.org/pypi/pytest>`__

Requirements for analyzing code coverage
----------------------------------------

-  `python-coveralls >=
   2.9.0 <https://pypi.python.org/pypi/python-coveralls>`__
-  `coverage >= 3.7.1 <https://pypi.python.org/pypi/coverage>`__
-  `pytest-cov >= 2.3.1 <https://pypi.python.org/pypi/pytest-cov>`__

Automatic testing
-----------------

The test suit is run automatically after each commit on:

-  Ubuntu 12.04 using drone.io |icon3|
-  Ubuntu 14.04 using CircleCI |CircleCI|
-  OSX-64 using TravisCI |icon1|
-  Windows using AppveyorCI |icon2|.

Code coverage is |Coverage Status|

Dependencies are monitored by versioneye |icon11|

Automatic builds
----------------

`Conda <http://conda.pydata.org/docs/intro.html>`__ packages are built
on CircleCI(Linux), Drone.io(Linux), TravisCI(MacOS) and
AppveyorCI(Windows). Source setuptools packages and wheels are built on
Linux for all systems. Binary setuptools packages are built for Windows
and MacOSX.

-  Conda packages |Anaconda-Server Badge|
-  Setuptools packages

Builds are controlled by Git tags. Tags like 1.0.2a4 are considered test
builds and are uploaded to
`testpypi <https://testpypi.python.org/pypi?:action=display&name=pydna>`__
and to Anaconda.org with a "test" label. These are only meant to test
the finished packages and are not meant to be used.

Tags like 1.0.3 are considered final builds and are built and uploaded
to `Anaconda.org <https://anaconda.org/BjornFJohansson/pydna>`__ under
the "main" label and to the regular
`pypi <https://pypi.python.org/pypi/pydna>`__ server.

Changelog
---------

See the `change
log <https://raw.githubusercontent.com/BjornFJohansson/pydna/py3/CHANGELOG.md>`__
for recent changes.

.. |icon| image:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/pydna.resized.png
   :target: https://pypi.python.org/pypi/pydna/
.. |GitHub stars| image:: https://img.shields.io/github/stars/BjornFJohansson/pydna.svg
   :target: https://github.com/BjornFJohansson/pydna/stargazers
.. |GitHub issues| image:: https://img.shields.io/github/issues/BjornFJohansson/pydna.svg
   :target: https://github.com/BjornFJohansson/pydna/issues
.. |abstr| image:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/BMC_resized.png
   :target: http://www.biomedcentral.com/1471-2105/16/142/abstract
.. |Documentation Status| image:: https://readthedocs.org/projects/pydna/badge/?version=latest
   :target: http://pydna.readthedocs.io/?badge=latest
.. |icon3| image:: https://drone.io/github.com/BjornFJohansson/pydna/status.png
   :target: https://drone.io/github.com/BjornFJohansson/pydna/latest
.. |CircleCI| image:: https://circleci.com/gh/BjornFJohansson/pydna/tree/py3dev.svg?style=shield
   :target: https://circleci.com/gh/BjornFJohansson/pydna/tree/py3dev
.. |icon1| image:: https://travis-ci.org/BjornFJohansson/pydna.svg
   :target: https://travis-ci.org/BjornFJohansson/pydna
.. |icon2| image:: https://ci.appveyor.com/api/projects/status/qdtk9biw5o0cae7u?svg=true
   :target: https://ci.appveyor.com/project/BjornFJohansson/pydna
.. |Coverage Status| image:: https://coveralls.io/repos/github/BjornFJohansson/pydna/badge.svg?branch=py3
   :target: https://coveralls.io/github/BjornFJohansson/pydna?branch=py3
.. |icon11| image:: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb/badge.svg
   :target: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb
.. |Anaconda-Server Badge| image:: https://anaconda.org/bjornfjohansson/pydna/badges/version.svg
   :target: https://anaconda.org/bjornfjohansson/pydna
