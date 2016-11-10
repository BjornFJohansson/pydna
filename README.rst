|icon| pydna
============

Planning genetic constructs with many parts, such as recombinant
metabolic pathways is usually done manually using a DNA sequence editor,
a task which quickly becomes unfeasible as scale and complexity of the
constructions increase.

The Pydna python package provide a human-readable formal description of
cloning and assembly strategies which also allows for automatic computer
simulation and verification.

Pydna provides simulation of:

-  Restriction digestion
-  Ligation
-  PCR
-  Primer design
-  Gibson assembly
-  Homologous recombination
-  Gel electrophoresis of DNA (**NEW Feature**)

Pydna was designed to provide a form of executable documentation
describing a subcloning or DNA assembly experiment. The pydna code
unambiguously describe a sub cloning experiment, and can be executed to
yield the sequence of the of the resulting DNA molecule. A cloning
strategy expressed in pydna is complete, unambiguous and stable. Pydna
has been designed to be understandable for biologists with some basic
understanding of Python.

Pydna can formalize planning and sharing of cloning strategies and is
especially useful for complex or combinatorial DNA molecule
constructions.

Look at some assembly strategies made in the Jupyter notebook format
`here <http://nbviewer.ipython.org/github/BjornFJohansson/ypk-xylose-pathways/blob/master/index.ipynb>`__.

There is also an open access paper in BMC Bioinformatics describing
pydna:

|abstr|

Please make a reference to the above paper if you publish work where
pydna was used.

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

Notably, homologous recombination and Gibson assembly between linear DNA
fragments can be easily simulated without any additional information
besides the primary sequence of the fragments.

Gel electrophoresis of DNA fragments can be simulated using the gel.py
module by `Bruno Silva <https://github.com/bruno2git>`__:

.. figure:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/gel.png
   :alt: simulated agarose gel

   alt text

Look at an example notebook with a gel simulation
`here <http://nbviewer.jupyter.org/github/BjornFJohansson/pydna/blob/master/scripts/gel_inline_ex.ipynb>`__.

The nine lines of Python below, simulates the construction of a
recombinant plasmid. DNA sequences are downloaded from Genbank by
accession numbers that are guaranteed to be stable over time.

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

Pydna is also be useful to automate the simulation of `sub
cloning <http://en.wikipedia.org/wiki/Subcloning>`__ experiments using
python. This is helpful to generate examples for teaching purposes.

Read the `documentation <http://pydna.readthedocs.io/index.html>`__ or
the
`cookbook <https://www.dropbox.com/sh/4re9a0wk03m95z4/AABpu4zwq4IuKUvK0Iy9Io0Fa?dl=0>`__
with example files for further information.

Please post a message in the `google
group <https://groups.google.com/d/forum/pydna>`__ for pydna if you have
problems, questions or comments. Feedback in the form of questions,
comments or criticism is very welcome!

Automatic testing and builds
----------------------------

The test suit is run automatically after each commit on:

-  OSX-64 using travis |icon1|
-  Windows using appveyor |icon2|.
-  Ubuntu using drone |icon3|

Documentation is built and displayed at readthedocs |Documentation
Status|

Code coverage is |Coverage Status|

Dependencies are monitored by versioneye |icon11|

Github stars: |GitHub stars|

Github issues: |GitHub issues|

Minimal installation requirements
---------------------------------

Pydna is currently developed on and for Python 3.5. Pydna versions
before 1.0.0 were compatible with python 2.7 only. The list below is the
minimal requirements for installing pydna. Biopython has c-extensions,
but the other modules are pure python.

-  `Python 3 <http://www.python.org>`__
-  `biopython >= 1.65 <http://pypi.python.org/pypi/biopython>`__
-  `networkx >= 1.8.1 <http://pypi.python.org/pypi/networkx>`__
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

-  `pint >=0.7.2 <https://pypi.python.org/pypi/pint>`__
-  `scipy <https://www.scipy.org>`__
-  `numpy <http://www.numpy.org>`__
-  `matplotlib <http://matplotlib.org>`__
-  `mpldatacursor <https://pypi.python.org/pypi/mpldatacursor>`__

The pydna conda package installs all optional requirements, currently:

-  `ipython>=4 <https://pypi.python.org/pypi/ipython>`__
-  `jupyter>=1.0.0 <https://pypi.python.org/pypi/jupyter>`__
-  `scipy>=0.16.0 <https://pypi.python.org/pypi/scipy>`__
-  `numpy>=1.10.1 <https://pypi.python.org/pypi/numpy>`__
-  `matplotlib>=1.5.0 <https://pypi.python.org/pypi/matplotlib>`__
-  `mpldatacursor>=0.6.1 <https://pypi.python.org/pypi/mpldatacursor>`__
-  `pint >=0.7.2 <https://pypi.python.org/pypi/pint>`__

Requirements for running tests
------------------------------

-  `nose>=1.3.4 <https://pypi.python.org/pypi/nose>`__
-  `coverage>=3.7.1 <https://pypi.python.org/pypi/coverage>`__

Installation using conda on Anaconda
------------------------------------

The absolutely best way of installing and using pydna is to use the free
`Anaconda <https://store.continuum.io/cshop/anaconda>`__ python
distribution.

Once Anaconda is installed, the conda package manager can be used to
install pydna. Pydna and its dependencies are available from the
`conda-forge <https://anaconda.org/conda-forge>`__ and
`BjornFJohansson <https://anaconda.org/bjornfjohansson>`__
`Anaconda.org <https://anaconda.org>`__ channels. The first step is to
add the channels:

::

    conda config --append channels conda-forge
    conda config --append channels BjornFJohansson

Then pydna can be installed by simply:

::

    conda install pydna

This works on Windows, MacOSX and Linux, and installs all necessary and
optional dependencies automatically in one go.

Installation using pip
----------------------

The second best way of installing pydna is with pip, the officially
`recommended <http://python-packaging-user-guide.readthedocs.org/en/latest>`__
tool.

Pip installs the minimal installation requirements automatically, but
not the optional requirements (see above). These have to be installed
manually.

Linux:
~~~~~~

::

    bjorn@bjorn-UL30A:~/pydna$ sudo pip install pydna

Windows:
~~~~~~~~

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

Windows dependencies
~~~~~~~~~~~~~~~~~~~~

Sometimes dependencies can be difficult to install on windows, as a C
compiler is necessary. If dependencies have to be installed separately,
this can be done using these binary installers for Windows:

+-----------------+----------------------------------------------------------+
| Dependency      | link                                                     |
+=================+==========================================================+
| Python (32,64)  | http://www.python.org/download                           |
+-----------------+----------------------------------------------------------+
| Biopython (32)  | http://biopython.org/wiki/Download                       |
+-----------------+----------------------------------------------------------+
| Biopython (64)  | http://www.lfd.uci.edu/~gohlke/pythonlibs/#biopython     |
+-----------------+----------------------------------------------------------+
| numpy (32,64)   | http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy         |
+-----------------+----------------------------------------------------------+
| networkx        | http://www.lfd.uci.edu/~gohlke/pythonlibs/#networkx      |
| (32,64)         |                                                          |
+-----------------+----------------------------------------------------------+
| pint            | http://www.lfd.uci.edu/~gohlke/pythonlibs/Pint-0.6-py2.p |
|                 | y3-none-any.whl                                          |
+-----------------+----------------------------------------------------------+
| scipy (32,64)   | http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy         |
+-----------------+----------------------------------------------------------+
| matplotlib      | http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib    |
| (32,64)         |                                                          |
+-----------------+----------------------------------------------------------+
| ipython>=4.0    | http://www.lfd.uci.edu/~gohlke/pythonlibs/#ipython       |
+-----------------+----------------------------------------------------------+
| jupyter         | http://www.lfd.uci.edu/~gohlke/pythonlibs/#jupyter       |
+-----------------+----------------------------------------------------------+

Source Code
-----------

Pydna is developed on
`Github <https://github.com/BjornFJohansson/pydna>`__.

Changelog
---------

See the `change
log <https://raw.githubusercontent.com/BjornFJohansson/pydna/py3/CHANGELOG.md>`__
for recent changes.

.. |icon| image:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/pydna.resized.png
   :target: https://pypi.python.org/pypi/pydna/
.. |abstr| image:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/BMC_resized.png
   :target: http://www.biomedcentral.com/1471-2105/16/142/abstract
.. |icon1| image:: https://travis-ci.org/BjornFJohansson/pydna.svg
   :target: https://travis-ci.org/BjornFJohansson/pydna
.. |icon2| image:: https://ci.appveyor.com/api/projects/status/qdtk9biw5o0cae7u?svg=true
   :target: https://ci.appveyor.com/project/BjornFJohansson/pydna
.. |icon3| image:: https://drone.io/github.com/BjornFJohansson/pydna/status.png
   :target: https://drone.io/github.com/BjornFJohansson/pydna/latest
.. |Documentation Status| image:: https://readthedocs.org/projects/pydna/badge/?version=latest
   :target: http://pydna.readthedocs.io/?badge=latest
.. |Coverage Status| image:: https://coveralls.io/repos/github/BjornFJohansson/pydna/badge.svg?branch=py3
   :target: https://coveralls.io/github/BjornFJohansson/pydna?branch=py3
.. |icon11| image:: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb/badge.svg
   :target: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb
.. |GitHub stars| image:: https://img.shields.io/github/stars/BjornFJohansson/pydna.svg
   :target: https://github.com/BjornFJohansson/pydna/stargazers
.. |GitHub issues| image:: https://img.shields.io/github/issues/BjornFJohansson/pydna.svg
   :target: https://github.com/BjornFJohansson/pydna/issues
