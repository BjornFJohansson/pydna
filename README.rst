|icon| pydna
============

|icon1|

|icon2|

|icon3|

|icon4|

|icon5|

|icon6|

|icon7|

|icon8|

|icon9|

|icon10|

|icon11|

Planning genetic constructs with many parts, such as recombinant
metabolic pathways is usually done manually using a DNA sequence editor,
a task which quickly becomes unfeasible as scale and complexity of the
constructions increase.

The Pydna python package provide a human-readable formal description of
cloning and assembly strategies which also allows for automatic computer
simulation and verification.

Pydna provides simulation of:

-  restriction digestion
-  ligation
-  PCR
-  primer design
-  Gibson assembly
-  homologous recombination

A cloning strategy expressed in pydna is complete, unambiguous and
stable. Pydna has been designed to be understandable for biologists with
limited programming skills.

Pydna formalize planning and sharing of cloning strategies and is
especially useful for complex or combinatorial DNA molecule
constructions.

Look at some assembly strategies made in the IPython notebook format
`here <http://nbviewer.ipython.org/github/BjornFJohansson/ypk-xylose-pathways/blob/master/index.ipynb>`__.

There at the open access BMC Bioinformatics publication describing
pydna:

|abstr|

Double stranded DNA sequence classes that make cut and paste cloning and
PCR very simple is provided.

See an example of pydna usage at the command line below:

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
other than the primary sequence of the fragments.

Most pydna functionality is implemented as methods for the double
stranded DNA sequence record classes Dseq and Dseqrecord, which are
subclasses of the `Biopython <http://biopython.org/wiki/Main_Page>`__
`Seq <http://biopython.org/wiki/Seq>`__ and
`SeqRecord <http://biopython.org/wiki/SeqRecord>`__ classes.

Pydna was designed to provide a form of executable documentation
describing a subcloning or DNA assembly experiment. The pydna code
unambiguously describe a sub cloning experiment, and can be executed to
yield the sequence of the of the resulting DNA molecule.

Pydna was designed to semantically imitate how sub cloning experiments
are typically documented in Scientific literature. Pydna code describing
a sub cloning is reasonably compact and meant to be easily readable.

The nine lines of Python below, simulates the construction of a
recombinant plasmid. DNA sequences are downloaded from Genbank by
accession numbers that are guaranteed to be stable.

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

Pydna might also be useful to automate the simulation of `sub
cloning <http://en.wikipedia.org/wiki/Subcloning>`__ experiments using
python. This could be helpful to generate examples for teaching
purposes. Read the
`documentation <https://pydna.readthedocs.org/en/latest>`__ or the
`cookbook <https://www.dropbox.com/sh/4re9a0wk03m95z4/AABpu4zwq4IuKUvK0Iy9Io0Fa?dl=0>`__
with example files for further information.

An `on-line <http://pydna-shell.appspot.com>`__ shell running Python
with pydna is available for simple experimentation. It is slower than
rinning pydna on your own computer.

Please post a message in the `google
group <https://groups.google.com/d/forum/pydna>`__ for pydna if you have
problems, questions or comments. Feedback in the form of questions,
comments or criticism is very welcome! ## Installation requirements

This package was developed on and for Python 2.7. Other versions have
not been tested.

-  `Python 2.7 <http://www.python.org>`__
-  `biopython >= 1.65 <http://pypi.python.org/pypi/biopython>`__
-  `networkx >= 1.8.1 <http://pypi.python.org/pypi/networkx>`__
-  `appdirs >=1.3.0 <https://pypi.python.org/pypi/appdir>`__
-  `prettytable>=0.7.2 <https://pypi.python.org/pypi/PrettyTable>`__

Requirements for running tests
------------------------------

-  `nose>=1.3.4 <https://pypi.python.org/pypi/nose>`__
-  `coverage>=3.7.1 <https://pypi.python.org/pypi/coverage>`__

Optional Requirements
---------------------

-  `ipython<=3.2.1 <https://pypi.python.org/pypi/ipython>`__

Pydna has been designed to be used from the IPython notebook. If you
have IPython installed, there are functions in pydna for importing
ipython notebooks as modules among other things.

Python 3
--------

This code has not been tried with Python 3. If there is sufficient
interest, there might be a Python 3 version in the future.

Installation using conda on Anaconda
------------------------------------

The best way of using Python in general is to use a free distribution
such as `Anaconda <https://store.continuum.io/cshop/anaconda>`__

There is a `conda <https://anaconda.org/bjornfjohansson/pydna>`__
package available for pydna, which is easily installed at the command
line using the conda package manager.

::

    conda install -c https://conda.anaconda.org/bjornfjohansson pydna

This works on Windows, MacOSX and Linux, and installs all dependencies
automatically in one go.

Installation using pip
----------------------

The second best way of installing pydna is with pip. Pip is the
officially
`recommended <http://python-packaging-user-guide.readthedocs.org/en/latest>`__
tool for installation of Python packages from PyPi. Pip installs
dependencies automatically.

Linux:
~~~~~~

::

    bjorn@bjorn-UL30A:~/Dropbox/pydna$ sudo pip install pydna

Windows:
~~~~~~~~

::

    C:\> pip install pydna

If you do not have pip, you can get it by following these
`instructions <http://www.pip-installer.org/en/latest/installing.html>`__

Installation from Source
------------------------

If you install from source, you need to install the dependencies
separately (listed above). Download one of the source installers from
the pypi site and extract the file. Open the pydna source code directory
(containing the setup.py file) in terminal and type:

::

    python setup.py install

Installation from binary distributions
--------------------------------------

There is a 64 bit windows executable and a windows wheel
`here <https://ci.appveyor.com/project/BjornFJohansson/pydna/build/artifacts>`__.
Note that these will not install required dependencies (see below).

Windows dependencies
~~~~~~~~~~~~~~~~~~~~

Sometimes the dependecies can be difficult to install on windows,
especially Biopython as a C compiler is necessary. If dependencies have
to be installed separately, this can be done using the binary installers
for Windows:

+--------------------+--------------------------------------------------------+
| Dependency         | link                                                   |
+====================+========================================================+
| Python (32,64)     | http://www.python.org/download                         |
+--------------------+--------------------------------------------------------+
| Biopython (32)     | http://biopython.org/wiki/Download                     |
+--------------------+--------------------------------------------------------+
| Biopython (64)     | http://www.lfd.uci.edu/~gohlke/pythonlibs/#biopython   |
+--------------------+--------------------------------------------------------+
| networkx (32,64)   | http://www.lfd.uci.edu/~gohlke/pythonlibs/#networkx    |
+--------------------+--------------------------------------------------------+

Source Code Repository
----------------------

Pydna is developed on
`Github <https://github.com/BjornFJohansson/pydna>`__

TODO
----

-  [ ] IPython 4 (Jupyter) support
-  [ ] Add agarose gel electrophoresis simulation

.. |icon| image:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/pydna.resized.png
   :target: https://pypi.python.org/pypi/pydna/
.. |icon1| image:: https://travis-ci.org/BjornFJohansson/pydna.svg
   :target: https://travis-ci.org/BjornFJohansson/pydna
.. |icon2| image:: https://ci.appveyor.com/api/projects/status/qdtk9biw5o0cae7u?svg=true
   :target: https://ci.appveyor.com/project/BjornFJohansson/pydna
.. |icon3| image:: https://drone.io/github.com/BjornFJohansson/pydna/status.png
   :target: https://drone.io/github.com/BjornFJohansson/pydna/latest
.. |icon4| image:: https://binstar.org/bjornfjohansson/pydna/badges/build.svg
   :target: https://binstar.org/bjornfjohansson/pydna/builds
.. |icon5| image:: https://binstar.org/bjornfjohansson/pydna/badges/version.svg
   :target: https://binstar.org/bjornfjohansson/pydna
.. |icon6| image:: https://coveralls.io/repos/BjornFJohansson/pydna/badge.svg?branch=master
   :target: https://coveralls.io/r/BjornFJohansson/pydna?branch=master
.. |icon7| image:: https://readthedocs.org/projects/pydna/badge/?version=latest
   :target: https://readthedocs.org/projects/pydna/?badge=latest
.. |icon8| image:: https://img.shields.io/pypi/v/pydna.png
   :target: https://pypi.python.org/pypi/pydna
.. |icon9| image:: https://img.shields.io/github/stars/BjornFJohansson/pydna.svg
   :target: https://github.com/BjornFJohansson/pydna/stargazers
.. |icon10| image:: https://img.shields.io/pypi/dm/pydna.png
   :target: https://pypi.python.org/pypi/pydna
.. |icon11| image:: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb/badge.svg
   :target: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb
.. |abstr| image:: https://raw.githubusercontent.com/BjornFJohansson/pydna/master/BMC_resized.png
   :target: http://www.biomedcentral.com/1471-2105/16/142/abstract
