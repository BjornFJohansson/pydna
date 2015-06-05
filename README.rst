=====
pydna
=====

.. image:: https://travis-ci.org/BjornFJohansson/pydna.svg 
    :target: https://travis-ci.org/BjornFJohansson/pydna

.. image:: https://ci.appveyor.com/api/projects/status/qdtk9biw5o0cae7u?svg=true
    :target: https://ci.appveyor.com/project/BjornFJohansson/pydna

.. image:: https://coveralls.io/repos/BjornFJohansson/pydna/badge.svg?branch=master  
    :target: https://coveralls.io/r/BjornFJohansson/pydna?branch=master 
  
.. image:: https://readthedocs.org/projects/pydna/badge/?version=latest
    :target: https://readthedocs.org/projects/pydna/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/pydna.png
    :target: https://pypi.python.org/pypi/pydna/
    :alt: Downloads

.. image:: https://img.shields.io/github/stars/BjornFJohansson/pydna.svg
    :target: https://github.com/BjornFJohansson/pydna
    :alt: Github
    
.. image:: https://img.shields.io/pypi/dm/pydna.png
    :target: https://pypi.python.org/pypi/pydna/
    :alt: Latest Version

.. image:: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb/badge.svg?style=flat(Dependency Status)!
    :target: https://www.versioneye.com/user/projects/553174c010e714f9e50010bb
    :alt: versioneye

Pydna provide functions for molecular biology using python.
Double stranded DNA sequence classes that make cut and paste
cloning and PCR very simple is provided (see example below). 
`Publication <http://www.biomedcentral.com/1471-2105/16/142/abstract>`_ describing pydna.

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

Notably, homologous recombination and Gibson assembly between linear
DNA fragments can be easily simulated.

Most functionality is implemented as methods for the double stranded
DNA sequence record classes Dseq and Dseqrecord, which are subclasses
of the `Biopython <http://biopython.org/wiki/Main_Page>`_
`Seq <http://biopython.org/wiki/Seq>`_
and
`SeqRecord <http://biopython.org/wiki/SeqRecord>`_ classes.

Pydna was designed to provide a form of executable documentation
describing a subcloning or DNA assembly experiment. The pydna code
unambiguously describe a sub cloning experiment, and can be executed
to yield the sequence of the of the resulting DNA molecule.

Pydna was designed to semantically imitate how sub cloning experiments are
typically documented in Scientific literature. Pydna code describing a
sub cloning is reasonably compact and meant to be easily readable.

The nine lines of Python below, simulates the construction of a recombinant
plasmid. DNA sequences are downloaded from Genbank by accession numbers that
are guaranteed to be stable.

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


Pydna might also be useful to automate the simulation of
`sub cloning <http://en.wikipedia.org/wiki/Subcloning>`_ experiments using
python. This could be helpful to generate examples for teaching purposes. Read
the `documentation <https://pydna.readthedocs.org/en/latest/>`_ or the
`cookbook <https://www.dropbox.com/sh/4re9a0wk03m95z4/AABpu4zwq4IuKUvK0Iy9Io0Fa?dl=0>`_ with example files
for further information.

An `on-line <http://pydna-shell.appspot.com/>`_ shell running Python with
pydna is available for experimentation.

Please post a message in the `google group <https://groups.google.com/d/forum/pydna>`_
for pydna if you have problems, questions or comments.

Feedback in the form of questions, comments or criticism is very welcome!

=======   ========== =====================================================================
version   date       comment
=======   ========== =====================================================================
0.9.3     2015-06-03 Shelve does not work under MacOS under certain conditions. 
                     This release tries to solve this by not specifying file extensions
                     for the cache files. Two functions are added, pydna. 

0.9.2     2015-05-28 pydna_data_dir is encoded to a string in __init__.py instead of 
                     unicode. The Popen module does not accept environment variables that 
                     are not strings.

0.9.1     2015-05-26 fixed critical error in the calculation of seguid and cseguid 
                     checksums

0.9.0     2015-05-26 seguid and cseguid are now url safe so they can be part of urls and
                     file names.
                     Dseqrecord.locus is an alias of Dseqrecord.name
                     Dseqrecord.accession is an alias of Dseqrecord.id
                     Dseqrecord.definition is an alias of Dseqrecord.description					 
                     changed how circular assembly products are identified to use cseguid.
                     removed proxy handling when proxy not set in download module.
                     added CHANGELOG.md, currently empty.
                     environment variable datadir is now pydna_data_dir.
                     removed environmental variable pydna_dna_dir.
                     if Dseqrecord is initiated with a name property that is longer than 
                     16 characters, it is truncated to 16 chars and a warning is issued. 
                     Default Dseqrecord name property is "na".
                     Default Dseqrecord id property is "-".
                     Default Dseqrecord description property is "@".
                     Dseqrecord __eq__ and __ne__ methods defined.
                     Dseqrecord.write now overwrites an old sequence with the same 
                     filename if the primary sequence is the same.
                     Dseqrecord.read now only looks in current working directory.
                     fixed ipynb_import test code.
                    
0.8.4     2015-04-17 Bugfix for parsing text files with unicode characters.

0.8.3     ?          ?   

0.8.2     ?          ?

0.8.1     2015-03-07 Bugfix for windows. The data directory was not created.

0.8.0	  2015-02-06 Mapping reads added.

0.7.2	  2014-11-21 First public release with the changes from 0.7.0 and 0.7.1.
					 Added a Pretty_str class to beautify output of strings in
					 the IPython shell. 

0.7.1     not public Short linkers can be incorporated in PCR primers in the 
                     assembly_primers function.

0.7.0     not public Caching to speed up Amplify, Assembly, download and the 
                     Desqrecord synced method. The data is stored in four shelf
                     files in the users application directory.
                     
                     amplify.shelf
                     assembly.shelf
                     genbank.shelf
                     synced.shelf                     
                     
                     The location is os specific.
                     See the documentation of appdirs 
                     https://pypi.python.org/pypi/appdirs/1.4.0

0.6.6                new function nopcr.

0.6.5     2014-07-31 bugfix: cutting an amplicon object now preserves features 
                     Changed requirement for NetworkX to 1.8.1

0.6.4     2014-07-09 The pcr function and Anneal class can now deal with primers 
                     with ambiguous codons like R = A or G. In the resulting PCR
                     product, the ambiguous nucleotides are preserved in the tails
                     i.e. the primer part not annealing. The annealing part will 
                     have the sequence corresponding to the template.  

0.6.3     2014-07-06 Dseqrecord.add_feature can now take a string or some other
                     sequence as input. The assembly primers function can now produce 
                     primers for a circular assembly.

0.6.2     2014-06-13 Dseqrecord gained three new methods:

                     isorf() method returning True or False.

                     List_features() method returns a list of all features as a
                     formatted ASCII table.

                     Extract_feature() extracts a feature in the form os a new
                     Dseqrecord object.

                     Changes to how the primer_design functions work, especially
                     assembly primers.

0.6.1     2014-04-25 Fixed a bug in the Dseqrecord synced method and removed the
                     utils synced function.

0.6.0     2014-04-18 Bugfixes and improvements in documentation.

0.5.0     2013-12-16 Changes to how the amplify and assembly modules work
                     the Amplicon and Assembly classes are now subclasses of
                     Dseqrecord.

0.2.2     2013-11-05 bugfix: changed the handling of compound features
                     to fit with the new version of BioPython (1.62) which is
                     now a requirement.

0.2.1     2013-08-18 ---

0.1.8     2013-06-02 bugfix: changed the SeqFeatures added to PCR products in the
                     amplify module to a dict of list of strings instead of
                     a dict of strings.

0.1.7     2013-05-29 Changed the code in amplify.Amplicon to handle features
                     spanning the origin of circular sequences.

0.1.6     2013-04-22 Changed the behaviour of the find method of the Dseq object
                     to find substrings that span the origin. Slicing for circular
                     Dseq objects now works slightly different.

0.1.5     2013-04-18 Changed the setup.py script to permit installation
                     of the source installer without access to a c compiler.

0.1.4     2013-04-10 Cleaned up some docstrings
                     Renamed Drecord -> Dseqrecord to be more consistent with
                     Dseq and Biopython Seq/SeqRecord.

                     Changed name of keyword argument for read and parse.
                     ds=True returns Dseqrecord(s) while ds=False returns
                     SeqRecords.

0.1.3     2013-04-09 pydna created from Python-dna.
=======   ========== =====================================================================

System Requirements
===================

- `Python 2.7 <http://www.python.org>`_.
- `Biopython >= 1.65 <http://pypi.python.org/pypi/biopython>`_.
- `networkx >= 1.8.1 <http://pypi.python.org/pypi/networkx>`_.
- `appdirs >=1.3.0 <https://pypi.python.org/pypi/appdir>`_.
- `prettytable>=0.7.2 <https://pypi.python.org/pypi/PrettyTable>`_.



Python 2.x
----------

This package was developed on and for Python 2.7. Other versions have not been tested.

Python 3.x
----------

This code has not been tried with Python 3. If there
is sufficient interest, there might be a Python 3 version in the future.

Installation
============

PIP
---

The best way of installing pydna is with pip. Pip is the
officially `recommended <http://python-packaging-user-guide.readthedocs.org/en/latest/>`_ tool
for installaion of Python packages from PyPi.
Pip installs dependencies automatically.

Linux:
::

 bjorn@bjorn-UL30A:~/Dropbox/pydna$ sudo pip install pydna

Windows:
::

 C:\> pip install pydna

If you do not have pip, you can get it by following
these `instructions <http://www.pip-installer.org/en/latest/installing.html>`_.


Source
------

If you install from source, you need to install the dependencies (listed above).
Download one of the source installers from the pypi site and extract the file.
Open the pydna source code directory (containing the setup.py file) in
terminal and type:

python setup.py install

Binary distribution
-------------------

There are no binary distributions available.


Windows
-------

If dependencies have to be installed separately, this can be done using the
binary installers for Windows for those who are not comfortable at the
command line:

================ ========================================================
Dependency       Hyperlink
================ ========================================================
Python (32,64)   <http://www.python.org/download/>
Biopython (32)   <http://biopython.org/wiki/Download>
Biopython (64)   <http://www.lfd.uci.edu/~gohlke/pythonlibs/#biopython>
networkx (32,64) <http://www.lfd.uci.edu/~gohlke/pythonlibs/#networkx>
================ ========================================================


Source Code Repository
----------------------

Pydna is hosted by [Github](https://github.com/BjornFJohansson/pydna)


Distribution Structure
======================

README.txt          -- This file.

LICENSE.txt         -- What you can do with the code.

setup.py            -- Installation file.

run_tests.py        -- run tests by "python run_tests.py"<enter>

pydna/              -- The code.

docs/               -- Documentation and cookbook.

scripts/            -- Miscellaneous and perhaps useful scripts and examples.

tests/              -- Testing code.

Todo
====

* Add identification of each fragment in the Contig.small_figure method.
