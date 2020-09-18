# ![icon](https://raw.githubusercontent.com/bjornFJohansson/pydna/master/docs/pics/pydna.resized.png) pydna

| ![Tests & Coverage](https://github.com/BjornFJohansson/pydna/workflows/Tests%20&%20Coverage/badge.svg)                                  |[![codecov](https://codecov.io/gh/BjornFJohansson/pydna/branch/master/graph/badge.svg)](https://codecov.io/gh/BjornFJohansson/pydna)    | [![PyPI version](https://badge.fury.io/py/pydna.svg)](https://badge.fury.io/py/pydna)                                                   |[![Anaconda-Server Badge0](https://anaconda.org/bjornfjohansson/pydna/badges/version.svg)](https://anaconda.org/bjornfjohansson/pydna)  | [![Google group : pydna](https://img.shields.io/badge/Google%20Group-pydna-blue.svg)](https://groups.google.com/g/pydna)        |
|-----------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| [![Documentation Status](https://readthedocs.org/projects/pydna/badge/?version=latest)](http://pydna.readthedocs.io/?badge=latest)      |[![GitHub issues](https://img.shields.io/github/issues/BjornFJohansson/pydna.svg)](https://github.com/BjornFJohansson/pydna/issues)     | [![Anaconda-Server Badge2](https://anaconda.org/bjornfjohansson/pydna/badges/license.svg)](https://anaconda.org/bjornfjohansson/pydna)  |[![GitHub stars](https://img.shields.io/github/stars/BjornFJohansson/pydna.svg)](https://github.com/BjornFJohansson/pydna/stargazers)   |                                                                                                                                 |


Planning genetic constructs with many parts and assembly steps, such as recombinant
metabolic pathways :petri_dish:, are often difficult to **properly** document as is evident from the
state of such documentation in the scientific literature :radioactive:.


The pydna python package provide a human-readable formal descriptions of :dna: cloning and genetic assembly
strategies in Python :snake: which allow for simulation and verification.


Pydna can be though of as [executable documentation](https://en.wikipedia.org/wiki/Literate_programming) for cloning.


A cloning strategy expressed in pydna is **complete**, **unambiguous** and **stable**.


Pydna provides simulation of:

- Restriction digestion
- Ligation
- PCR
- Primer design
- Gibson assembly
- Golden gate assembly
- Homologous recombination
- Gel electrophoresis of DNA with generation of gel images

Virtually any sub-cloning experiment can be described in pydna, and its execution yield
the sequence of the of intermediate and final resulting DNA molecule(s).

Pydna has been designed to be understandable for biologists with only some basic understanding of Python.

Pydna can formalize planning and sharing of cloning strategies and is especially useful for complex or combinatorial
DNA molecule constructions.

Look at some assembly strategies made in the Jupyter notebook
format [here](http://nbviewer.ipython.org/github/BjornFJohansson/ypk-xylose-pathways/blob/release/index.ipynb).

There is an open access paper in BMC Bioinformatics describing pydna:

[![abstr](https://raw.githubusercontent.com/bjornFJohansson/pydna/release/docs/pics/BMC_resized.png)](http://www.biomedcentral.com/1471-2105/16/142/abstract)

Please reference the above paper when using pydna.

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)


## Usage

Most pydna functionality is implemented as methods for the double stranded DNA sequence record
classes Dseq and Dseqrecord, which are subclasses of the [Biopython](http://biopython.org/wiki/Main_Page) [Seq](http://biopython.org/wiki/Seq) and [SeqRecord](http://biopython.org/wiki/SeqRecord) classes.

These classes make cut and paste cloning and PCR very simple:

    >>> from pydna.dseq import Dseq
    >>> seq = Dseq("GGATCCAAA","TTTGGATCC",ovhg=0)
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

Notably, homologous recombination and Gibson assembly between linear DNA fragments
can be easily simulated without any additional information besides the primary sequence of the fragments.

Gel electrophoresis of DNA fragments can be simulated using the included gel module
by [Bruno Silva](https://github.com/bruno2git):

![alt text](https://raw.githubusercontent.com/bjornFJohansson/pydna/release/docs/pics/gel.png "simulated agarose gel")

Look at an example notebook with a gel simulation [here](http://nbviewer.jupyter.org/github/BjornFJohansson/pydna/blob/release/scripts/gel_inline_ex.ipynb).

Pydna can be very compact. The nine lines of Python below, simulates the construction of a recombinant plasmid.
DNA sequences are downloaded from Genbank by accession numbers that are guaranteed to be stable over time.

    from pydna.genbank import Genbank
    gb = Genbank("myself@email.com") # Tell Genbank who you are!
    gene = gb.nucleotide("X06997") # Kluyveromyces lactis LAC12 gene for lactose permease.
    from pydna.parsers import parse_primers
    primer_f,primer_r = parse_primers(''' >760_KlLAC12_rv (20-mer)
                                          ttaaacagattctgcctctg

                                          >759_KlLAC12_fw (19-mer)
                                          aaatggcagatcattcgag ''')
    from pydna.amplify import pcr
    pcr_prod = pcr(primer_f,primer_r, gene)
    vector = gb.nucleotide("AJ001614") # pCAPs cloning vector
    from Bio.Restriction import EcoRV
    lin_vector = vector.linearize(EcoRV)
    rec_vec =  ( lin_vector + pcr_prod ).looped()

Pydna can automate the simulation of [sub cloning](http://en.wikipedia.org/wiki/Subcloning) experiments using
python. This is helpful to generate examples for teaching purposes.

Read the documentation (below) or the [cookbook](https://www.dropbox.com/sh/4re9a0wk03m95z4/AABpu4zwq4IuKUvK0Iy9Io0Fa?dl=0) with example files
for further information.

Please post a message in the [google group](https://groups.google.com/d/forum/pydna)
for pydna if you need help or have problems, questions or comments :sos:.

Feedback & suggestions are very welcome!
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Who is using pydna?

[An Automated Protein Synthesis Pipeline with Transcriptic and Snakemake](http://blog.booleanbiotech.com/transcriptic_protein_synthesis_pipeline.html)

[Pyviko: an automated Python tool to design gene knockouts in complex viruses with overlapping genes](https://www.ncbi.nlm.nih.gov/pubmed/28061810)
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Documentation

Documentation is built using [Sphinx](http://www.sphinx-doc.org/) from [docstrings](https://www.python.org/dev/peps/pep-0257/)
in the code and displayed at readthedocs [![Documentation Status](https://readthedocs.org/projects/pydna/badge/?version=latest)](http://pydna.readthedocs.io/?badge=latest)

The [numpy](www.numpy.org) [docstring format](https://github.com/numpy/numpy/blob/release/doc/HOWTO_DOCUMENT.rst.txt) is used.
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Installation using conda on Anaconda

The absolutely best way of installing and using pydna is to use the
free [Anaconda](https://store.continuum.io/cshop/anaconda) or [Miniconda](http://conda.pydata.org/miniconda.html) python distributions.

Anaconda is a large download (about 400 Mb) while Miniconda is about 40-50 Mb.

Once Anaconda (or Miniconda) is installed, the conda package manager can be used to install pydna.
Pydna and its dependencies are available from the [BjornFJohansson](https://anaconda.org/bjornfjohansson) package channel at
[Anaconda.org](https://anaconda.org).

The first step is to add the channel by typing the command below followed by return:

    conda config --append channels BjornFJohansson

Then pydna can be installed by typing the command below followed by return:

    conda install pydna

This works on Windows, MacOSX and Linux, and installs all necessary and optional dependencies automatically (see below).
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Installation using pip

The second best way of installing pydna is with pip, the
officially [recommended](http://python-packaging-user-guide.readthedocs.org/en/latest) tool.

Pip is included in recent Python versions.

Pip installs the minimal installation requirements automatically, but not the optional requirements (see below).
This will probably not work directly on windows, as biopython is not directly installable.

### Linux:

    bjorn@bjorn-UL30A:~/pydna$ sudo pip install pydna

### Windows:

You should be able to pip install pydna from the Windows terminal as biopython now can be installed with pip as well.

    C:\> pip install pydna

By default python and pip are not on the PATH. You can re-install Python and select this option, or give
the full path for pip. Try something like this, depending on where your copy of Python is installed:

    C:\Python37\Scripts\pip install pydna
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Installation from Source

If you install from source, you need to install all dependencies separately (listed above).
Download one of the source installers from the pypi site or from Github and extract the file.
Open the pydna source code directory (containing the setup.py file) in
terminal and type:

    python setup.py install
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Source Code

Pydna is developed on [Github](https://github.com/BjornFJohansson/pydna) :octocat:.
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Minimal installation requirements

Pydna is currently developed on and for Python 3.6 or 3.7. Pydna versions before 1.0.0 were compatible with python 2.7 only.
The list below is the minimal requirements for installing pydna. Biopython has c-extensions, but the other modules are pure python.

- [Python 3.6, 3.7 or 3.8](http://www.python.org)
- [biopython >= 1.65](http://pypi.python.org/pypi/biopython)
- [networkx >= 1.8.1](http://pypi.python.org/pypi/networkx)
- [pyparsing >= 2.1.10](https://pypi.python.org/pypi/pyparsing)
- [appdirs >=1.3.0](https://pypi.python.org/pypi/appdirs)
- [prettytable>=0.7.2](https://pypi.python.org/pypi/PrettyTable)
- [requests](https://pypi.org/project/requests)
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Optional Requirements

Pydna has been designed to be used from the Jupyter notebook. If [IPython](https://ipython.org/)
and [Jupyter](http://jupyter.org/) are installed, importing ipython notebooks as modules among are
supported among other things.

If the modules listed below are installed, gel simulation functionality will be available.

- [numpy](http://www.numpy.org)
- [scipy](https://www.scipy.org)
- [matplotlib](http://matplotlib.org)
- [mpldatacursor](https://pypi.python.org/pypi/mpldatacursor)
- [pint >= 0.7.2](https://pypi.python.org/pypi/pint)

The pydna conda package installs the optional requirements listed above as well as:

- [ipython](https://pypi.python.org/pypi/ipython)
- [jupyter](https://pypi.python.org/pypi/jupyter)
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Requirements for running tests

- [pytest>=6.0.2](https://pypi.python.org/pypi/pytest)
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Requirements for analyzing code coverage

- [python-coveralls >= 2.9.0](https://pypi.python.org/pypi/python-coveralls)
- [coverage >= 3.7.1](https://pypi.python.org/pypi/coverage)
- [pytest-cov >= 2.3.1](https://pypi.python.org/pypi/pytest-cov)
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Automatic testing

The test suit is run automatically after each commit on Linux, macOS and Windows using a GitHub action.
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Changelog

See the [change log](docs/CHANGELOG.md) for recent changes.
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)]
## Release process

There are three github actions associated with this package:

- pydna_test_and_coverage_workflow.yml
- pydna_setuptools_build_workflow.yml
- pydna_conda_build_workflow.yml

The pydna_test_and_coverage_workflow.yml is triggered on all pushed non-tagged commits.

The two other commits are triggered by publishing a github release.




:microbe:  :portugal:
