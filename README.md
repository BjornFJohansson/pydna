# ![icon](https://raw.githubusercontent.com/bjornFJohansson/pydna/master/docs/pics/pydna.resized.png) pydna

| [![Tests & Coverage](https://github.com/BjornFJohansson/pydna/actions/workflows/pydna_test_and_coverage_workflow.yml/badge.svg?branch=dev_bjorn)](https://github.com/BjornFJohansson/pydna/actions/workflows/pydna_test_and_coverage_workflow.yml) | [![codecov](https://codecov.io/gh/BjornFJohansson/pydna/branch/master/graph/badge.svg)](https://codecov.io/gh/BjornFJohansson/pydna/branch/master) | [![PyPI version](https://badge.fury.io/py/pydna.svg)](https://badge.fury.io/py/pydna)                                                  | [![Google group : pydna](https://img.shields.io/badge/Google%20Group-pydna-blue.svg)](https://groups.google.com/g/pydna)              |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| [![Documentation Status](https://readthedocs.org/projects/pydna/badge/?version=latest)](http://pydna.readthedocs.io/?badge=latest)                                                                                                                 | [![GitHub issues](https://img.shields.io/github/issues/BjornFJohansson/pydna.svg)](https://github.com/BjornFJohansson/pydna/issues)                | [![Anaconda-Server Badge2](https://anaconda.org/bjornfjohansson/pydna/badges/license.svg)](https://anaconda.org/bjornfjohansson/pydna) | [![GitHub stars](https://img.shields.io/github/stars/BjornFJohansson/pydna.svg)](https://github.com/BjornFJohansson/pydna/stargazers) |



Planning genetic constructs with many parts and assembly steps, such as recombinant
metabolic pathways :petri_dish:, are often difficult to **properly** document as is evident from the poor
state of documentation in the scientific literature :radioactive:.


The pydna python package provide a human-readable formal descriptions of :dna: cloning and genetic assembly
strategies in Python :snake: which allow for simulation and verification. Pydna can be used as [executable documentation](https://en.wikipedia.org/wiki/Literate_programming) for cloning.


A cloning strategy expressed in pydna is **complete**, **unambiguous** and **stable**.


Pydna provides simulation of:

- Primer design
- PCR
- Restriction digestion
- Ligation
- Gel electrophoresis of DNA with generation of gel images
- Homologous recombination
- Gibson assembly
- Golden gate assembly (in progress)


Virtually any sub-cloning experiment can be described in pydna, and its execution yield
the sequences of intermediate and final DNA molecules.

Pydna has been designed with the goal of being understandable for biologists with only some basic understanding of Python.

Pydna can formalize planning and sharing of cloning strategies and is especially useful for complex or combinatorial
DNA molecule constructions.

Start by looking at the [cookbook](https://github.com/BjornFJohansson/pydna/blob/master/docs/cookbook/cookbook.ipynb).

Some simple examples can be found  [here](https://github.com/MetabolicEngineeringGroupCBMA/pydna-examples#pydna-examples).

For more elaborate use, look at some assembly strategies of D-xylose metabolic pathways [MetabolicEngineeringGroupCBMA/ypk-xylose-pathways](https://github.com/MetabolicEngineeringGroupCBMA/ypk-xylose-pathways#pereira-et-al-2016).

See below for documentation.

![----]( http://bit.ly/coloredline)



## Usage

Most pydna functionality is implemented as methods for the double stranded DNA sequence record
classes Dseq and Dseqrecord, which are subclasses of the [Biopython](http://biopython.org/wiki/Main_Page) [Seq](http://biopython.org/wiki/Seq) and [SeqRecord](http://biopython.org/wiki/SeqRecord) classes.

These classes make PCR primer design, PCR simulation and cut-and-paste cloning very simple:

[![example](https://raw.githubusercontent.com/BjornFJohansson/pydna/master/docs/example.png)](https://github.com/BjornFJohansson/pydna/blob/master/docs/example.ipynb)

As the example above shows, pydna keeps track of sticky ends and features.


Pydna can be very compact. The eleven lines of Python below simulates the construction of a recombinant plasmid.
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

Another use case for pydna in the automatic generation of [sub cloning](http://en.wikipedia.org/wiki/Subcloning) examples for teaching purposes. These examples

Feedback & suggestions are very welcome! Please post a message in the [google group](https://groups.google.com/d/forum/pydna) for pydna if you need help or have problems, questions or comments :sos:.

![----]( http://bit.ly/coloredline)


## Who is using pydna?

Taylor, L. J., & Strebel, K. (2017).
Pyviko: an automated Python tool to design gene knockouts in complex viruses with overlapping genes.
BMC Microbiology, 17(1), 12.
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28061810)


Wang, Y., Xue, H., Pourcel, C., Du, Y., & Gautheret, D. (2021).
2-kupl: mapping-free variant detection from DNA-seq data of matched samples.
In Cold Spring Harbor Laboratory (p. 2021.01.17.427048). [DOI](https://doi.org/10.1101/2021.01.17.427048)
[PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8180056)


[An Automated Protein Synthesis Pipeline with Transcriptic and Snakemake](http://blog.booleanbiotech.com/transcriptic_protein_synthesis_pipeline.html)


and other projects on [github](https://github.com/BjornFJohansson/pydna/network/dependents?package_id=UGFja2FnZS01MjQ2MjYzNQ%3D%3D)

![----]( http://bit.ly/coloredline)

There is an open access paper in BMC Bioinformatics describing pydna:

[![abstr](https://raw.githubusercontent.com/bjornFJohansson/pydna/master/docs/pics/BMC_resized.png)](http://www.biomedcentral.com/1471-2105/16/142/abstract)

Please reference the above paper:


Pereira, F., Azevedo, F., Carvalho, Â., Ribeiro, G. F., Budde, M. W., & Johansson, B. (2015). Pydna: a simulation and documentation tool for DNA assembly strategies using python. BMC Bioinformatics, 16(142), 142.


When using pydna.

![----]( http://bit.ly/coloredline)

## Documentation :page_with_curl:

Documentation is built using [Sphinx](http://www.sphinx-doc.org/) from [docstrings](https://www.python.org/dev/peps/pep-0257/)
in the code and displayed at readthedocs [![Documentation Status](https://readthedocs.org/projects/pydna/badge/?version=latest)](http://pydna.readthedocs.io/?badge=latest).
The [numpy](www.numpy.org) [docstring format](https://github.com/numpy/numpy/blob/release/doc/HOWTO_DOCUMENT.rst.txt) is used.

![----]( http://bit.ly/coloredline)

## Installation using pip

Pip is included in recent Python versions and is the
officially [recommended](http://python-packaging-user-guide.readthedocs.org/en/latest) tool.

Pip installs the minimal installation requirements automatically, but not the optional requirements (see below).

    pip install pydna

or use the --pre switch to get the latest version of pydna.

    pip install pydna --pre

for optional functionality do:

    pip install pydna[clipboard,download,express,gel]

Remove options inside the square brackets as required, but be sure not to leave spaces as pip will not recognize the options. See below under "Optional dependencies".

### Windows:

You should be able to pip install pydna from the Windows terminal as biopython now can be installed with pip as well.

    C:\> pip install pydna

By default python and pip are not on the PATH. You can re-install Python and select this option during installation, or give the full path for pip. Try something like this, depending on where your copy of Python is installed:

    C:\Python37\Scripts\pip install pydna

![----]( http://bit.ly/coloredline)

## Source Code

Pydna is developed on [Github](https://github.com/BjornFJohansson/pydna) :octocat:.
I am happy to collaborate on new features or bugfixes.

![----]( http://bit.ly/coloredline)

## Minimal installation dependencies

The list below is the minimal requirements for installing pydna.
Biopython and pydivsufsort has c-extensions, but the other modules are pure python.

- [Python 3.8, 3.9, 3.10, 3.11 or 3.12](http://www.python.org)
- [appdirs](https://pypi.python.org/pypi/appdirs)
- [biopython](http://pypi.python.org/pypi/biopython)
- [networkx](http://pypi.python.org/pypi/networkx)
- [prettytable](https://pypi.python.org/pypi/PrettyTable)
- [pydivsufsort](https://pypi.python.org/pypi/pydivsufsort)
- [pyfiglet](https://pypi.python.org/pypi/pyfiglet)

Pydna is importable even without pyfiglet.

## Optional dependencies

These can be installed `pip install pydna[clipboard,gel,download,express]`
where `[clipboard,gel,download,express]` is the list of options available. Any
combination of the words inside the square brackets are allowed, but no white space.


### `clipboard`

Enables the `pydna.dseqrecord.Dseqrecord.copy_gb_to_clipboard()` and `pydna.dseqrecord.Dseqrecord.copy_fasta_to_clipboard()`

These methods will put a copy the sequence on the clipboard in either Genbank (gb) or fasta format.


| Dependency                                          | Function in pydna                                      |
| --------------------------------------------------- | ------------------------------------------------------ |
| [pyperclip](https://pypi.python.org/pypi/pyperclip) | copy sequence to clipboard                             |

### download

Pyparsing enables the `pydna.genbankfixer.gbtext_clean()` function that can automatically
correct malformed sequence files in Genbank format. These are often found online, so this option also installs requests to enable the  `pydna.genbankfixer.download.download_text()` function which can be used to get cleaned up text from a URL.


| Dependency                                          | Function in pydna                                      |
| --------------------------------------------------- | ------------------------------------------------------ |
| [pyparsing](https://pypi.python.org/pypi/pyparsing) | fix corrupt Genbank files with pydna.genbankfixer      |
| [requests](https://pypi.org/project/requests)       | download sequences with pydna.download                 |

### express

This option enables the `pydna.utils.cai()` function and the `cai()` method
available from subclasses of `pydna.seqrecord.SeqRecord`, such as
`pydna.dseqrecord.Dseqrecord`.

| [cai2](https://pypi.python.org/pypi/cai2)           | codon adaptation index calculations in several modules |

### gel

Scipy, matplotlib and pillow (PIL) enable the generation of gel images. Numpy is also
needed, but usually installed as a dependency of biopython.


| Dependency                                          | Function in pydna                                      |
| --------------------------------------------------- | ------------------------------------------------------ |
| [scipy](https://www.scipy.org)                      | gel simulation with pydna.gel                          |
| [matplotlib](http://matplotlib.org)                 | “                                                      |
| [pillow](https://github.com/python-pillow/Pillow)   | “                                                      |


## Requirements for running tests, coverage and profiling

- [pytest](https://pypi.org/project/pytest)
- [pytest-cov](https://pypi.org/project/pytest-cov)
- [pytest-doctestplus](https://pypi.org/project/pytest-doctestplus)
- [pytest-profiling](https://pypi.org/project/pytest-profiling)
- [coverage](https://pypi.org/project/coverage)
- [nbval](https://pypi.org/project/nbval)
- [requests-mock](https://pypi.org/project/requests-mock)

for instance by `pip install pytest pytest-cov pytest-doctestplus pytest-profiling coverage nbval requests-mock`

Running the entire test suite also require:

- scipy
- matplotlib
- pillow
- pyparsing
- requests
- cai2

That can be installed by `pip install pydna[clipboard,gel,download,express]`

or by `pip install pyparsing requests cai2 scipy matplotlib pillow`


![----]( http://bit.ly/coloredline)

## Contributing

Please direct pull requests towards the `develop` branch.


### Local development

1. Use [Poetry](https://pypi.org/project/poetry) to install dependencies and activate virtual environment.

    ```bash
    # If you want the virtual environment to be created in this folder
    poetry config virtualenvs.in-project true

    # Install dependencies (extras are required for tests to pass)
    poetry install --all-extras

    # Activate virtual environment
    poetry shell
    ```

2. Make your changes.
3. Add the necessary tests in `tests/`.
4. Run the tests from the root directory with `python run_test.py`.

#### Building the documentation locally

Below the commands to run a local sphinx server that auto-updated when files are changed.

```
# Install docs dependency group
poetry install --with docs

# Start the sphinx server to see docs live by default at http://127.0.0.1:8000/
sphinx-autobuild --watch src/ docs docs/_build/html

```

## Automatic testing & Release process

See the [releases](https://github.com/BjornFJohansson/pydna/releases) for changes and releases.

There are two github actions for this package:

- `.github/workflows/pydna_test_and_coverage_workflow.yml`
- `.github/workflows/pydna_pypi_build_workflow.yml`

The test_and_coverage workflow is triggered on all pushed commits for all branches except the `master` branch. This workflow run tests, doctests and a series of Jupyter notebooks using pytest on Linux, Windows and macOS with all
supported python versions.

The build workflow builds a PyPI packages using poetry. This workflow is triggered by publishing a Github release manually from the Github web interface.

![----]( http://bit.ly/coloredline)

## Building a PyPI package with [Poetry](https://pypi.org/project/poetry)

1. Commit changes to git
2. Tag the commit according to the [Semantic Versioning](https://semver.org) format, for example "v2.0.1a3". Do not forget the "v" or poetry will not recognize the tag.

        git tag v2.0.1a3

3. Pydna uses the poetry [poetry-dynamic-versioning](https://pypi.org/project/poetry-dynamic-versioning) plugin.

        poetry dynamic-versioning # This sets the version number in the source files

4. Verify the version

        poetry version

5. Build package:

        poetry build # run this command in the root directory where the pyproject.toml file is located

6. Verify the filename of the files in the dist/ folder, they should match

7. Publish to pypi

        poetry publish

![----]( http://bit.ly/coloredline)

## History

Pydna was made public in 2012 on [Google code](https://code.google.com/archive/p/pydna).

:microbe:

:portugal:
