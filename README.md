# ![icon](https://raw.githubusercontent.com/bjornFJohansson/pydna/master/docs/pics/pydna.resized.png) pydna

| [![Tests & Coverage](https://github.com/BjornFJohansson/pydna/actions/workflows/pydna_test_and_coverage_workflow.yml/badge.svg?branch=dev_bjorn)](https://github.com/BjornFJohansson/pydna/actions/workflows/pydna_test_and_coverage_workflow.yml) | [![codecov](https://codecov.io/gh/BjornFJohansson/pydna/branch/master/graph/badge.svg)](https://codecov.io/gh/BjornFJohansson/pydna/branch/master) | [![PyPI version](https://badge.fury.io/py/pydna.svg)](https://badge.fury.io/py/pydna)                                                  | [![Google group : pydna](https://img.shields.io/badge/Google%20Group-pydna-blue.svg)](https://groups.google.com/g/pydna)              |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| [![Documentation Status](https://github.com/BjornFJohansson/pydna/actions/workflows/publish-docs.yml/badge.svg)](https://github.com/BjornFJohansson/pydna/actions/workflows/publish-docs.yml)                                                      | [![GitHub issues](https://img.shields.io/github/issues/BjornFJohansson/pydna.svg)](https://github.com/BjornFJohansson/pydna/issues)                | [![Anaconda-Server Badge2](https://anaconda.org/bjornfjohansson/pydna/badges/license.svg)](https://anaconda.org/bjornfjohansson/pydna) | [![GitHub stars](https://img.shields.io/github/stars/BjornFJohansson/pydna.svg)](https://github.com/BjornFJohansson/pydna/stargazers) |


Planning genetic constructs with many parts and assembly steps, such as recombinant
metabolic pathways :petri_dish:, are often difficult to **properly** document as is evident from the poor
state of documentation in the scientific literature :radioactive:.


The pydna python package provide a human-readable formal descriptions of :dna: cloning and genetic assembly
strategies in Python :snake: which allow for simulation and verification.
Pydna can be used as [executable documentation](https://en.wikipedia.org/wiki/Literate_programming) for cloning.

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

Some simple examples can be found [here](https://github.com/MetabolicEngineeringGroupCBMA/pydna-examples#pydna-examples).

For more elaborate use, look at some assembly strategies of D-xylose metabolic pathways [MetabolicEngineeringGroupCBMA/ypk-xylose-pathways](https://github.com/MetabolicEngineeringGroupCBMA/ypk-xylose-pathways#pereira-et-al-2016).

See below for documentation.

## Acknowledgement 🤝

If you use pydna in your research, please reference the paper:

Pereira, F., Azevedo, F., Carvalho, Â., Ribeiro, G. F., Budde, M. W., & Johansson, B. (2015). Pydna: a simulation and documentation tool for DNA assembly strategies using python. BMC Bioinformatics, 16(142), 142. [doi:10.1186/s12859-015-0544-x](https://doi.org/10.1186/s12859-015-0544-x)

## Documentation and usage 📚

Full documentation of all modules and classes can be found at [https://bjornfjohansson.github.io/pydna](https://bjornfjohansson.github.io/pydna).

Most pydna functionality is implemented as methods for the double stranded DNA sequence record
classes Dseq and Dseqrecord, which are subclasses of the [Biopython](http://biopython.org/wiki/Main_Page)
[Seq](http://biopython.org/wiki/Seq) and [SeqRecord](http://biopython.org/wiki/SeqRecord) classes.

These classes make PCR primer design, PCR simulation and cut-and-paste cloning very simple:

[![example](https://raw.githubusercontent.com/BjornFJohansson/pydna/master/docs/example.png)](https://github.com/BjornFJohansson/pydna/blob/master/docs/example.ipynb)

As the example above shows, pydna keeps track of sticky ends and features.

Pydna can be very compact. The eleven lines of Python below simulates the construction of a recombinant plasmid.
DNA sequences are downloaded from Genbank by accession numbers that are guaranteed to be stable over time.

```python
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
```

## Installation 📦

By default, pydna is installed with minimal dependencies, but there are optional dependencies for additional functionality.

<details>
<summary>Click here to see optional dependencies</summary>

### `clipboard`

Enables the `pydna.dseqrecord.Dseqrecord.copy_gb_to_clipboard()` and `pydna.dseqrecord.Dseqrecord.copy_fasta_to_clipboard()`

These methods will put a copy the sequence on the clipboard in either Genbank (gb) or fasta format.


| Dependency                                          | Function in pydna          |
| --------------------------------------------------- | -------------------------- |
| [pyperclip](https://pypi.python.org/pypi/pyperclip) | copy sequence to clipboard |

### `download`

Pyparsing enables the `pydna.genbankfixer.gbtext_clean()` function that can automatically
correct malformed sequence files in Genbank format. These are often found online, so this
option also installs requests to enable the  `pydna.genbankfixer.download.download_text()` function which can be used to get cleaned up text from a URL.


| Dependency                                          | Function in pydna                                 |
| --------------------------------------------------- | ------------------------------------------------- |
| [pyparsing](https://pypi.python.org/pypi/pyparsing) | fix corrupt Genbank files with pydna.genbankfixer |
| [requests](https://pypi.org/project/requests)       | download sequences with pydna.download            |

### `express`

This option enables the `pydna.utils.cai()` function and the `cai()` method
available from subclasses of `pydna.seqrecord.SeqRecord`, such as
`pydna.dseqrecord.Dseqrecord`.

| [cai2](https://pypi.python.org/pypi/cai2)           | codon adaptation index calculations in several modules |

### `gel`

Scipy, matplotlib and pillow (PIL) enable the generation of gel images. Numpy is also
needed, but usually installed as a dependency of biopython.


| Dependency                                        | Function in pydna             |
| ------------------------------------------------- | ----------------------------- |
| [scipy](https://www.scipy.org)                    | gel simulation with pydna.gel |
| [matplotlib](http://matplotlib.org)               | “                             |
| [pillow](https://github.com/python-pillow/Pillow) | “                             |

</details>

### Installing with pip 🐍

```bash
# use the `--pre` flag to get the latest version of pydna.
pip install --pre --upgrade pydna

# to install the optional dependencies, you can use the following command:
pip install --pre --upgrade pydna[clipboard,download,express,gel]
```

Remove options inside the square brackets as required, but be sure not to leave spaces as pip will not recognize the options. See below under "Optional dependencies".

### Installing with poetry 🧙‍♂️

If your project uses [poetry](https://python-poetry.org/) to manage dependencies, you can install pydna with the following commands:

```bash
# Basic installation
poetry add pydna
# With optional dependencies (ommit the options you don't want)
poetry add pydna --extras "clipboard download express gel"

# If you already have it installed and you want to add or remove optional
# dependencies, you have to uninstall and install again
poetry remove pydna
poetry add pydna --extras "express gel"
```

## Contributing and feedback 🛠️

Feedback & suggestions are very welcome! Please create an issue with your question, comment or suggestion. Please include the version of pydna you are using and code to reproduce the issue if possible.

If you don't have a github account, you can get in touch through the [google group](https://groups.google.com/d/forum/pydna) for pydna.

Below are the instructions for developers who want to contribute to pydna. Please direct pull requests towards the `dev_bjorn` branch.

### Fork the repository and set up a dev branch 🍴

Fork the entire repository (not just the `master` branch by unticking the "Copy the `master` branch only" box)

<img src="docs/_static/create-fork.png" alt="Fork repository" width="400">

Create your branch starting from `dev_bjorn`, and if your changes are related to an issue, call the branch `issue_<number>`.

```bash
# Clone the repository
git clone https://github.com/<your-username>/pydna.git

# Change to the repository directory
cd pydna

# Go to the dev_bjorn branch
git checkout -b dev_bjorn

# Pull the current version of dev_bjorn
git pull origin dev_bjorn

# Create your own branch
git checkout -b issue_<number>
```

### Local development 💻

1. Use [Poetry](https://python-poetry.org/docs/#installation) to install dependencies and activate virtual environment.

    ```bash
    # If you want the virtual environment to be created in this folder
    # (this is now the default, see `poetry.toml`)
    poetry config virtualenvs.in-project true

    # Install dependencies (extras are required for tests to pass)
    poetry install --all-extras

    # Activate virtual environment
    poetry shell
    ```
2. Make your changes.
3. Add the necessary tests in `tests/`.
4. Run the tests from the root directory with `python run_test.py`.
   > **TIP:** You can run a particular test file with `pytest -vs test_file.py` (`-v` for verbose and `-s` to see print statements in the test). If you want to run just a single test, you can use `pytest -vs -k test_name`, where `test_name` is the name of the test function.
5. Install `pre-commit` hooks if you haven't by running `pre-commit install`. `pre-commit` should be available in the environment, since it is installed by poetry.
   > **TIP:** The hooks are a series of checks that will be run before you commit your code. If any of the checks fail, the commit will not be allowed. Some of them auto-fix the code (e.g., `black` formatting), so you can simply do `git add .` and commit again.
6. Push the changes to your fork
7. From your fork,make a PR towards the branch `dev_bjorn` in the original repository.

### Continuous integration 🤖

The test_and_coverage workflow is triggered on all pushed commits for all branches except the `master` branch. This workflow run tests, doctests and a series of Jupyter notebooks using pytest on Linux, Windows and macOS with all
supported python versions.

### Building the documentation locally 📚

Documentation is built using [Sphinx](http://www.sphinx-doc.org/) from [docstrings](https://www.python.org/dev/peps/pep-0257/)
using a GitHub [action](https://github.com/BjornFJohansson/pydna/actions/workflows/publish-docs.yml).
The [numpy](www.numpy.org) [docstring format](https://numpy.org/doc/stable/dev/howto-docs.html#docstring-intro) is used.

Below the commands to run a local sphinx server that auto-updated when files are changed.

```bash
# Install docs dependency group
poetry install --with docs

# Start the sphinx server to see docs live by default at http://127.0.0.1:8000/
sphinx-autobuild --watch src/ docs docs/_build/html

```

## Release process 🚀

See the [releases](https://github.com/BjornFJohansson/pydna/releases) for changes and releases.

The build workflow builds a PyPI packages using poetry. This workflow is triggered by publishing a Github release manually from the Github web interface.

![----]( http://bit.ly/coloredline)

## History 📜

Pydna was made public in 2012 on [Google code](https://code.google.com/archive/p/pydna).

:microbe:

:portugal:

## Who is using pydna? 🧪

Taylor, L. J., & Strebel, K. (2017).
Pyviko: an automated Python tool to design gene knockouts in complex viruses with overlapping genes.
BMC Microbiology, 17(1), 12.
[PubMed](https://www.ncbi.nlm.nih.gov/pubmed/28061810)

Wang, Y., Xue, H., Pourcel, C., Du, Y., & Gautheret, D. (2021).
2-kupl: mapping-free variant detection from DNA-seq data of matched samples.
In Cold Spring Harbor Laboratory (p. 2021.01.17.427048). [DOI](https://doi.org/10.1101/2021.01.17.427048)
[PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8180056)

[ShareYourCloning](https://shareyourcloning.org), a web application for designing and documenting DNA cloning strategies.

[An Automated Protein Synthesis Pipeline with Transcriptic and Snakemake](http://blog.booleanbiotech.com/transcriptic_protein_synthesis_pipeline.html)

and other projects on [github](https://github.com/BjornFJohansson/pydna/network/dependents?package_id=UGFja2FnZS01MjQ2MjYzNQ%3D%3D)
