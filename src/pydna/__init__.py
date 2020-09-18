#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""
:copyright: Copyright 2013 - 2019 by Björn Johansson. All rights reserved.
:license:   This code is part of the pydna package, governed by the
            license in LICENSE.txt that should be included as part of this package.

pydna
=====
Pydna is a python package providing code for simulation of the creation of recombinant DNA molecules
using `molecular biology <https://en.wikipedia.org/wiki/Molecular_biology>`_ techniques.

Provided:
  1. PCR simulation
  2. Assembly simulation based on shared identical sequences
  3. Primer design for amplification of a given sequence
  4. Automatic design of primer tails for Gibson assembly or homologous recombination.
  5. Restriction digestion and cut&paste cloning
  6. Agarose gel simulation
  7. Download sequences from Genbank
  8. Parsing various sequence formats including the capacity to handle broken Genbank format

pydna package layout
--------------------

The most important modules and how to import functions or classes from them are listed below.
Class names starts with a capital letter, functions with a lowercase letter:

::

      from pydna.module import function
      from pydna.module import Class

      Example: from pydna.gel import Gel

      pydna
         ├── amplify
         │         ├── Anneal
         │         └── pcr
         ├── assembly
         │          └── Assembly
         ├── design
         │        ├── assembly_fragments
         │        └── primer_design
         ├── download
         │          └── download_text
         ├── dseqrecord
         │            └── Dseqrecord
         ├── gel
         │     └── Gel
         ├── genbank
         │         ├── genbank
         │         └── Genbank
         ├── parsers
         │         ├── parse
         │         └── parse_primers
         └── readers
                   ├── read
                   └── read_primers



How to use the documentation
----------------------------
Documentation is available as docstrings provided in the source code for each module.
These docstrings can be inspected by reading the source code directly.
See further below on how to obtain the code for pydna.

In the python shell, use the built-in ``help`` function to view a function's docstring::

  >>> from pydna import readers
  >>> help(readers.read)
  ... # doctest: +SKIP

The doctrings are also used to provide an automaticly generated reference manual available online at
`read the docs <https://pydna.readthedocs.io>`_.

Docstrings can be explored using `IPython <http://ipython.org/>`_, an advanced Python shell with
TAB-completion and introspection capabilities. To see which functions are available in `pydna`,
type `pydna.<TAB>` (where `<TAB>` refers to the TAB key).
Use `pydna.open_config_folder?<ENTER>`to view the docstring or `pydna.open_config_folder??<ENTER>` to view the source code.

In the `Spyder IDE <https://github.com/spyder-ide/spyder>`_ it is possible to place the cursor immediately before
the name of a module,class or function and press ctrl+i to bring up docstrings in a separate window in Spyder

Code snippets are indicated by three greater-than signs::

    >>> x=41
    >>> x=x+1
    >>> x
    42

pydna source code
-----------------

The pydna source code is `available on Github <https://github.com/BjornFJohansson/pydna>`_.

How to get more help
--------------------

Please join the `Google grooup <https://groups.google.com/forum/#!forum/pydna>`_ for pydna, this is the preferred
location for help. If you find bugs in pydna itself, open an issue at the `Github repository <https://github.com/BjornFJohansson/pydna/issues>`_.

Examples of pydna in use
------------------------

See this repository for a collection of `examples <https://github.com/BjornFJohansson/pydna-examples>`_.

"""


from pydna import _version
import os as _os
import sys as _sys
import subprocess as _subprocess
import logging as _logging
import logging.handlers as _handlers
import appdirs as _appdirs
import configparser as _configparser
import prettytable as _prettytable


__author__ = "Björn Johansson"
__copyright__ = "Copyright 2013 - 2018 Björn Johansson"
__credits__ = ["Björn Johansson", "Mark Budde"]
__license__ = "BSD"
__maintainer__ = "Björn Johansson"
__email__ = "bjorn_johansson@bio.uminho.pt"
__status__ = "Development"  # "Production" #"Prototype"


# create config directory
_os.environ["pydna_config_dir"] = _os.getenv(
    "pydna_config_dir", _appdirs.user_config_dir("pydna")
)

try:
    _os.makedirs(_os.environ["pydna_config_dir"])
except OSError:
    if not _os.path.isdir(_os.environ["pydna_config_dir"]):
        raise

# set path for the pydna.ini file
_ini_path = _os.path.join(_os.environ["pydna_config_dir"], "pydna.ini")

# initiate a config parser instance
_parser = _configparser.ConfigParser()

# if a pydna.ini exists, it is read
if _os.path.exists(_ini_path):
    _parser.read(_ini_path)
else:  # otherwise it is created with default settings
    with open(_ini_path, "w", encoding="utf-8") as f:  # TODO needs encoding?
        _parser["main"] = {
            "loglevel": str(_logging.WARNING),
            "email": "someone@example.com",
            "data_dir": _appdirs.user_data_dir("pydna"),
            "log_dir": _appdirs.user_log_dir("pydna"),
            "cached_funcs": "pydna.genbank.genbank.nucleotide",
            "ape": "put/path/to/ape/here",
            "primers": "put/path/to/primers/here",
            "enzymes": "put/path/to/enzymes/here",
        }
        _parser.write(f)

# pydna related environmental variables are set from pydna.ini if they are not set already
_mainsection = _parser["main"]
_os.environ["pydna_ape"] = _os.getenv(
    "pydna_ape", _mainsection.get("ape", "put/path/to/ape/here")
)
_os.environ["pydna_cached_funcs"] = _os.getenv(
    "pydna_cached_funcs", _mainsection.get("cached_funcs", "none")
)

_os.environ["pydna_data_dir"] = _os.getenv(
    "pydna_data_dir", _mainsection.get("data_dir", _appdirs.user_data_dir("pydna"))
)
_os.environ["pydna_email"] = _os.getenv(
    "pydna_email", _mainsection.get("email", "someone@example.com")
)
_os.environ["pydna_enzymes"] = _os.getenv(
    "pydna_enzymes", _mainsection.get("enzymes", "put/path/to/enzymes/here")
)
_os.environ["pydna_log_dir"] = _os.getenv(
    "pydna_log_dir", _mainsection.get("log_dir", _appdirs.user_log_dir("pydna"))
)
_os.environ["pydna_loglevel"] = _os.getenv(
    "pydna_loglevel", _mainsection.get("loglevel", str(_logging.WARNING))
)
_os.environ["pydna_primers"] = _os.getenv(
    "pydna_primers", _mainsection.get("primers", "put/path/to/primers/here")
)


# create log directory if not present
_os.makedirs(_os.environ["pydna_log_dir"], exist_ok=True)  # changes to file system ####
_logmsg = "Log directory {}".format(_os.environ["pydna_log_dir"])

# create logger
_logger = _logging.getLogger("pydna")
_logger.setLevel(int(_os.environ["pydna_loglevel"]))
_hdlr = _handlers.RotatingFileHandler(
    _os.path.join(_os.environ["pydna_log_dir"], "pydna.log"),
    mode="a",
    maxBytes=10 * 1024 * 1024,
    backupCount=10,
    encoding="utf-8",
)
_formatter = _logging.Formatter("%(asctime)s %(levelname)s %(funcName)s %(message)s")
_hdlr.setFormatter(_formatter)
_logger.addHandler(_hdlr)
_logger.info(_logmsg)
_logger.info("Environmental variable pydna_ape          = %s", _os.environ["pydna_ape"])
_logger.info(
    "Environmental variable pydna_cached_funcs = %s", _os.environ["pydna_cached_funcs"]
)
_logger.info(
    "Environmental variable pydna_data_dir     = %s", _os.environ["pydna_data_dir"]
)
_logger.info(
    "Environmental variable pydna_email        = %s", _os.environ["pydna_email"]
)
_logger.info(
    "Environmental variable pydna_log_dir      = %s", _os.environ["pydna_log_dir"]
)
_logger.info(
    "Environmental variable pydna_loglevel     = %s", _os.environ["pydna_loglevel"]
)
_logger.info(
    "Environmental variable pydna_primers      = %s", _os.environ["pydna_primers"]
)
# create cache directory if not present

_os.makedirs(
    _os.environ["pydna_data_dir"], exist_ok=True
)  # changes to file system ####

# find out if optional dependecies for gel module are in place


def _missing_modules_for_gel():
    import importlib
    from importlib import util

    _missing = []
    for _optm in ["scipy", "numpy", "matplotlib", "mpldatacursor", "pint"]:
        _missing.extend([_optm] if not util.find_spec(_optm) else [])
    del importlib
    del util
    return _missing


_missing = _missing_modules_for_gel()

if _missing:
    _logger.warning(
        "gel simulation will NOT be available. Missing modules: %s", ", ".join(_missing)
    )
else:
    _logger.info("gel simulation is available, optional dependencies were found.")


__version__ = _version.version
del _version
_sys.modules.pop("pydna._version", None)

_logger.info("__version__ = %s", __version__)


class _PydnaWarning(Warning):
    """Pydna warning.

    Pydna uses this warning (or subclasses of it), to make it easy to
    silence all warning messages:

    >>> import warnings
    >>> from pydna import _PydnaWarning
    >>> warnings.simplefilter('ignore', _PydnaWarning)

    Consult the warnings module documentation for more details.
    """

    pass


class _PydnaDeprecationWarning(_PydnaWarning):
    """pydna deprecation warning.

    Pydna uses this warning instead of the built in DeprecationWarning
    since those are ignored by default since Python 2.7.

    To silence all our deprecation warning messages, use:

    >>> import warnings
    >>> from pydna import _PydnaDeprecationWarning
    >>> warnings.simplefilter('ignore', _PydnaDeprecationWarning)

    Code marked as deprecated will be removed in a future version
    of Pydna. This can be discussed in the Pydna google group:
    https://groups.google.com/forum/#!forum/pydna

    """

    pass


def open_current_folder():
    """Calling this function opens the current working directory
    in the default file manager. The location for this folder is
    given by the :func:`os.getcwd` function"""
    return _open_folder(_os.getcwd())


_logger.info("Current working directory = os.getcwd() = %s", _os.getcwd())


def open_cache_folder():
    """Calling this function opens the pydna cache folder in
    the default file manager. The location for this folder is stored
    in the *pydna_data_dir* environmental variable."""
    return _open_folder(_os.environ["pydna_data_dir"])


def open_config_folder():
    """Calling this function opens the pydna configuration folder in
    the default file manager. The location for this folder is stored
    in the *pydna_config_dir* environmental variable.

    The `pydna.ini` file can be edited to make pydna quicker to use.
    See the documentation of the :class:configparser.ConfigParser´ class.

    Below is the content of a typical `pydna.ini` file on a Linux
    system.

    ::

        [main]
        loglevel=30
        email=myemail@example.org
        data_dir=/home/bjorn/.local/share/pydna
        log_dir=/home/bjorn/.cache/pydna/log
        ape=tclsh /home/bjorn/.ApE/apeextractor/ApE.vfs/lib/app-AppMain/AppMain.tcl
        cached_funcs=Genbank_nucleotide
        primers=/home/bjorn/Dropbox/wikidata/PRIMERS.txt
        enzymes=/home/bjorn/Dropbox/wikidata/RestrictionEnzymes.txt

    The email address is set to someone@example.com by default. If you change this
    to you own address, the :func:`pydna.genbank.genbank` function can be used
    to download sequences from Genbank directly without having to explicitly add
    the email address.

    Pydna can cache results from the following functions or methods:

    - :func:`pydna.genbank.Genbank.nucleotide`   Genbank_nucleotide
    - :func:`pydna.amplify.Anneal`               amplify_Anneal
    - :func:`pydna.assembly.Assembly`            assembly_Assembly
    - :func:`pydna.download.download_text`       download.download_text
    - :func:`pydna.dseqrecord.Dseqrecord.synced` Dseqrecord_synced

    These can be added separated by a comma to the cached_funcs entry in **pydna.ini**
    file or the pydna_cached_funcs environment variable.

    """
    return _open_folder(_os.environ["pydna_config_dir"])


def open_log_folder():
    return _open_folder(_os.environ["pydna_log_dir"])


def _open_folder(pth):
    if _sys.platform == "win32":
        _subprocess.run(["start", pth], shell=True)
    elif _sys.platform == "darwin":
        _subprocess.run(["open", pth])
    else:
        try:
            _subprocess.run(["xdg-open", pth])
        except OSError:
            return "no cache to open."


def get_env():
    """Calling this function prints a an ascii table containing all environmental
    variables set by `pydna` and their values. Pydna related variables have names
    that starts with `pydna_`

    """
    from pydna._pretty import pretty_str as _pretty_str

    _table = _prettytable.PrettyTable(["Variable", "Value"])
    _table.set_style(_prettytable.DEFAULT)
    _table.align["Variable"] = "l"  # Left align
    _table.align["Value"] = "l"  # Left align
    _table.padding_width = 1  # One space between column edges and contents
    for k, v in sorted(_os.environ.items()):
        if k.startswith("pydna"):
            _table.add_row([k, v])
    return _pretty_str(_table)


def logo():
    """This function returns the ascii-art logotype of pydna.

    >>> import pydna
    >>> print(pydna.logo())
                     _
                    | |
     ____  _   _  __| |___   __ ___
    |  _ \| | | |/ _  |  _ \(____ |
    | |_| | |_| ( (_| | | | / ___ |
    |  __/ \__  |\____|_| |_\_____|
    |_|   (____/
    <BLANKLINE>
    """
    from pydna._pretty import pretty_str as _pretty_str
    import textwrap

    return _pretty_str(
        textwrap.dedent(
            r"""
                     _
                    | |
     ____  _   _  __| |___   __ ___
    |  _ \| | | |/ _  |  _ \(____ |
    | |_| | |_| ( (_| | | | / ___ |
    |  __/ \__  |\____|_| |_\_____|
    |_|   (____/
    """[
                1:
            ]
        )
    )


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
