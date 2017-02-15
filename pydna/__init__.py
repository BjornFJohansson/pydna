#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os               as _os
import sys              as _sys
import subprocess       as _subprocess
import errno            as _errno
import glob             as _glob
import logging          as _logging
import logging.handlers as _handlers
import appdirs          as _appdirs
import configparser     as _configparser
import prettytable      as _prettytable
from pydna._pretty      import pretty_str  as _pretty_str
'''
# pydna

The pydna package.

:copyright: Copyright 2013 - 2016 by Björn Johansson. All rights reserved.
:license:   This code is part of the pydna distribution and governed by its
            license.  Please see the LICENSE.txt file that should have been included
            as part of this package.

'''

__author__       = "Björn Johansson"
__copyright__    = "Copyright 2013 - 2015 Björn Johansson"
__credits__      = ["Björn Johansson", "Mark Budde"]
__license__      = "BSD"
__maintainer__   = "Björn Johansson"
__email__        = "bjorn_johansson@bio.uminho.pt"
__status__       = "Development" # "Production" #"Prototype"
from pydna._version import get_versions as _get_versions
__version__      = _get_versions()['version'][:5]
__long_version__ = _get_versions()['version']
del _get_versions


'''
Pydna caches results from the assembly dsdna and amplify
modules. pydna sets an environmental variable "pydna_cache"
which can have three different values:

"cached"        A cache directory if created.
                if possible, cached results are returned.
                if cached results are not available,
                new results are created and cached.
                This is the default.

"nocache"       Results are not written or read from cache.
                No cache directory is created.

"refresh"       new results are made and old are overwritten if they
                exist.

"compare"       The compare functionality is not implemented yet.
                Results are made new and compared with cached results.
                If cached results are not available an exception is raised.

                If new results are not identical with cached, a warning is raised
                and details written to a log located in the data_dir
                The log entry should give the name of the calling script
                if possible and as much details as possible.
                http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python/
'''



# create config directory
_os.environ["pydna_config_dir"] = _os.getenv("pydna_config_dir", _appdirs.user_config_dir("pydna"))
try:
    _os.makedirs( _os.environ["pydna_config_dir"] )
except OSError:
    if not _os.path.isdir( _os.environ["pydna_config_dir"] ):
        raise

# set path for the pydna.ini file
_ini_path = _os.path.join( _os.environ["pydna_config_dir"], "pydna.ini" )

# initiate a config parser instance
_parser = _configparser.ConfigParser()

# if a pydna.ini exists, it is read
if _os.path.exists(_ini_path):
    _parser.read(_ini_path)
else: # otherwise it is created with default settings
    with open(_ini_path, 'w') as f:
        _parser["main"] = { 'loglevel': str(_logging.WARNING),
                            'email'   : "someone@example.com",
                            'data_dir': _appdirs.user_data_dir("pydna"),
                            'log_dir' : _appdirs.user_log_dir("pydna"),
                            'cache'   : 'nocache',
                            'ape'     : 'put/path/to/ape/here',
                            'primers' : ''}
        _parser.write(f)

# Seven pydna related environmental variables are set from pydna.ini if they are not set already
_mainsection = _parser["main"]
_os.environ["pydna_loglevel"] = _os.getenv("pydna_loglevel", _mainsection.get("loglevel",str(_logging.WARNING)))
_os.environ["pydna_email"]    = _os.getenv("pydna_email",    _mainsection.get("email","someone@example.com"))
_os.environ["pydna_data_dir"] = _os.getenv("pydna_data_dir", _mainsection.get("data_dir",_appdirs.user_data_dir("pydna")))
_os.environ["pydna_log_dir"]  = _os.getenv("pydna_log_dir",  _mainsection.get("log_dir",_appdirs.user_log_dir("pydna")))
_os.environ["pydna_cache"]    = _os.getenv("pydna_cache",    _mainsection.get("cache", 'nocache'))
_os.environ["pydna_ape"]      = _os.getenv("pydna_ape",      _mainsection.get("ape",'put/path/to/ape/here'))
_os.environ["pydna_primers"]  = _os.getenv("pydna_primers",  _mainsection.get("primers", ''))



# Check sanity of pydna_cache variable
if _os.environ["pydna_cache"] not in ("cached", "nocache", "refresh", "compare"):
    raise Exception("cache (os.environ['pydna_cache']) is not cached, nocache, refresh or compare")

# create log directory if not present
try:
    _os.makedirs( _os.environ["pydna_log_dir"] )
    _logmsg = "Created log directory {}".format(_os.environ["pydna_log_dir"])
except OSError:
    if _os.path.isdir( _os.environ["pydna_log_dir"] ):
        _logmsg = "Log directory {} found.".format(_os.environ["pydna_log_dir"])
    else:
        raise

# create logger
_logger = _logging.getLogger("pydna")
_logger.setLevel( int(_os.environ["pydna_loglevel"]) )
_hdlr = _handlers.RotatingFileHandler(_os.path.join( _os.environ["pydna_log_dir"] , 'pydna.log'), mode='a', maxBytes=10*1024*1024, backupCount=10, encoding='utf-8')
_formatter = _logging.Formatter('%(asctime)s %(levelname)s %(funcName)s %(message)s')
_hdlr.setFormatter(_formatter)
_logger.addHandler(_hdlr)
_logger.info(_logmsg)
_logger.info('Assigning environmental variable pydna_data_dir = %s', _os.environ["pydna_data_dir"] )

# create cache directory if not present
try:
    _os.makedirs( _os.environ["pydna_data_dir"] )
    _logger.info("Created data directory %s", _os.environ["pydna_data_dir"] )
except OSError:
    if _os.path.isdir( _os.environ["pydna_data_dir"] ):
        _logger.info("data directory %s found", _os.environ["pydna_data_dir"])
    else:
        raise



from pydna.amplify    import Anneal
from pydna.amplify    import pcr
from pydna.amplify    import nopcr

from pydna.assembly   import Assembly

from pydna.genbank    import Genbank
from pydna.genbank    import genbank

from pydna.download   import download_text

from pydna.dseq       import Dseq
from pydna.dseqrecord import Dseqrecord

from pydna.readers    import read
from pydna.readers    import read_primer

from pydna.parsers    import parse
from pydna.parsers    import parse_primers

from pydna.editor             import Editor
from pydna.common_sub_strings import common_sub_strings

from pydna.design      import cloning_primers
from pydna.design      import assembly_primers
from pydna.design      import integration_primers

from pydna.design      import primer_design
from pydna.design      import assembly_fragments

from pydna.utils              import eq
from pydna.utils              import shift_origin
from pydna.utils              import pairwise
from pydna.utils              import cseguid

from pydna.genbankfixer       import gbtext_clean



# find out if optional dependecies for gel module are in place
_missing_modules_for_gel = []

try:
    import scipy
    del scipy
except ImportError:
    _missing_modules_for_gel.append("scipy")
try:
    import numpy
    del numpy
except ImportError:
    _missing_modules_for_gel.append("numpy")
try:
    import matplotlib
    del matplotlib
except ImportError:
    _missing_modules_for_gel.append("matplotlib")
try:
    import mpldatacursor
    del mpldatacursor
except ImportError:
    _missing_modules_for_gel.append("mpldatacursor")
try:
    import pint
    del pint
except ImportError:
    _missing_modules_for_gel.append("pint")

if _missing_modules_for_gel:
    _logger.warning("gel simulation will NOT be available. Missing modules: %s",
                     ", ".join(_missing_modules_for_gel))
else:
    _logger.info("gel simulation will be available.")
    from pydna.gel import Gel

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
    return _open_folder( _os.getcwd() )
    
def open_cache_folder():
    _open_folder( _os.environ["pydna_data_dir"] )

def open_config_folder():
    _open_folder( _os.environ["pydna_config_dir"] )

def open_log_folder():
    _open_folder( _os.environ["pydna_log_dir"] )

def _open_folder(pth):
    if _sys.platform=='win32':
        _subprocess.Popen(['start', pth], shell=True)
    elif _sys.platform=='darwin':
        _subprocess.Popen(['open', pth])
    else:
        try:
            _subprocess.Popen(['xdg-open', pth])
        except OSError:
            return "no cache to open."
            
def delete_cache(categories=[ "amplify*", "assembly*", "genbank*", "web*", "synced*" ]):
    msg = "cache file deletion\n"
    for cat in categories:
        files = _glob.glob(_os.path.join( _os.environ["pydna_data_dir"], cat))
        for file_ in files:
            msg += file_
            try:
                _os.remove( file_ )
                msg += " deleted.\n"
            except OSError as e:
                if e._errno == _errno.ENOENT:
                    msg += " no file to delete.\n"
    return _pretty_str(msg)

def set_nocache():
    _os.environ["pydna_cache"]="nocache"
def set_cached():
    _os.environ["pydna_cache"]="cached"
def set_refresh():
    _os.environ["pydna_cache"] ="refresh"
def set_compare():
    _os.environ["pydna_cache"] ="compare"

def get_env():
    _table = _prettytable.PrettyTable(["Variable", "Value"])
    _table.set_style(_prettytable.DEFAULT)
    _table.align["Variable"] = "l" # Left align
    _table.align["Value"] = "l" # Left align
    _table.padding_width = 1 # One space between column edges and contents
    for k,v in _os.environ.items():
        if k.startswith("pydna"):
            _table.add_row([k,v])
    return _pretty_str(_table)
    
logo=_pretty_str("                 _             \n"       
                 "                | |            \n"            
                 " ____  _   _  __| |___   __ ___\n"
                 "|  _ \| | | |/ _  |  _ \(____ |\n"
                 "| |_| | |_| ( (_| | | | / ___ |\n"
                 "|  __/ \__  |\____|_| |_\_____|\n"
                 "|_|   (____/                   \n")
