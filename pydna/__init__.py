#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013,2014, 2015 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''
    pydna
    ~~~~~

    The pydna package.

    :copyright: Copyright 2013 - 2015 by Björn Johansson. All rights reserved.
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
from ._version import get_versions
__version__      = get_versions()['version'][:5]
__long_version__ = get_versions()['version']
del get_versions


'''
Pydna caches results from the assembly2 dsdna and amplify
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

import os            as _os
import sys           as _sys
import errno         as _errno
import subprocess    as _subprocess
import appdirs       as _appdirs
from configparser import SafeConfigParser as _SafeConfigParser

_os.environ["pydna_config_dir"] = _os.getenv("pydna_config_dir") or _appdirs.user_config_dir("pydna")

# create config directory
try:
    _os.makedirs( _os.environ["pydna_config_dir"] )
except OSError:
    if not _os.path.isdir( _os.environ["pydna_config_dir"] ):
        raise

_ini_path = _os.path.join( _os.environ["pydna_config_dir"], "pydna.ini" )

_parser = _SafeConfigParser()

if _os.path.exists(_ini_path):
    _parser.read(_ini_path)
else:
    with open(_ini_path, 'w') as f:
        _parser.add_section('main')
        _parser.set('main','email', "someone@example.com")
        _parser.set('main','data_dir', _appdirs.user_data_dir("pydna"))
        _parser.set('main','log_dir',  _appdirs.user_log_dir("pydna"))
        _parser.set('main','cache','cached')
        _parser.set('main','ape','put/path/to/ape/here')
        _parser.set('main','primers','')
        _parser.write(f)


_os.environ["pydna_email"]    = _os.getenv("pydna_email")    or _parser.get("main", "email")
_os.environ["pydna_data_dir"] = _os.getenv("pydna_data_dir") or _parser.get("main", "data_dir")
_os.environ["pydna_log_dir"]  = _os.getenv("pydna_log_dir")  or _parser.get("main", "log_dir")
_os.environ["pydna_cache"]    = _os.getenv("pydna_cache")    or _parser.get("main", "cache")

_os.environ["pydna_ape"]      = _os.getenv("pydna_ape")      or _parser.get("main", "ape")
_os.environ["pydna_primers"]  = _os.getenv("pydna_primers")  or _parser.get("main", "primers")

#_os.environ["pydna_cache"] = _os.getenv("pydna_cache") or "cached"

if _os.environ["pydna_cache"] not in ("cached", "nocache", "refresh", "compare"):
    raise Exception("cache (os.environ['pydna_cache']) is not cached, nocache, refresh or compare")

#_os.environ["pydna_data_dir"] = _os.getenv("pydna_data_dir") or _appdirs.user_data_dir("pydna")

# create data directory
try:
    _os.makedirs( _os.environ["pydna_data_dir"] )
except OSError:
    if not _os.path.isdir( _os.environ["pydna_data_dir"] ):
        raise

# create log directory
try:
    _os.makedirs( _os.environ["pydna_log_dir"] )
except OSError:
    if not _os.path.isdir( _os.environ["pydna_log_dir"] ):
        raise

# create logger
import logging as _logging
import logging.handlers as _handlers
_logger = _logging.getLogger("pydna")
#_logger.setLevel(_logging.DEBUG)
_hdlr = _handlers.RotatingFileHandler(_os.path.join( _os.environ["pydna_log_dir"] , 'pydna.log'), mode='a', maxBytes=10*1024*1024, backupCount=10, encoding='utf-8')
_formatter = _logging.Formatter('%(asctime)s %(levelname)s %(message)s')
_hdlr.setFormatter(_formatter)
_logger.addHandler(_hdlr)

_logger.info('Assigning environmental variable pydna_data_dir = {}'.format( _os.environ["pydna_data_dir"] ))


from pydna.amplify                                  import Anneal
from pydna.amplify                                  import pcr
from pydna.amplify                                  import nopcr
from pydna.assembly                                 import Assembly
from pydna.download                                 import Genbank
from pydna.download                                 import genbank
from pydna.download                                 import Web
from pydna.download                                 import parse_url
from pydna.download                                 import read_url
from pydna.dsdna                                    import Dseq
from pydna.dsdna                                    import Dseqrecord
from pydna.dsdna                                    import parse
from pydna.dsdna                                    import parse2
from pydna.dsdna                                    import read

#from pydna.editor                                   import Editor
from pydna.findsubstrings_suffix_arrays_python      import common_sub_strings
from pydna.primer_design                            import cloning_primers
from pydna.primer_design                            import assembly_primers
from pydna.primer_design                            import print_primer_pair
from pydna.primer_design                            import integration_primers
from pydna.utils                                    import copy_features
from pydna.utils                                    import eq
from pydna.utils                                    import shift_origin
from pydna.utils                                    import pairwise
from pydna.utils                                    import cseguid
from pydna.primer_design                            import Primer
from pydna._pretty                                  import pretty_str as _pretty_str


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

if _missing_modules_for_gel:
    _logger.warning("gel simulation will not be available. Missing modules: {}"
        .format(", ".join(_missing_modules_for_gel)))
else:
    from pydna.gel import Gel




#numpy>=1.10.1
#matplotlib>=1.5.0
#scipy>=0.16.0
#mpldatacursor>=0.6.1

#from pydna.gel                                      import gen_sample
#from pydna.gel                                      import weight_standards
#from pydna.gel                                      import weight_standard_sample
#from pydna.gel                                      import lin_div_Qty
#from pydna.gel                                      import random_Dseqs
#from pydna.gel                                      import ureg
#from pydna.gel                                      import Q_

try:
    del primer_design
except NameError:
    pass

try:
    del assembly
except NameError:
    pass

try:
    del dsdna
except NameError:
    pass

try:
    del editor
except NameError:
    pass

try:
    del amplify
except NameError:
    pass

try:
    del utils
except NameError:
    pass

try:
    del download
except NameError:
    pass

try:
    del _simple_paths8
except NameError:
    pass

try:
    del py_rstr_max
except NameError:
    pass

try:
    del findsubstrings_suffix_arrays_python
except NameError:
    pass

def delete_cache(delete=[ "amplify", "assembly", "genbank", "web", "synced" ]):
    msg = ""
    for file_ in delete:
        msg += file_
        try:
            _os.remove( _os.path.join( _os.environ["pydna_data_dir"], file_) )
            msg += " deleted.\n"
        except OSError as e:
            if e._errno == _errno.ENOENT:
                msg += " no file to delete.\n"
    return _pretty_str(msg)

def open_cache_folder():
    _open_folder( _os.environ["pydna_data_dir"] )

def open_config_folder():
    _open_folder( _os.environ["pydna_config_dir"] )

def open_log_folder():
    _open_folder( _os.environ["pydna_log_dir"] )



def _open_folder(pth):
    if _sys.platform=='win32':
        _subprocess.Popen(['start', pth], shell= True)
    elif _sys.platform=='darwin':
        _subprocess.Popen(['open', pth])
    else:
        try:
            _subprocess.Popen(['xdg-open', pth])
        except OSError:
            return "no cache to open."

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
