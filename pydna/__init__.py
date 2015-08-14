#!/usr/bin/env python
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

__author__       = u"Björn Johansson"
__copyright__    = u"Copyright 2013 - 2015 Björn Johansson"
__credits__      = [u"Björn Johansson", u"Mark Budde"]
__license__      = u"BSD"
__maintainer__   = u"Björn Johansson"
__email__        = u"bjorn_johansson@bio.uminho.pt"
__status__       = u"Development" # "Production" #"Prototype"
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

import os
import sys
import shutil
import errno
import subprocess
import appdirs

os.environ["pydna_cache"] = os.getenv("pydna_cache") or "cached"

if os.environ["pydna_cache"] not in ("cached", "nocache", "refresh", "compare"):
    raise Exception("cache (os.environ['pydna_cache']) is not either cached, nocache, refresh or compare")

os.environ["pydna_data_dir"] = os.getenv("pydna_data_dir") or appdirs.user_data_dir("pydna").encode(sys.getfilesystemencoding())

# create data directory

try:
    os.makedirs( os.environ["pydna_data_dir"] )
except OSError:
    if not os.path.isdir( os.environ["pydna_data_dir"] ):
        raise


# create logger
import logging
import logging.handlers
logger = logging.getLogger("pydna")
hdlr = logging.handlers.RotatingFileHandler(os.path.join( os.environ["pydna_data_dir"] , 'pydna.log'), mode='a', maxBytes=0, backupCount=0, encoding='utf-8', delay=0)
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)
logger.info('Assigning environmental variable pydna_data_dir = {}'.format( os.environ["pydna_data_dir"] ))


from pydna.amplify                                  import Anneal
from pydna.amplify                                  import pcr
from pydna.amplify                                  import nopcr
from pydna.assembly                                 import Assembly
from pydna.download                                 import Genbank
from pydna.download                                 import web
from pydna.dsdna                                    import Dseq
from pydna.dsdna                                    import Dseqrecord
from pydna.dsdna                                    import parse
from pydna.dsdna                                    import parse2
from pydna.dsdna                                    import read
from pydna.editor                                   import Editor
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
from pydna.pretty                                   import pretty_str

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
            os.remove( os.path.join( os.environ["pydna_data_dir"], file_) )
            msg += " deleted.\n"
        except OSError, e:
            if e.errno == errno.ENOENT:
                msg += " no file to delete.\n"
    return pretty_str(msg)

def open_cache():
    if sys.platform=='win32':
        subprocess.Popen(['start', os.environ["pydna_data_dir"]], shell= True)
    elif sys.platform=='darwin':
        subprocess.Popen(['open', os.environ["pydna_data_dir"]])
    else:
        try:
            subprocess.Popen(['xdg-open', os.environ["pydna_data_dir"]])
        except OSError:
            return "no cache to open."


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
