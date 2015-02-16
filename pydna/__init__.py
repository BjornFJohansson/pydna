#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2013,2014 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''
    Python-dna
    ~~~~~~~~~~

    The Python-dna package.

    :copyright: Copyright 2013,2014 by Björn Johansson. All rights reserved.
    :license:   This code is part of the Python-dna distribution and governed by its
                license.  Please see the LICENSE.txt file that should have been included
                as part of this package.

'''

__version__      = "0.8.0"
__date__         = "2015-02-16"
__author__       = "Björn Johansson"
__copyright__    = "Copyright 2013, 2014, 2015 Björn Johansson"
__credits__      = ["Björn Johansson", "Mark Budde"]
__license__      = "BSD"
__maintainer__   = "Björn Johansson"
__email__        = "bjorn_johansson@bio.uminho.pt"
__status__       = "Development" # "Production" #"Prototype"



'''
Pydna caches results from the assembly2 dsdna and amplify
modules. pydna sets an environmental variable "pydna_cache"
which can have three different values:

"cached"        if possible, cached results are returned
                if chashed results are not available
                new results are created and cached.
                This is the default.

"nocache"       results are not cached or read from cache

"refresh"       new results are made and old are overwritten if they
                exist.

"compare"       The compare functionality is not implemented yet.
                Results are made new and compared with cached results.
                If cached results are not available an exception is raised.

                If new results are not identical with cached, a warning is raised
                and details written to a log located in the data_dir
                The log entry should give the name of the calling script if possible and as
                much details as possible.
                http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python/
'''

import os
import appdirs
import logging
import logging.handlers

os.environ["pydna_cache"]  = "cached"

# set data directory depending on environment
if os.getenv("DRONE") or os.getenv("CI"):
    datadir = os.path.join(os.getcwd(),"DATA")
else:
    datadir = appdirs.user_data_dir("pydna")

# create data directory
try:
    os.mkdir(datadir)
except OSError:
    if not os.path.isdir( datadir ):
        raise
user_name  = "username not set"
user_email = "user email not set"

# create logger
logger = logging.getLogger("pydna")

hdlr = logging.handlers.RotatingFileHandler(os.path.join(datadir, 'pydna.log'), mode='a', maxBytes=0, backupCount=0, encoding='utf-8', delay=0)

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

hdlr.setFormatter(formatter)

logger.addHandler(hdlr)

logger.setLevel(logging.WARNING)

logger.info('Assigning environmental variable datadir = {}'.format(datadir))

os.environ["datadir"] = datadir

if not os.getenv("pydna_dna_dir"):
    os.environ["pydna_dna_dir"] = u''

if not os.getenv("pydna_dna_dirs"):
    os.environ["pydna_dna_dirs"] = u''


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


