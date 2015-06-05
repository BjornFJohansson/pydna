#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
import os
import sys

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
dna_dir  = "/home/bjorn/Dropbox/Public/pydna-DNA-assembly/constructs"
dna_dirs = [x[0] for x in os.walk(dna_dir)]

sys.path = dna_dirs + sys.path
os.environ["pydna_dna_dirs"] = os.pathsep.join(dna_dirs)
os.environ["PYTHONPATH"] = os.environ["pydna_dna_dirs"]

import pydna
'''

# I've implemented a workaround that modifies sys.path in my interpreter's startup script.
# This was an extension of my solution to another issue regarding shell environment
# variable inheritance [1]. It will take the user's PYTHONPATH as specified in either
# /etc/profile or ~/.profile and and insert it into the sys.path. If you don't need
# system environment variables transferred to Spyder's os.environ, just delete/comment
# out that portion of the script. My script is as follows:

'''
import os
import sys
import subprocess

envstr = subprocess.check_output('source /etc/profile; source ~/.profile; printenv', shell=True)
env = dict([a.split('=') for a in envstr.strip().split('\n')])
os.environ.update(env)

if env.has_key('PYTHONPATH'):
    pythonpath = env['PYTHONPATH'].split(':')
    pythonpath.reverse()
    for path in pythonpath:
        sys.path.insert(2, path)
'''

#os.environ["PYTHONPATH"] = os.pathsep.join(x[0] for x in os.walk(dna_dir) if not x[0].endswith("archive"))
#pYPKa_A_KlLAC4
#sys.path.append(os.path.expanduser('~/code/eol_hsrl_python'))
#os.environ['HSRL_INSTRUMENT']='gvhsrl'
#os.environ['HSRL_CONFIG']=os.path.expanduser('~/hsrl_config')
