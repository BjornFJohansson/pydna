#!/bin/bash                 # This “shebang” tells what program to use to interpret the script.

conda config --add channels ulmo
conda config --add channels synthicity 

$PYTHON setup.py install    # Python command to install the script.
