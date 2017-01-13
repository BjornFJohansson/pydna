#!/bin/bash                 # This “shebang” tells what program to use to interpret the script.
echo "build.sh running"
$PYTHON setup.py build
$PYTHON setup.py install    # Python command to install the script.
