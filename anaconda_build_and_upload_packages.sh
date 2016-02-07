#!/bin/bash

echo "Building anaconda packages"

conda build .

#conda convert ~/anaconda/conda-bld/linux-64/pydna*.tar.bz2 -p all --output-dir /home/ubuntu/src/github.com/BjornFJohansson/miniconda/conda-bld/

conda convert ~/anaconda/conda-bld/linux-64/pydna*.tar.bz2 -p all --output-dir ~/anaconda/conda-bld/

anaconda upload ~/anaconda/conda-bld/linux-64/pydna* ~/anaconda/conda-bld/
anaconda upload ~/anaconda/conda-bld/linux-32/pydna* ~/anaconda/conda-bld/
anaconda upload ~/anaconda/conda-bld/win-32/pydna* ~/anaconda/conda-bld/
anaconda upload ~/anaconda/conda-bld/win-64/pydna* ~/anaconda/conda-bld/
anaconda upload ~/anaconda/conda-bld/osx-64/pydna* ~/anaconda/conda-bld/

$SHELL



# http://docs.binstar.org/conda.html
# If you have previously generated TOKEN (check Token Generation) then you may run upload process also in this way:
# $ binstar -t ${TOKEN} upload /home/USERNAME/anaconda/conda-bld/linux64/conda_gc_test-1.2.1-py27_3.tar.bz2
# http://docs.binstar.org/conda.html
