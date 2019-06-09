#!/usr/bin/env bash

#export TWINE_REPOSITORY="https://test.pypi.org"
#condalabel="test"
condalabel="main"







echo "=============================================================="
echo "BASH_VERSION" $BASH_VERSION
echo $(git --version)
tagname="$(git describe --abbrev=0 --tags)"
tag="$(git rev-list $tagname | head -n 1)"
com="$(git rev-parse HEAD)"
dirty=$(git describe --tags --dirty --always)
msg=$(git log -1 --pretty=%B)
echo "=============================================================="
echo "Establish git variables:"
echo "=============================================================="
echo "Current commit hash : $com"
echo "Dirty tag           : $dirty"
echo "Commit msg          : $msg"
echo "=============================================================="
echo "Environment variables:"
echo "=============================================================="
echo "CI                   = $CI"
echo "APPVEYOR             = $APPVEYOR"
echo "CIRCLECI             = $CIRCLECI"
echo "TRAVIS               = $TRAVIS"
echo "CODESHIP             = $CI_NAME"
echo "APPVEYOR_REPO_BRANCH = $APPVEYOR_REPO_BRANCH"
echo "TRAVIS_BRANCH        = $TRAVIS_BRANCH"
echo "CI_BRANCH            = $CI_BRANCH"
echo "=============================================================="
echo "Build information:"
echo "=============================================================="















if [[ "$com" = "$tag" ]]&&[[ $dirty = $tagname ]]
then
    echo "Tagged commit: $tagname"
    PEP440="^([1-9]*!)?(0|[1-9]*)(\.(0|[1-9]*))*((a|b|rc)(0|[1-9]*))?(\.post(0|[1-9]*))?(\.dev(0|[1-9]*))?$"
    if [[ $tagname =~ $PEP440 ]]
    then
        echo "Git tag is a canonical PEP440 release version number"
        echo "deploy a setuptools package to pypi."
        echo "deploy conda packages to anaconda.org in the BjornFJohansson channel"
        tagged_commit=true
    else
        echo "Git tag is *NOT* a canonical PEP440 release version number"
        echo "git tag ($tagname) was not recognized"
        echo "or"
        echo "$dirty != $tagname"
        exit 0
    fi
elif [[ $msg = *"skip"* ]]
then
    echo "'skip' found in commit msg: '$msg'"
    echo "tests and builds skipped."
    echo "=============================================================="
    exit 0 
else
    echo "'skip' not found in commit msg: '$msg'"
    echo "but commit not tagged or tag dirty"
    echo "test suit will be run."
    tagged_commit=false

unset VIRTUAL_ENV
fi















if [[ $CI = true ]]||[[ $CI = True ]]
then
    echo "====================Running on CI server======================"
    [[ ! -z "$TWINE_USERNAME" ]] && echo "TWINE_USERNAME is set" || echo "TWINE_USERNAME is empty"
    [[ ! -z "$TWINE_PASSWORD" ]] && echo "TWINE_PASSWORD is set" || echo "TWINE_PASSWORD is empty"
    [[ ! -z "$ANACONDATOKEN"  ]] && echo "ANACONDATOKEN is set"  || echo "ANACONDATOKEN is empty"
    if [[ $TRAVIS = true ]]
    then
        branch=$TRAVIS_BRANCH
        echo "Running on TRAVIS, download Miniconda for MacOSX"
        miniconda="wget -q http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $APPVEYOR = true ]]||[[ $APPVEYOR = True ]]
    then
        branch=$APPVEYOR_REPO_BRANCH
        echo "Running on APPVEYOR, use installed Miniconda for Windows"
        miniconda="source appveyor_source_file.sh"
    elif [[ $CIRCLECI = true ]]
    then
        branch=$CI_BRANCH
        miniconda="wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $CI_NAME = codeship ]]
    then
        miniconda="wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        echo "CIRCLECI = $CIRCLECI"
        echo "CI_NAME  = $CI_NAME"
        exit 1
    fi
    echo "execute: $miniconda"
    $miniconda
    if [[ -f Miniconda_latest.sh ]]
    then
        bash Miniconda_latest.sh -b -p $HOME/miniconda
        export PATH="$HOME/miniconda/bin:$PATH"
        rm Miniconda_latest.sh
    fi
    conda config --set always_yes yes --set show_channel_urls yes --set anaconda_upload no
    conda update -yq conda
    conda update -yq pip
    conda config --append channels conda-forge 
    conda config --append channels BjornFJohansson
else
    echo "Not running on CI server, probably running on local computer"
    branch=$(git rev-parse --symbolic-full-name --abbrev-ref HEAD)
    local_computer=true
fi











echo "====================create conda environments================="
conda env create -f python36.yml
conda env create -f python37.yml












if [[ $tagged_commit = true ]]
then
    echo "========build conda package and setuptools package(s)========="
    conda install -yq -c anaconda conda-build
    conda install -yq -c anaconda conda-verify
    echo "conda-build version used:"
    conda-build -V
    rm -rf dist
    rm -rf build
    rm -rf tests/htmlcov
    pth2="$(conda build . --output --py 3.6)"
    pth3="$(conda build . --output --py 3.7)"
    echo "========build path(s)========================================="
    echo $pth2
    echo $pth3
    source activate python36
    conda build --python 3.6 --no-include-recipe --dirty .
    source activate python37
    conda build --python 3.7 --no-include-recipe --dirty .
    echo "========conda upload(s)======================================="
    if [[ $CI = true ]]||[[ $CI = True ]]
    then
        echo "========anaconda upload using ANACONDATOKEN====================="
        anaconda -t $ANACONDATOKEN upload $pth2 --label $condalabel --force
        anaconda -t $ANACONDATOKEN upload $pth3 --label $condalabel --force
    else
        echo "========anaconda upload========================================="
        anaconda upload $pth2 --label $condalabel --force
        anaconda upload $pth3 --label $condalabel --force
    fi
    if [[ $TRAVIS = true ]] # MacOSX on Travis
    then
        source activate python36
        python setup.py build bdist_wheel
        source activate python37
        python setup.py build bdist_wheel
    elif [[ $APPVEYOR = true ]]||[[ $APPVEYOR = True ]] # Windows on appveyor
    then
        source activate python36
        python setup.py build bdist_wheel 
        source activate python37
        python setup.py build bdist_wheel
        appveyor PushArtifact dist/*        
    elif [[ $CI_NAME = codeship ]]  # Linux on codeship
    then
        source activate python36
        python setup.py build bdist_wheel 
        source activate python37
        python setup.py build bdist_wheel 
        twine upload dist/*.whl --skip-existing
    elif [[ $local_computer = true ]]
    then
        echo "Local linux: python setup.py sdist --formats=zip bdist_wheel"
        source activate python36
        python setup.py build bdist_wheel 
        source activate python37
        python setup.py build bdist_wheel 
        twine upload dist/*.whl --skip-existing
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        echo "CIRCLECI = $CIRCLECI"
        exit 1
    fi
else
    echo "====================tests================="
    source activate python36
    python run_test.py
    source activate python37
    python run_test.py
fi
