#!/usr/bin/env bash
echo "=============================================================="
echo "BASH_VERSION" $BASH_VERSION
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
echo "=============================================================="
echo "Build information:"
echo "=============================================================="
if [[ "$com" = "$tag" ]]&&[[ $dirty = $tagname ]]
then
    echo "Tagged commit      : $tagname"
    tagged_commit=true
    re_final="^[0-9]\.[0-9]\.[0-9]$"
    re_alpha="^[0-9]\.[0-9]\.[0-9]a[0-9]+$"
    if [[ $tagname =~  $re_final ]]
    then
        echo "Tag indicate Final release"
        echo "deploy to pypi and anaconda.org with label 'main'."
        pypiserver="pypi"
        condalabel="main"
    elif  [[ $tagname =~ $re_alpha ]]
    then
        echo "Release tag indicate Alpha release"
        echo "deploy to testpypi and anaconda.org with label 'test'."
        pypiserver="testpypi"
        condalabel="test"
    else
        echo "Build cancelled because"
        echo "Release tag ($tagname) was not recognized"
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
echo "=============================================================="
if [[ $CI = true ]]||[[ $CI = True ]]
then
    echo "Running on CI server"
    echo "Creating a .pypirc file for setuptools"
    echo "[server-login]
    username: $pypiusername
    password: $pypipassword

    [distutils]
    index-servers=
        pypi
        testpypi

    [testpypi]
    repository = https://test.upload.pypi.org/legacy/
    username = $pypiusername
    password = $pypipassword

    [pypi]
    repository = https://pypi.python.org/pypi
    username = $pypiusername
    password = $pypipassword" > $HOME/.pypirc

    if [[ $TRAVIS = true ]]
    then
        echo "Running on TRAVIS, download Miniconda for MacOSX"
        miniconda="wget -q http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $APPVEYOR = true ]]||[[ $APPVEYOR = True ]]
    then
        echo "Running on APPVEYOR, use installed Miniconda for Windows"
        miniconda="source appveyor_source_file.sh"
        #miniconda="wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Windows-x86_64.exe -O Miniconda_latest.sh"
    elif [[ $CIRCLECI = true ]]
    then
        miniconda="wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        echo "CIRCLECI = $CIRCLECI"
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
    conda update -yq conda
    conda update -yq pip
    conda install conda-verify -yq
    conda install jinja2 -yq
    conda config --add channels BjornFJohansson
else
    echo "Not running on CI server, probably running on local computer"
fi
if [[ $tagged_commit = true ]]
then
    echo "build conda package and setuptools package(s)"
    conda install -yq conda-build
    conda-build -V
    conda create -yq -n pydnapipbuild   python=3.5 anaconda-client urllib3 twine pypandoc pandoc
    conda create -yq -n pydnacondabuild python=3.5 anaconda-client pypandoc pandoc nbval
    rm -rf dist
    rm -rf build
    rm -rf tests/htmlcov
    source activate pydnacondabuild
    which python
    python --version
    pth="$(conda build . --output)"
    echo $pth
    #conda info -a
    conda build .
    if [[ $CI = true ]]||[[ $CI = True ]]
    then
        anaconda -t $TOKEN upload $pth --label $condalabel --force
    else
        anaconda upload $pth --label $condalabel --force
    fi
    source activate pydnapipbuild
    conda upgrade -yq pip

    if [[ $TRAVIS = true ]]
    then
        python setup.py build bdist_wheel
        twine upload -r $pypiserver dist/pydna*.whl --skip-existing
    elif [[ $APPVEYOR = true ]]||[[ $APPVEYOR = True ]]
    then
        python setup.py build bdist_wininst bdist_msi
        twine upload -r $pypiserver dist/pydna*.exe --skip-existing
        twine upload -r $pypiserver dist/pydna*.msi --skip-existing
        appveyor PushArtifact dist/*
    elif [[ $CIRCLECI = true ]] # Linux
    then
        python setup.py sdist --formats=zip
        twine upload -r $pypiserver dist/pydna*.zip --skip-existing
    elif [[ $(uname) = "Linux" ]]
    then
        echo "Local linux: python setup.py sdist --formats=gztar,zip bdist_wheel"
        python setup.py sdist --formats=gztar,zip bdist_wheel
        twine upload -r $pypiserver dist/pydna*.zip --skip-existing
        twine upload -r $pypiserver dist/pydna*.gz  --skip-existing
        twine upload -r $pypiserver dist/pydna*.whl --skip-existing
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        echo "CIRCLECI = $CIRCLECI"
        exit 1
    fi
    ls dist
else
    echo "create test environment"
    conda env create -f test_environment.yml -q
    source activate testenv
    which python
    python --version
    python run_test.py
fi
