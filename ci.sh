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
echo "DRONE                = $DRONE"
echo "TRAVIS               = $TRAVIS"
echo "APPVEYOR             = $APPVEYOR"
echo "CIRCLECI             = $CIRCLECI"
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
        echo "Release tag indicate Final release"
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
        exit 1
    fi
else
    echo "Commit not tagged or tag dirty"
    echo "Run test suite, No build or install"
    tagged_commit=false
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
    repository = https://testpypi.python.org/pypi
    username = $pypiusername
    password = $pypipassword

    [pypi]
    repository = https://pypi.python.org/pypi
    username = $pypiusername
    password = $pypipassword" > $HOME/.pypirc
    if [[ $DRONE = true ]]
    then
        miniconda="wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $TRAVIS = true ]]
    then
        miniconda="wget -q http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $APPVEYOR = true ]]||[[ $APPVEYOR = True ]]
    then
        miniconda=""
        CONDA_ROOT="C:\Miniconda35x64"
        PATH="$CONDA_ROOT;$CONDA_ROOT/Scripts;$CONDA_ROOT/Library/bin;$PATH"
        # Miniconda is installed by default on APPVEYOR
    elif [[ $CIRCLECI = true ]]
    then
        miniconda="wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "DRONE    = $DRONE"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        echo "CIRCLECI = $CIRCLECI"
        exit 1
    fi
    $miniconda
    if [[ -f Miniconda_latest.sh ]]
    then
        bash Miniconda_latest.sh -b -p $HOME/miniconda
        export PATH="$HOME/miniconda/bin:$PATH"
    fi
    which python
    python --version
    conda update -yq conda
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels BjornFJohansson
else
    echo "Not running on CI server, probably running on local computer"
fi
#if [[ $(uname) = *"NT"* ]]
#then
#    source=""
#else
#    source=source
#fi
if [[ $tagged_commit = true ]]
then
    echo "build conda package and setuptools package(s)"
    conda install -yq conda-build
    conda create -q -y -n pydnapipbuild   python=3.5 anaconda-client
    conda create -q -y -n pydnacondabuild python=3.5 anaconda-client
    rm -rf dist
    rm -rf build
    rm -rf tests/htmlcov
    source activate pydnacondabuild
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
    conda install -yq urllib3 twine
    #conda install -y -q -c conda-forge pandoc=1.18
    #pandoc --from=markdown --to=rst --output=README.rst README.md
    #git add README.rst
    #git commit -m "processed README.md --> README.rst"
    #git tag -d $tagname
    #git tag $tagname
    if [[ $DRONE=true ]]
    then
        python setup.py sdist --formats=gztar,zip bdist_wheel
    elif [[ $TRAVIS=true ]]
    then
        python setup.py bdist_dmg
    elif [[ $APPVEYOR=true ]]||[[ $APPVEYOR=True ]]
    then
        python setup.py bdist_wininst
    elif [[ $CIRCLECI=true ]]
    then
        python setup.py sdist --formats=gztar,zip bdist_wheel
    elif [[ $(uname) = "Linux" ]]
    then
        python setup.py sdist --formats=gztar,zip bdist_wheel
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "DRONE    = $DRONE"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        echo "CIRCLECI = $CIRCLECI"
        exit 1
    fi
    twine upload -r $pypiserver dist/* --skip-existing

else
    echo "create test environment"
    conda env create -f test_environment.yml -q
    source activate testenv
    python run_test.py
fi




