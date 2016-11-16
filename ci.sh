echo "Establish git variables"
tagname="$(git describe --abbrev=0 --tags)"
tag="$(git rev-list $tagname | head -n 1)"
com="$(git rev-parse HEAD)"
branch="$(git rev-parse --abbrev-ref HEAD)"
echo "Branch             : $branch"
echo "Current commit hash: $com"
if [[ "$com" = "$tag" ]]
then
    echo "Tagged commit : $tagname"
    tagged_commit=true
    re_final="^[0-9]\.[0-9]\.[0-9]$"
    re_alpha="^[0-9]\.[0-9]\.[0-9]a[0-999]$"
    if [[ $tagname =~  $re_final ]]&&[[ "$branch" = "py3" ]]
    then
        echo -e "Release tag and branch indicate Final release\ndeploy to pypi and anaconda.org with label 'main'. \nThis is only done from the py3 branch"
        pypiserver="pypi"
        condalabel="main"
    elif  [[ $tagname =~ $re_alpha ]]&&[[ "$branch" = "py3dev" ]]
    then
        echo -e "Release tag and branch indicate Alpha release\ndeploy to testpypi and anaconda.org with label 'test'. \nThis is only done from the py3dev branch"
        pypiserver="testpypi"
        condalabel="test"
    else
        echo -e "Release tag ($tagname) was not recognized or branch ($branch) was not py3 or py3dev"
        exit 1
    fi
else
    echo "Commit not tagged"
    tagged_commit=false
fi
if [[ $CI = true ]]
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
        miniconda="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $TRAVIS = true ]]
    then
        miniconda="http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $APPVEYOR = true ]]
    then
        miniconda=""
        # Miniconda is installed by default on APPVEYOR
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "DRONE    = $DRONE"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        exit 1
    fi
    wget -q $miniconda
    if [[ -f Miniconda_latest.sh ]]
    then
        bash Miniconda_latest.sh -b -p $HOME/miniconda
        export PATH="$HOME/miniconda/bin:$PATH"
    fi
    conda update -yq conda
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels BjornFJohansson
else
    echo "Not running on CI server, probably running on local computer"
fi
if [[ $tagged_commit = true ]]
then
    echo "build conda package and setuptools package(s)"
    conda install -yq conda-build anaconda-client
    if [ "$branch" = "py2" ]
    then
        conda create -q -y -n pipbuild   python=2.7
        conda create -q -y -n pydnabuild python=2.7
    elif [ "$branch" = "py3" ]||[ "$branch" = "py3dev" ]
    then
        conda create -q -y -n pipbuild   python=3.5
        conda create -q -y -n pydnabuild python=3.5
    fi
    #conda info --envs
    source activate pydnabuild
    pth="$(conda build . --output)"
    echo $pth
    #conda info -a
    conda build .
    ###########################################################################anaconda -t $TOKEN upload $pth --label $condalabel --force
    source activate pipbuild
    conda upgrade -yq pip
    #pip install setuptools wheel twine
    pip install twine
    #conda install -y -q -c conda-forge pandoc=1.18
    #pandoc --from=markdown --to=rst --output=README.rst README.md
    #git add README.rst
    #git commit -m "processed README.md --> README.rst"
    #git tag -d $tagname
    #git tag $tagname
    if [[ $DRONE=true ]]
    then
        python setup.py build sdist --formats=gztar,zip bdist_wheel
    elif [[ $TRAVIS=true ]]
    then
        python setup.py build bdist_dmg
    elif [[ $APPVEYOR=true ]]
    then
        python setup.py build bdist_wininst
    elif [[ $(uname) = "Linux" ]]
    then
        python setup.py build sdist --formats=gztar,zip bdist_wheel
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "DRONE    = $DRONE"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        exit 1
    fi
    exit 1
    ##################################################################################twine upload -r $pypiserver dist/*
else
    echo "Commit not tagged"
    echo "No build or install, only run test suite"
    echo "create test environment"
    conda env create -f test_environment.yml -q
    source activate testenv
    python run_test.py
fi