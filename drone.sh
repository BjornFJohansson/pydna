echo "Establish git variables"
tagname="$(git describe --abbrev=0 --tags)"
tag="$(git rev-list $tagname | head -n 1)"
com="$(git rev-parse HEAD)"
branch="$(git rev-parse --abbrev-ref HEAD)"

echo "Branch        : $branch"
echo "Current commit: $com"

if [[ $DRONE=true ]]
then
    echo "Running on DRONE"
    echo "Create a .pypirc file"
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
    # installing Miniconda
    wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh
    bash Miniconda_latest.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    conda update -yq conda
    conda config --add channels BjornFJohansson
    conda config --add channels conda-forge
else
    echo "Not running on DRONE"
fi

if [[ "$com" = "$tag" ]]
then
    echo "Tagged commit : $tagname"
    echo "build conda package and setuptools module distributions"
    conda install -yq conda-build anaconda-client
    if [ "$branch" = "py2" ]
    then
        conda create -q -y -n pydnabuild python=2.7
        conda create -q -y -n pipbuild   python=2.7
    elif [ "$branch" = "py3" ]||[ "$branch" = "py3dev" ]
    then
        conda create -q -y -n pydnabuild python=3.5
        conda create -q -y -n pipbuild   python=3.5
    fi
    source activate pydnabuild
    pth="$(conda build . --output)"
    echo $pth
    conda info -a
    conda build .
    source activate pipbuild
    #conda install -y -q -c conda-forge pandoc=1.18
    #pandoc --from=markdown --to=rst --output=README.rst README.md
    #git add README.rst
    #git commit -m "processed README.md -->"
    conda upgrade pip -yq
    pip install twine
    python setup.py build sdist --formats=gztar,zip bdist_wheel
    re_final="^[0-9]\.[0-9]\.[0-9]$"
    re_alpha="^[0-9]\.[0-9]\.[0-9]a[0-999]$"
    if [[ $tagname =~  $re_final ]]||[ "$branch" = "py3" ]
    then
        echo "Release tag indicate Final release"
        echo "deploy to pypi and anaconda.org"
        echo "this is only done from the py3 branch"
        anaconda -t $TOKEN upload $pth
        #python setup.py register
        twine upload dist/* --skip-existing        
    elif  [[ $tagname =~ $re_alpha ]]||[ "$branch" = "py3dev" ]
    then
        echo "Release tag indicate Alpha release"
        echo "deploy to testpypi and anaconda.org with label 'test'"
        echo "this is only done from the py3dev branch"
        anaconda -t $TOKEN upload $pth --label test --force
        #python setup.py register -r testpypi
        twine upload -r testpypi dist/*
    else
        echo "Release tag was not recognized"
        echo "or branch was not py3 or py3dev"
        echo "do nothing"
    fi
else
    echo "Commit not tagged"
    echo "No build or install, only run test suite"
    echo "create test environment"
    conda env create -f test_environment.yml -q
    source activate testenv
    python run_test.py
fi
