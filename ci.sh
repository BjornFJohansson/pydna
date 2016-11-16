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
else
    echo "Commit not tagged"
    tagged_commit=false
fi
printf "\n\n\n\n\n"
if [[ $CI=true ]]
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
    if [[ $DRONE=true ]]
    then
        $miniconda = "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $TRAVIS=true ]]
    then
        $miniconda = "http://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda_latest.sh"
    elif [[ $APPVEYOR=true ]]
    then
        $miniconda = ""
        # Miniconda is installed by default on APPVEYOR
    else
        echo "Running on CI server but none of the expected environment variables are set to true"
        echo "CI       = $CI"
        echo "DRONE    = $DRONE"
        echo "TRAVIS   = $TRAVIS"
        echo "APPVEYOR = $APPVEYOR"
        exit 1
    wget -q $miniconda
    if [[ -f Miniconda_latest.sh ]]
    then
        bash Miniconda_latest.sh -b -p $HOME/miniconda
        export PATH="$HOME/miniconda/bin:$PATH"
    fi
    conda update -yq conda
    conda config --add channels BjornFJohansson
    conda config --add channels conda-forge
else
    echo "Not running on CI server, probably running on local computer"
fi
if [[ $tagged_commit ]]
then
    echo "Tagged commit : $tagname"
    echo "build conda package and setuptools module distributions"
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
    source activate pipbuild
    conda install -y -q -c conda-forge pandoc=1.18
    pandoc --from=markdown --to=rst --output=README.rst README.md
    git add README.rst
    git commit -m "processed README.md -->README.rst"
    git tag $tagname
    



