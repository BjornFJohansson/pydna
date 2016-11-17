tagname="$(git describe --abbrev=0 --tags)"
tag="$(git rev-list $tagname | head -n 1)"
com="$(git rev-parse HEAD)"
dirty=$(git describe --tags --dirty --always)
msg=$(git log -1 --pretty=%B)
GITBRANCH="$(git rev-parse --abbrev-ref HEAD)"
TRAVISBRANCH=${TRAVIS_PULL_REQUEST_BRANCH:-$TRAVIS_BRANCH}
branch=${DRONE_BRANCH:-$TRAVISBRANCH}
branch=${branch:-$APPVEYOR_REPO_BRANCH}
branch=${branch:-$GITBRANCH}

echo "TRAVIS_PULL_REQUEST_BRANCH $TRAVIS_PULL_REQUEST_BRANCH"
echo "DRONE_BRANCH $DRONE_BRANCH"
echo "TRAVIS_BRANCH $TRAVIS_BRANCH"
echo "APPVEYOR_REPO_BRANCH $APPVEYOR_REPO_BRANCH"
echo "GITBRANCH $GITBRANCH"




echo "=============================================================="
echo "Establish git variables:"
echo "=============================================================="
echo "Branch              : $branch"
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
echo "=============================================================="
echo "Build information:"
echo "=============================================================="
if [[ "$com" = "$tag" ]]
then
    echo "Tagged commit      : $tagname"
    tagged_commit=true
    re_final="^[0-9]\.[0-9]\.[0-9]$"
    re_alpha="^[0-9]\.[0-9]\.[0-9]a[0-9]+$"
    if [[ $tagname =~  $re_final ]]&&[[ $branch = "py3" ]]&&[[ $dirty = $tagname ]]
    then
        echo "Release tag and branch indicate Final release"
        echo "deploy to pypi and anaconda.org with label 'main'."
        echo "This is only done from the py3 branch"
        pypiserver="pypi"
        condalabel="main"
    elif  [[ $tagname =~ $re_alpha ]]&&[[ $branch = "py3dev" ]]&&[[ $dirty = $tagname ]]
    then
        echo "Release tag and branch indicate Alpha release"
        echo "deploy to testpypi and anaconda.org with label 'test'."
        echo "This is only done from the py3dev branch"
        pypiserver="testpypi"
        condalabel="test"
    else
        echo "Build cancelled because"
        echo "Release tag ($tagname) was not recognized"
        echo "or"
        echo "branch ($branch) was not py3 or py3dev"
        echo "or"
        echo "$dirty != $tagname"
        exit 1
    fi
else
    echo "Commit not tagged"
    echo "Run test suite, No build or install"
    tagged_commit=false
fi
echo "=============================================================="
