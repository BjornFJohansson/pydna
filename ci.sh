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
        echo "Release tag and branch indicate Final release"
        echo "deploy to pypi and anaconda.org with label 'main'."
        pypiserver="pypi"
        condalabel="main"
    elif  [[ $tagname =~ $re_alpha ]]
    then
        echo "Release tag and branch indicate Alpha release"
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
