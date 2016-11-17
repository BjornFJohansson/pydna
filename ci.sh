echo -e "Establish git variables:\n=============================================================="
tagname="$(git describe --abbrev=0 --tags)"
tag="$(git rev-list $tagname | head -n 1)"
com="$(git rev-parse HEAD)"
dirty=$(git describe --tags --dirty --always)
msg=$(git log -1 --pretty=%B)

branch=${TRAVIS_PULL_REQUEST_BRANCH:-$TRAVIS_BRANCH}

echo "Branch              : $branch"

echo "Current commit hash : $com"
echo "Dirty tag           : $dirty"
echo "Commit msg          : $msg"
echo "=============================================================="
echo "CI                   = $CI"
echo "DRONE                = $DRONE"
echo "TRAVIS               = $TRAVIS"
echo "APPVEYOR             = $APPVEYOR"
echo "=============================================================="

echo "APPVEYOR_REPO_BRANCH = $APPVEYOR_REPO_BRANCH"

echo "$DRONE_COMMIT: the commit hash currently being built"
echo "$DRONE_BRANCH: the branch currently being built"
echo "$BUILD_ID: the current build number"
echo "$GIT_COMMIT: the commit hash currently being built"
echo "$GIT_BRANCH: the branch currently being built"
