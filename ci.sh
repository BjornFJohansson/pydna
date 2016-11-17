echo -e "Establish git variables:\n=============================================================="
tagname="$(git describe --abbrev=0 --tags)"
tag="$(git rev-list $tagname | head -n 1)"
com="$(git rev-parse HEAD)"
dirty=$(git describe --tags --dirty --always)
msg=$(git log -1 --pretty=%B)
GITBRANCH="$(git rev-parse --abbrev-ref HEAD)"
TRAVISBRANCH=${TRAVIS_PULL_REQUEST_BRANCH:-$TRAVIS_BRANCH}
branch=${DRONE_BRANCH:-$TRAVISBRANCH$APPVEYOR_REPO_BRANCH$GITBRANCH}
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
