 #!/usr/bin/env bash


#git describe --exact-match HEAD

tagname="$(git describe --abbrev=0 --tags)"

tag="$(git rev-list $tagname | head -n 1)"

com="$(git rev-parse HEAD)"

if [ $com==$tag ]; then
   echo deploy
else
   echo $tag
   echo $com
fi
