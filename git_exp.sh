 #!/usr/bin/env bash


echo "[server-login]
username: $username
password: $password

[distutils]
index-servers=
    pypi
    pypitest

[pypitest]
repository = https://testpypi.python.org/pypi
username = $username
password = $password

[pypi]
repository = https://pypi.python.org/pypi
username = $username
password = $password" > .pypirc

tagname="$(git describe --abbrev=0 --tags)"

tag="$(git rev-list $tagname | head -n 1)"

com="$(git rev-parse HEAD)"


if [ "$com" = "$tag" ]
then
   echo deploy
   #echo $tagname
   echo $tag
   echo $com
   #pip install wheel
   #python setup.py sdist sdist --formats=gztar,zip bdist_wheel upload

else
   echo no deploy
   echo $tag
   echo $com
fi



python setup.py sdist upload -r pypitest

