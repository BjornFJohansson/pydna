bash build.sh
rem "%PYTHON%" setup.py build
rem "%PYTHON%" setup.py install
if errorlevel 1 exit 1
