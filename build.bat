
conda config --add channels ulmo
conda config --add channels synthicity

"%PYTHON%" setup.py install
if errorlevel 1 exit 1
