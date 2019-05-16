#!/usr/bin/env bash
source /home/bjorn/anaconda3/bin/activate bjorn36
python --version # example way to see that your virtual env loaded
git branch
pytest . -vv -s --durations=10 -k test_module_ass
echo "press any key to close"
read -n1 slask 
