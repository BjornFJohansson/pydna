#!/bin/bash

/home/bjorn/anaconda3/envs/bjorn39/bin/sphinx-apidoc -f -o . ../src/pydna

echo `basename $0`

read
