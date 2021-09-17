#!/bin/bash

sphinx-apidoc -f -o . ../src/pydna

echo `basename $0`

read
