#!/bin/bash

rm -r _build

/home/bjorn/anaconda3/envs/bjorn37/bin/sphinx-build . _build/html

xdg-open _build/html/index.html

echo `basename $0`

read
