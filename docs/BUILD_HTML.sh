#!/bin/bash

sphinx-build -b html . _build/html

echo `basename $0`
