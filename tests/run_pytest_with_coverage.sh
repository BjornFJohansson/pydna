#!/usr/bin/env bash 
pytest . -v -s --cov=pydna --cov-report=html --cov-report=xml  #(in tests folder)
# coverage run -m py.test && coverage html  #(in tests folder)
read
