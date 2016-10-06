#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "1.0.0"

print('''

This script automatically builds, converts and uploads 
conda packages.

It is meant to be placed in the same folder as the meta.yaml


''')



import subprocess
import os

cwd = os.getcwd()

script_dir = os.path.dirname(os.path.realpath(__file__))

os.chdir(script_dir)

bashCommand = "conda build . --output "
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
pth = process.communicate()[0]
print(bashCommand)
print(pth)
response = input("build (y/n) ? ")

if response.lower().startswith("y"):
    bashCommand = "conda build . --no-anaconda-upload"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
    output = process.communicate()[0]
    print(output)
else:
    print("build skipped.")

response = input("convert (y/n) ? ")

if response.lower().startswith("y"):
    bashCommand = "conda convert {} -p all --output-dir conda-bld".format(pth)
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
    output = process.communicate()[0]
    print(output)
else:
    print("convert skipped.")
    
os.chdir("conda-bld")

fn = os.path.basename(pth)

for directory in os.listdir("."):
    
    bldpth = os.path.join(directory, fn)  
    print(bldpth)
    response = input("upload (y/n) ? ")
    if response.lower().startswith("y"):
        bashCommand = "anaconda upload {}".format(bldpth)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, universal_newlines=True)
        output = process.communicate()[0]
        print(output)
    else:
        print("upload skipped.")   

os.chdir(cwd)

print("done!")
