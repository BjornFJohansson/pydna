#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, re, shutil, pathlib

#for f in [f for f in os.listdir(".") if f.startswith("test_")]:
#    shutil.copy(f, str(pathlib.Path(f).with_suffix(".py_")))


for f in [f for f in os.listdir(".") if f.endswith(".py_")]:
    if not "test_parse.py_" in f:
        continue
    
    print(f)
    
    with open(f,"r") as file_:
        lines = file_.readlines()
    
    newlines=[]
    newline=""
    flag=False
    
    for line in lines:
        if "assert_equal" in line:    
            regex = "(\s+)assert_equal\((.+),(.+)\)"            
            m = re.search(regex, line)            
            newline = "{}assert {} == {}".format(m.groups()[0], m.groups()[1].strip(), m.groups()[2].strip())
            newlines.append(newline)
            flag=True
        else:
            newlines.append(line)
    if flag:        
        nf = str(pathlib.Path(f).with_suffix(".py"))
        print("".join(newlines))
        #import sys; sys.exit(42)
        with open(nf,"w+") as file_:
            file_.writelines(newlines)
    else:
        shutil.copy(f, str(pathlib.Path(f).with_suffix(".py")))