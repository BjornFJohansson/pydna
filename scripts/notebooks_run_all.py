#!/usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import os
import pathlib
import sys

try:
    from termcolor import colored
except ImportError:

    def colored(text, color):
        return text


thisdir = os.path.dirname(os.path.realpath(__file__))

exit_ = 0
exits = 0
notebooks_w_errors = []

for dirpath, dirnames, filenames in os.walk(thisdir):
    if dirpath.startswith((".", "_")):
        continue
    for file_ in filenames:
        if not file_.endswith(".ipynb"):
            continue
        pth = pathlib.Path(os.path.join(dirpath, file_))
        if any([f for f in pth.parts if f.startswith(("_", "."))]):
            continue
        os.chdir(dirpath)
        cmd = ["jupyter", "nbconvert", "--execute", "--inplace", file_]
        print(colored(dirpath, "blue"), colored(file_, "green"), end="")
        prc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        exit_ = prc.returncode
        if exit_ == 0:
            print(" OK")
        else:
            print(prc.stderr.decode("utf-8"))
            notebooks_w_errors.append(pth)
            exits += exit_
print("run_all.py finished with exit code {}".format(exit_))

for nb in notebooks_w_errors:
    print(nb)

sys.exit(exits)
