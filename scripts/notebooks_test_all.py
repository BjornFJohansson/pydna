#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import os
import pathlib
import sys

try:
    from termcolor import colored
except ImportError:

    def colored(text, color):
        return text


exit_ = 0
exits = 0
notebooks_w_errors = []

thisdir = os.path.dirname(os.path.realpath(__file__))


for dirpath, dirnames, filenames in sorted(os.walk(thisdir)):
    if dirpath.startswith((".", "_")):
        continue
    for file_ in filenames:
        if not file_.endswith(".ipynb"):
            continue
        if file_ == "index.ipynb":
            continue
        p = pathlib.Path(os.path.join(dirpath, file_))
        if any([f for f in p.parts if f.startswith(("_", "."))]):
            continue
        os.chdir(dirpath)
        print("\n\n")
        print(colored(dirpath, "blue"), colored(file_, "red"))
        exit_ = pytest.cmdline.main(["--verbose", "--capture=no", "--nbval", file_])
        if exit_:
            notebooks_w_errors.append(p)
            exits += exit_

print("finished with exit code {}".format(exit_))

print("Notebooks with errors:")
for nb in notebooks_w_errors:
    print(nb)

sys.exit(exits)
