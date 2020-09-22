#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import platform
import pytest
import pathlib
import pkg_resources

# pytest . --cov=pydna --cov-report=html --cov-report=xml --nbval
# --current-env  --capture=no --durations=10 --doctest-modules

#
try:
    from pyfiglet import Figlet
except ImportError:
    asciitext = print
else:

    def asciitext(*args, **kwargs):
        f = Figlet(font="doom")
        print(f.renderText(" ".join(args)), **kwargs)


def main():

    os.environ["pydna_data_dir"] = tempfile.mkdtemp(prefix="pydna_data_dir_")
    os.environ["pydna_log_dir"] = tempfile.mkdtemp(prefix="pydna_log_dir_")
    os.environ["pydna_config_dir"] = tempfile.mkdtemp(prefix="pydna_config_dir_")
    os.environ["pydna_loglevel"] = str(logging.DEBUG)

    asciitext("tests py {}".format(platform.python_version()))

    installed = {pkg.key for pkg in pkg_resources.working_set}

    args = []

    if "coveralls" in installed:
        print("coveralls-python is installed.")

        args = ["--cov=pydna",
                "--cov-report=html",
                "--cov-report=xml",
                "--import-mode=importlib"]
    else:
        print("coveralls-python NOT installed! (pip install coveralls)")

    if "nbval" in installed:
        print("nbval is installed.")
        args.append("--nbval")
        args.append("--current-env")
    else:
        print("nbval NOT installed! (pip install nbval)")

    mainargs = ["tests", "--capture=no", "--durations=10"] + args

    result_suite = pytest.cmdline.main(mainargs)


    from pydna import __file__ as pydnainit

    doctestdir = str(pathlib.Path(pydnainit).parent)

    asciitext("doctests py {}".format(platform.python_version()))
    doctestargs = [
        doctestdir,
        "--doctest-modules",
        "--capture=no",
        "--import-mode=importlib",
    ]
    result_doctest = pytest.cmdline.main(doctestargs)

    asciitext("done!")

    return result_doctest and result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
