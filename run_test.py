#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import platform
import pytest

# pytest . --cov=pydna --cov-report=html --cov-report=xml --nbval --current-env  --capture=no --durations=10 --doctest-modules

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

    asciitext("tests python {}".format(platform.python_version()))

    try:
        import coveralls
    except ImportError:
        print("coveralls-python NOT installed!")
        args = []
    else:
        print("coveralls-python is installed!")
        del coveralls
        args = ["--cov=pydna", "--cov-report=html", "--cov-report=xml"]
    try:
        # adds functionality to py.test to recognise and collect
        # Jupyter notebooks
        import nbval
    except ImportError:
        print("nbval NOT installed!")
    else:
        print("nbval is installed!")
        del nbval
        args.append("--nbval")
        args.append("--current-env")

    mainargs = [".", "--capture=no", "--durations=10"] + args
    cwd = os.getcwd()
    os.chdir("tests")
    result_suite = pytest.cmdline.main(mainargs)
    os.chdir(cwd)

    asciitext("doctests python {}".format(platform.python_version()))
    doctestargs = ["src/pydna", "--doctest-modules", "--capture=no"]
    result_doctest = pytest.cmdline.main(doctestargs)

    asciitext("done!")

    return result_doctest and result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
