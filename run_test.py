#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import pytest

def main():

    os.environ["pydna_data_dir"] = tempfile.mkdtemp(prefix="pydna_data_dir_")
    os.environ["pydna_log_dir"] = tempfile.mkdtemp(prefix="pydna_log_dir_")
    os.environ["pydna_config_dir"] = tempfile.mkdtemp(prefix="pydna_config_dir_")
    os.environ["pydna_loglevel"] = str(logging.DEBUG)

    args = [ "tests",
             "src",
             "--capture=no",
             "--durations=10",
             "--cov=pydna",
             "--cov-report=html",
             "--cov-report=xml",
             "--import-mode=importlib",
             "--nbval",
             "--current-env",
             "--doctest-modules",
             "--capture=no",
             "--import-mode=importlib",
             "-vvv"]

    result_suite = pytest.cmdline.main(args)

    return result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
