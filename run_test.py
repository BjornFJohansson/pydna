#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import pytest

# https://discuss.python.org/t/testing-doctests-breaks-src-layout/3728
# https://coverage.readthedocs.io/en/coverage-5.2.1/cmd.html


def main():

    os.environ["pydna_data_dir"] = tempfile.mkdtemp(
        prefix="pydna_data_dir_")

    os.environ["pydna_log_dir"] = tempfile.mkdtemp(
        prefix="pydna_log_dir_")

    os.environ["pydna_config_dir"] = tempfile.mkdtemp(
        prefix="pydna_config_dir_")

    os.environ["pydna_loglevel"] = str(logging.DEBUG)

    args = [
        "src",  # doctestdir
        "--cov=pydna",
        "--cov-report=html",
        "--cov-report=xml",
        "--capture=no",
        "--durations=10",
        "--import-mode=importlib",
        "--nbval",
        "--current-env",
        "--doctest-modules",
        "--capture=no",
        "--import-mode=importlib",
        "-vvv",
    ]

    result_suite_src = pytest.main(args)

    args = [
        "tests",  # test suite
        "--cov=pydna",
        "--cov-append",
        "--cov-report=html",
        "--cov-report=xml",
        "--capture=no",
        "--durations=10",
        "--import-mode=importlib",
        "--nbval",
        "--current-env",
        "--doctest-modules",
        "--capture=no",
        "--import-mode=importlib",
        "-vvv",
    ]

    result_suite_tests = pytest.main(args)

    return result_suite_tests and result_suite_src


if __name__ == "__main__":
    result = main()
    sys.exit(result)
