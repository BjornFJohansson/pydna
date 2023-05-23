#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import pytest


def main():
    """docstring."""

    os.environ["pydna_data_dir"] = tempfile.mkdtemp(prefix="pydna_data_dir_")
    os.environ["pydna_log_dir"] = tempfile.mkdtemp(prefix="pydna_log_dir_")
    os.environ["pydna_config_dir"] = tempfile.mkdtemp(prefix="pydna_config_dir_")
    os.environ["pydna_loglevel"] = str(logging.DEBUG)

    args = [
        "src",
        "tests",
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
        "-vvv",
        "--profile",  # profiling
    ]

    import pstats

    stats = pstats.Stats("./prof/combined.prof")
    stats.print_stats("pydna", 0.1)

    # Or alternatively
    # stats.print_stats("local_path", 20) # Only show 20 of the listings
    # stats.sort_stats('cumulative').print_stats('dir_name', 20) # Sort by cumulative time

    return int(pytest.main(args))


if __name__ == "__main__":
    sys.exit(main())
