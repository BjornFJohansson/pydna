#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import tempfile
import pytest

def main():

    if os.getenv("CI"):
        ci = os.getenv("DRONE") or os.getenv("TRAVIS") or os.getenv("APPVEYOR")
        print("Tests run on continuous integration server {}".format(ci))
        cwd = os.getcwd()
        print("current working directroy:", cwd)
        
        os.environ["pydna_data_dir"]    = os.path.join(cwd,"DATA")
        os.environ["pydna_log_dir"]     = os.path.join(cwd,"DATA")
        os.environ["pydna_config_dir"]  = os.path.join(cwd,"DATA")

        # create data directory if not present
        try:
            os.makedirs( os.environ["pydna_data_dir"] )
        except OSError:
            if os.path.isdir( os.environ["pydna_data_dir"] ):
                pass
            else:
                raise

        print('os.environ["pydna_data_dir"] = ',  os.environ["pydna_data_dir"])
        print('os.environ["pydna_log_dir"] = ',   os.environ["pydna_log_dir"])
        print('os.environ["pydna_config_dir"] = ',os.environ["pydna_config_dir"])
        print('')

    else:
        print("Tests run locally")
        os.environ["pydna_data_dir"] = tempfile.mkdtemp(prefix="pydna_data_dir_test_")
        os.environ["pydna_log_dir"]  = tempfile.mkdtemp(prefix="test_pydna_log_dir_test")

    print("      _             _                     _               _            _               _ _       _ ")
    print("     | |           | |                   | |             | |          | |             (_) |     | |")
    print("  ___| |_ __ _ _ __| |_   _ __  _   _  __| |_ __   __ _  | |_ ___  ___| |_   ___ _   _ _| |_ ___| |")
    print(" / __| __/ _` | '__| __| | '_ \| | | |/ _` | '_ \ / _` | | __/ _ \/ __| __| / __| | | | | __/ _ \ |")
    print(" \__ \ || (_| | |  | |_  | |_) | |_| | (_| | | | | (_| | | ||  __/\__ \ |_  \__ \ |_| | | ||  __/_|")
    print(" |___/\__\__,_|_|   \__| | .__/ \__, |\__,_|_| |_|\__,_|  \__\___||___/\__| |___/\__,_|_|\__\___(_)")
    print("                         | |     __/ |                                                             ")
    print("                         |_|    |___/                                                              ")

    try:
        import coveralls
    except ImportError:
        print("python-coveralls NOT installed!")
        args = []
    else:
        del coveralls
        args = ["--cov=pydna", "--cov-report=html", "--cov-report=xml"]    
    try:
        import nbval
    except ImportError:
        print("nbval NOT installed!")
    else:
        del nbval
        args.append("--nbval")   

    args = [".", "-v", "-s"] + args 
    cwd = os.getcwd()
    os.chdir("tests")
    pytest.cmdline.main(args)
    os.chdir(cwd)
    try:
        shutil.copy(os.path.join("tests","coverage.xml"), "coverage.xml")
    except FileNotFoundError:
        pass
    args = ["pydna", "--doctest-modules", "-v", "-s"]
    pytest.cmdline.main(args)

    print("                  _               _            _               _ _          __ _       _     _              _ _ ")
    print("                 | |             | |          | |             (_) |        / _(_)     (_)   | |            | | |")
    print("  _ __  _   _  __| |_ __   __ _  | |_ ___  ___| |_   ___ _   _ _| |_ ___  | |_ _ _ __  _ ___| |__   ___  __| | |")
    print(" | '_ \| | | |/ _` | '_ \ / _` | | __/ _ \/ __| __| / __| | | | | __/ _ \ |  _| | '_ \| / __| '_ \ / _ \/ _` | |")
    print(" | |_) | |_| | (_| | | | | (_| | | ||  __/\__ \ |_  \__ \ |_| | | ||  __/ | | | | | | | \__ \ | | |  __/ (_| |_|")
    print(" | .__/ \__, |\__,_|_| |_|\__,_|  \__\___||___/\__| |___/\__,_|_|\__\___| |_| |_|_| |_|_|___/_| |_|\___|\__,_(_)")
    print(" | |     __/ | ")
    print(" |_|    |___/  ")

if __name__ == '__main__':
    main()
