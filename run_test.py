#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import logging
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
        os.environ["pydna_data_dir"]    = tempfile.mkdtemp(prefix="pydna_data_dir_")
        os.environ["pydna_log_dir"]     = tempfile.mkdtemp(prefix="pydna_log_dir_")
        os.environ["pydna_config_dir"]  = tempfile.mkdtemp(prefix="pydna_config_dir_")
        os.environ["pydna_loglevel"]    = str( logging.DEBUG )
        
        
    print("                 _             ")
    print("                | |            ")
    print(" _ __  _   _  __| |_ __   __ _ ")
    print("| '_ \| | | |/ _` | '_ \ / _` |")
    print("| |_) | |_| | (_| | | | | (_| |")
    print("| .__/ \__, |\__,_|_| |_|\__,_|")
    print("| |     __/ |     ")
    print("|_|    |___/      ")
    print(" _            _   ")
    print("| |          | |  ")
    print("| |_ ___  ___| |_ ")
    print("| __/ _ \/ __| __|")
    print("| ||  __/\__ \ |_ ")
    print(" \__\___||___/\__|")
    print("           _ _       _ ")
    print("          (_) |     | |")
    print(" ___ _   _ _| |_ ___| |")
    print("/ __| | | | | __/ _ \ |")
    print("\__ \ |_| | | ||  __/_|")
    print("|___/\__,_|_|\__\___(_)")
    print("                       ")
    print("                       ")
    
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

    print("  __ _       _     _              _ _ ")
    print(" / _(_)     (_)   | |            | | |")
    print("| |_ _ _ __  _ ___| |__   ___  __| | |")
    print("|  _| | '_ \| / __| '_ \ / _ \/ _` | |")
    print("| | | | | | | \__ \ | | |  __/ (_| |_|")
    print("|_| |_|_| |_|_|___/_| |_|\___|\__,_(_)")


if __name__ == '__main__':
    main()
