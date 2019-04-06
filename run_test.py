#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import shutil
import logging
import tempfile
import platform
import pytest
try:
    from pyfiglet import Figlet
except ImportError: 
    asciitext = print
else:
    def asciitext(*args,**kwargs):
        f = Figlet(font='doom') 
        print(f.renderText(" ".join(args)), **kwargs)


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
    
    asciitext("test suite on python {}".format(platform.python_version()))
    
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
        import nbval
    except ImportError:
        print("nbval NOT installed!")
    else:
        print("nbval is installed!")
        del nbval
        args.append("--nbval") 
    
    mainargs = [".", "-vv", "-s", "--durations=10"] + args 
    cwd = os.getcwd()
    os.chdir("tests")
    result_suite = pytest.cmdline.main(mainargs)
    os.chdir(cwd)
    
    try:
        shutil.copy(os.path.join("tests","coverage.xml"), "coverage.xml")
        shutil.copy(os.path.join("tests",".coverage"),    ".coverage")
    except FileNotFoundError:
        pass

    asciitext("doc testson python {}".format(platform.python_version()))
    doctestargs = ["pydna", "--doctest-modules", "-vv", "-s"]
    result_doctest = pytest.cmdline.main(doctestargs)

    asciitext("done!")

    return result_doctest and result_suite

if __name__ == '__main__':
    result = main()
    sys.exit(result)
