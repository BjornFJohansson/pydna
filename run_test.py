#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import nose
import appdirs
import shutil
import errno

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True

if os.getenv("DRONE") or os.getenv("CI") or os.getenv("APPVEYOR"):
    print "Continouos integration"
    os.environ["pydna_data_dir"] = os.path.join(os.getcwd(),"..","..","DATA")
else:
    os.environ["pydna_data_dir"] = appdirs.user_data_dir("pydna_test").encode(sys.getfilesystemencoding())
    try:
        shutil.rmtree( os.environ["pydna_data_dir"] )
    except OSError, e:
        if e.errno == errno.ENOENT:
            print "no cache to delete."
    print "cache deleted."

def main():
    print "      _             _                     _               _            _               _ _       _ "
    print "     | |           | |                   | |             | |          | |             (_) |     | |"
    print "  ___| |_ __ _ _ __| |_   _ __  _   _  __| |_ __   __ _  | |_ ___  ___| |_   ___ _   _ _| |_ ___| |"
    print " / __| __/ _` | '__| __| | '_ \| | | |/ _` | '_ \ / _` | | __/ _ \/ __| __| / __| | | | | __/ _ \ |"
    print " \__ \ || (_| | |  | |_  | |_) | |_| | (_| | | | | (_| | | ||  __/\__ \ |_  \__ \ |_| | | ||  __/_|"
    print " |___/\__\__,_|_|   \__| | .__/ \__, |\__,_|_| |_|\__,_|  \__\___||___/\__| |___/\__,_|_|\__\___(_)"
    print "                         | |     __/ |                                                             "
    print "                         |_|    |___/                                                              "
    print

    cwd = os.getcwd()
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(os.path.join(dname, "tests"))
    nose.run(argv=[__file__, "--all-modules", "--verbosity=3", "--nocapture", "--with-doctest", "--doctest-options=+ELLIPSIS"])
    os.chdir(os.path.join(dname,"pydna"))
    nose.run(argv=[__file__, "--all-modules", "--verbosity=3", "--nocapture", "--with-doctest", "--doctest-options=+ELLIPSIS"])
    os.chdir(cwd)

    print "cache files", os.listdir( os.environ["pydna_data_dir"] )

    print "                  _               _            _               _ _          __ _       _     _              _ _ "
    print "                 | |             | |          | |             (_) |        / _(_)     (_)   | |            | | |"
    print "  _ __  _   _  __| |_ __   __ _  | |_ ___  ___| |_   ___ _   _ _| |_ ___  | |_ _ _ __  _ ___| |__   ___  __| | |"
    print " | '_ \| | | |/ _` | '_ \ / _` | | __/ _ \/ __| __| / __| | | | | __/ _ \ |  _| | '_ \| / __| '_ \ / _ \/ _` | |"
    print " | |_) | |_| | (_| | | | | (_| | | ||  __/\__ \ |_  \__ \ |_| | | ||  __/ | | | | | | | \__ \ | | |  __/ (_| |_|"
    print " | .__/ \__, |\__,_|_| |_|\__,_|  \__\___||___/\__| |___/\__,_|_|\__\___| |_| |_|_| |_|_|___/_| |_|\___|\__,_(_)"
    print " | |     __/ | "
    print " |_|    |___/  "
    print

if __name__ == '__main__':
    main()
