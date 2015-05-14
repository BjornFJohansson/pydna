#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import nose

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
os.environ["pydna_cache"]  = "nocache"

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
