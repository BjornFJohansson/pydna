#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, tempfile
import nose

def main():

    if os.getenv("DRONE") or os.getenv("CI") or os.getenv("APPVEYOR"):
        print("\n\nTests run on continuous integration server\n\n")
        print("cwd", os.path.join(os.getcwd()))
        print()
        
        os.environ["pydna_data_dir"]    = os.path.join(os.getcwd(),"DATA")
        os.environ["pydna_log_dir"]     = os.path.join(os.getcwd(),"DATA")
        os.environ["pydna_config_dir"]  = os.path.join(os.getcwd(),"DATA")

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
        print()

    else:
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
        args = []
    else:
        del coveralls
        args = ["--with-coverage", "--cover-package=pydna", "--cover-erase"]    

    nose.run(argv=["--verbosity=3", 
                   "--nocapture", 
                   "--with-doctest", 
                   "--doctest-options=+ELLIPSIS"]+args)

    #print(("cache files", os.listdir( os.environ["pydna_data_dir"] )))

    print("                  _               _            _               _ _          __ _       _     _              _ _ ")
    print("                 | |             | |          | |             (_) |        / _(_)     (_)   | |            | | |")
    print("  _ __  _   _  __| |_ __   __ _  | |_ ___  ___| |_   ___ _   _ _| |_ ___  | |_ _ _ __  _ ___| |__   ___  __| | |")
    print(" | '_ \| | | |/ _` | '_ \ / _` | | __/ _ \/ __| __| / __| | | | | __/ _ \ |  _| | '_ \| / __| '_ \ / _ \/ _` | |")
    print(" | |_) | |_| | (_| | | | | (_| | | ||  __/\__ \ |_  \__ \ |_| | | ||  __/ | | | | | | | \__ \ | | |  __/ (_| |_|")
    print(" | .__/ \__, |\__,_|_| |_|\__,_|  \__\___||___/\__| |___/\__,_|_|\__\___| |_| |_|_| |_|_|___/_| |_|\___|\__,_(_)")
    print(" | |     __/ | ")
    print(" |_|    |___/  ")
    print("")

if __name__ == '__main__':
    print("script executed!!!!")
    main()
