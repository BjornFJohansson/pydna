pydna
===============================

A Custom Jupyter Widget Library

Installation
------------

To install use pip:

    $ pip install reprwidget
    $ jupyter nbextension enable --py --sys-prefix reprwidget


For a development installation (requires npm),

    $ git clone https://github.com/MEC/pydna.git
    $ cd pydna
    $ pip install -e .
    $ jupyter nbextension install --py --symlink --sys-prefix reprwidget
    $ jupyter nbextension enable --py --sys-prefix reprwidget
