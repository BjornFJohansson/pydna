#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""A DNA sequence widget."""

import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QGridLayout
from PyQt5.QtWidgets import QPlainTextEdit
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QDesktopWidget
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QSize
import os
from multiprocessing import Process


class SequenceWindow(QMainWindow):
    """docstring."""

    def __init__(self, record):
        QMainWindow.__init__(self)
        self.setMinimumSize(QSize(1000, 100))

        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

        self.setWindowTitle(record.name)
        centralWidget = QWidget()
        self.setCentralWidget(centralWidget)
        gridLayout = QGridLayout()
        centralWidget.setLayout(gridLayout)
        font = QFont("Consolas")
        font.setPointSize(16)
        centralWidget.setFont(font)
        textw = QPlainTextEdit()

        text = str(record.seq)
        html = ""
        for fl in record.features[0].location.parts:
            html = (f"{text[:fl.start]}<span style='background-color:#ffff00;'>"
                    f"{text[fl.start:fl.end]}</span>"
                    f"{text[fl.end:]}")

        textw.appendHtml(html)
        textw.setReadOnly(True)
        gridLayout.addWidget(textw, 0, 0)

def run_app(record):
    if os.fork() != 0:
        return
    app = QtWidgets.QApplication(sys.argv)
    app.mainWin = SequenceWindow(record)
    app.mainWin.show()
    app.exec_()


if __name__ == '__main__':

    t = """
        LOCUS       myDNA                     12 bp    DNA     linear       02-JAN-2023
        DEFINITION  .
        ACCESSION
        VERSION
        SOURCE      .
          ORGANISM  .
        COMMENT
        COMMENT     ApEinfo:methylated:0
        FEATURES             Location/Qualifiers
              misc_feature    5..8
                              /locus_tag="New Feature"
                              /label="New Feature"
                              /ApEinfo_label="New Feature"
                              /ApEinfo_fwdcolor="cyan"
                              /ApEinfo_revcolor="green"
                              /ApEinfo_graphicformat="arrow_data {{0 0.5 0 1 2 0 0 -1 0
                              -0.5} {} 0} width 5 offset 0"
        ORIGIN
                1 atccGATCgg tc
        //
        """
    from pydna.readers import read

    s = read(t)

    s.features

    p = Process(target=run_app, args=(s,))
    p.start()
    p.join()
