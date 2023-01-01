#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
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
from PyQt5.QtGui import QTextCharFormat
from PyQt5.QtGui import QColor
from PyQt5.QtGui import QTextFormat
from PyQt5.QtCore import QSize


class SequenceWindow(QMainWindow):
    """docstring."""

    def __init__(self, text):
        QMainWindow.__init__(self)
        self.setMinimumSize(QSize(1000, 100))

        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

        self.setWindowTitle("sequence.name")
        centralWidget = QWidget(self)
        self.setCentralWidget(centralWidget)
        gridLayout = QGridLayout(self)
        centralWidget.setLayout(gridLayout)
        font = QFont("Consolas")
        font.setPointSize(16)
        centralWidget.setFont(font)

        textw = QPlainTextEdit()
        textw.appendPlainText(text)
        textw.appendPlainText(text)
        textw.appendPlainText(text)

        textw.setReadOnly(True)
        textw.document().findBlock(0)
        fmt = QTextCharFormat()
        fmt.setBackground(QColor("yellow"))
        gridLayout.addWidget(textw, 0, 0)

def run_app(text):
    app = QtWidgets.QApplication(sys.argv)
    mainWin = SequenceWindow(text)
    mainWin.show()
    app.exec_()


if __name__ == "__main__":
    run_app("GATC" * 5)













#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spyder python template."""

# import sys
# from PyQt5.QtWidgets import QApplication, QWidget, QLineEdit, QPlainTextEdit, QVBoxLayout
# from PyQt5.QtCore import QRegExp
# from PyQt5.QtGui import QColor, QRegExpValidator, QSyntaxHighlighter, QTextCharFormat

# class SyntaxHighlighter(QSyntaxHighlighter):
#     def __init__(self, parnet):
#         super().__init__(parnet)
#         self._highlight_lines = {}

#     def highlight_line(self, line_num, fmt):
#         if isinstance(line_num, int) and line_num >= 0 and isinstance(fmt, QTextCharFormat):
#             self._highlight_lines[line_num] = fmt
#             block = self.document().findBlockByLineNumber(line_num)
#             self.rehighlightBlock(block)

#     def clear_highlight(self):
#         self._highlight_lines = {}
#         self.rehighlight()

#     def highlightBlock(self, text):
#         blockNumber = self.currentBlock().blockNumber()
#         fmt = self._highlight_lines.get(blockNumber)
#         if fmt is not None:
#             self.setFormat(0, len(text), fmt)

# class AppDemo(QWidget):
#     def __init__(self):
#         super().__init__()
#         self.resize(1200, 800)

#         mainLayout = QVBoxLayout()

#         validator = QRegExpValidator(QRegExp(r'[0-9 ]+'))

#         self.lineEdit = QLineEdit()
#         self.lineEdit.setStyleSheet('font-size: 30px; height: 50px;')
#         self.lineEdit.setValidator(validator)
#         self.lineEdit.textChanged.connect(self.onTextChanged)
#         mainLayout.addWidget(self.lineEdit)

#         self.textEditor = QPlainTextEdit()
#         self.textEditor.setStyleSheet('font-size: 30px; color: green')
#         mainLayout.addWidget(self.textEditor)

#         for i in range(1, 21):
#             self.textEditor.appendPlainText('Line {0}'.format(i))

#         self.highlighter = SyntaxHighlighter(self.textEditor.document())
#         self.setLayout(mainLayout)

#     def onTextChanged(self, text):
#         fmt = QTextCharFormat()
#         fmt.setBackground(QColor('yellow'))

#         self.highlighter.clear_highlight()

#         try:
#             lineNumber = int(text) - 1
#             self.highlighter.highlight_line(lineNumber, fmt)
#         except ValueError:
#             pass

# app = QApplication(sys.argv)
# demo = AppDemo()
# demo.show()
# sys.exit(app.exec_())