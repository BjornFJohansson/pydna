#!/usr/bin/env python
# -*- coding: utf-8 -*-

import PyPDF2
import sys
import subprocess
import re
import time

pdf = "SynthesisReport_3931076.pdf"

pdf = PyPDF2.PdfFileReader(open(pdf, "rb"))

text = u""

for i in range(0, pdf.getNumPages()):
    print str(i)
    extractedText = pdf.getPage(i).extractText()
    text += extractedText

print text

"""

first primer after

QC
Report
1

then
name
sequence (ends with number indicating length)
PicCht3.FWD  <-name
cggctGAATTCATTAATGCT    <-seqeunce
AGAAGTAAC (29)     <- length

second primer and following after "MALDI"

MALDI
2
PicCht3.FWD
cggctGAATTCATTAATGCT
AGAAGTAAC (29)

"""


sys.exit(42)

if __name__ == "__main__":
    pass

    # Number 65027
    # Name Bárbara Filipa Cerqueira Bernardino
    # Email A65027@alunos.uminho.pt
    # Course Engenharia Biológica
    # Enrolled SAUM No

    # regex = u"Number(.+?)Name(.+?)Email(.+?)Course(.+?)Enrolled SAUM(Yes|No)"
    # regex = u"Número(.+?)Nome(.+?)Email(.+?)Curso(.+?)Inscrito SAUM(Sim|Não)"

    match = re.findall(regex, text)

    with open("alunos.txt", "w") as f:
        for m in match:
            f.write(u"{}   {}\n".format(m[0], m[1]).encode("utf8"))

    from pyparsing import Word, Literal, printables, LineStart, SkipTo, Combine

    name = Word(printables).setResultsName("name")
    seq_start = Literal("5'").suppress()
    seq_stop = Literal("3'").suppress()
    sequence = Combine(seq_start + SkipTo(seq_stop)).setResultsName("seq")
    mwg_primer = LineStart() + name + SkipTo(LineStart()) + sequence

    result = mwg_primer.scanString(raw_string)

    seqlist = [data for data, dataStart, dataEnd in result]

    number += len(seqlist)

    fasta_string = ""

    for data in seqlist:
        number -= 1
        s = data.seq.strip("-").replace("\n", "").replace(" ", "")
        fasta_string += ">{number}_{name} ({length}-mer)\n{seq}\n\n".format(
            number=number, name=data.name, length=len(s), seq=s
        )

    fasta_string
