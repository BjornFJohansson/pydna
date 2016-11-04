"""
GenBank FlatFile Parser

Parses GenBank FlatFile (GFF) into a JSON representation that's quickly usable in other scripts / etc.

Copyright (c) 2011 Anselm Levskaya (http://anselmlevskaya.com)
Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php) license.  
"""
import re
import json
from optparse import OptionParser
from pyparsing import *

#for debugging, print location and return token unaltered
def printf(s,l,t):
    print(l,": ",t)
    return t

#===============================================================================
# GenBank Grammar
#===============================================================================

#===============================================================================
# GenBank LOCUS Entry Parser

# LOCUS string, the first line of any genbank-sh file
# unlike many shite parsers, this should work for NCBI, ApE, and NTI style gb files
GoodLocus =    Literal("LOCUS") + \
               Word(alphas+nums+'-_().'+'\\').setResultsName("name") + \
               Word(nums)+Suppress(Literal('bp')) + \
               Word(alphas+'-').setResultsName("moleculetype") + \
               (CaselessLiteral("linear")|CaselessLiteral("circular")).setResultsName("topology") + \
               Suppress(Optional(Word(alphas))) + \
               Word(alphas+nums+'-').setResultsName("date")

# Older versions of ApE don't include a LOCUS name! Need separate def for this case:
BrokenLocus =  Literal("LOCUS").setResultsName("name") + \
               Word(nums)+Suppress(Literal('bp')) + \
               Word(alphas+'-').setResultsName("moleculetype") + \
               (CaselessLiteral("linear")|CaselessLiteral("circular")).setResultsName("topology") + \
               Suppress(Optional(Word(alphas))) + \
               Word(alphas+nums+'-').setResultsName("date")

LocusEntry = (GoodLocus|BrokenLocus)

GoodLocus =  ( Literal("LOCUS") +
               Word(alphas+nums+'-_().'+'\\').setResultsName("name") +
               Word(nums).setResultsName("size")+Suppress(CaselessLiteral('bp')) +
               Word(alphas+'-').setResultsName("seqtype") +
               (CaselessLiteral("linear")|CaselessLiteral("circular")).setResultsName("topology") +
               Optional(Word(alphas)).setResultsName("divcode") +
               Word(alphas+nums+'-').setResultsName("date") )

# Older versions of ApE don't include a LOCUS name! Need separate def for this case:
BrokenLocus1 =( Literal("LOCUS").setResultsName("name") +
                Word(nums).setResultsName("size")+Suppress(CaselessLiteral('bp')) +
                Word(alphas+'-').setResultsName("seqtype") +
                (CaselessLiteral("linear")|CaselessLiteral("circular")).setResultsName("topology") +
                Optional(Word(alphas)).setResultsName("divcode") +
                Word(alphas+nums+'-').setResultsName("date") )
            
# LOCUS       YEplac181	5741 bp 	DNA	SYN
BrokenLocus2 =( Literal("LOCUS") +
                Word(alphas+nums+'-_().'+'\\').setResultsName("name") +
                Word(nums).setResultsName("size")+Suppress(CaselessLiteral('bp')) +
                Word(alphas+'-').setResultsName("seqtype") +
                Optional(CaselessLiteral("linear")|CaselessLiteral("circular")).setResultsName("topology") +
                Optional(Word(alphas)).setResultsName("divcode") )

LocusEntry = (GoodLocus|BrokenLocus1|BrokenLocus2)

#===============================================================================
# Generic Entry

# this catches everything but the FEATURES and SEQUENCE entries, really should add parsing code for
# ACCESSION, COMMENTS, REFERENCE, ORGANISM, etc.
# (Though these entries are generally useless when it comes to hacking on DNA)

# All entries in a genbank file headed by an all-caps title with no space between start-of-line and title
CapWord = Word("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
# after titled line, all subsequent lines have to have at least one space in front of them
# this is how we split up the genbank record
SpacedLine =  White(min=1) + CharsNotIn("\n") + LineEnd()
#HeaderLine = CapWord + CharsNotIn("\n") + LineEnd()
GenericEntry =  Group(CapWord + Combine(CharsNotIn("\n") + LineEnd() +\
                             ZeroOrMore( SpacedLine ))).setResultsName("generics",listAllMatches=True)


#===============================================================================
# Definition Entry
#SuppressedSpacedLine =  Suppress(White(min=1)) + CharsNotIn("\n") + LineEnd()
#DefinitionEntry =  Suppress(Literal("DEFINITION")) + Combine(CharsNotIn("\n") + LineEnd() + ZeroOrMore( SuppressedSpacedLine ))

#===============================================================================
# GenBank Feature Table Parser

#==== Genbank Location String Parser
#
# a string of slices w. functional modifiers that go at most two levels deep
# single position is just a number i.e. 23423
# slice is N1..N2  w. N1<N2
# i.e.
# 23..88  --> seq[23:89] in python syntax (genbank uses inclusive slicing)
# 234..555 
# complement(234..434) --> rc(seq[234:435])
# join(23..343,454..666,777..999) --> seq[23:344]+seq[454:667]+seq[777:1000]
# complement(join(23..343,454..666,777..999))
# join(complement(34..123),complement(333..565))
#
# additionally the slices can have ambiguous locs like <454..999 or 232..>331
# also note the dumb 34.38 fuzzy slice notation
# i.e. <45..900  says the real feature starts "somewhere" before pos 45
#       45.48 says feature somewhere between bases 45->48
# lot of weird annotations best avoided because they connote ~useless knowledge for synthetic design
#
# if you don't know where something is, don't use it or guess and move on

LPAREN = Suppress("(")
RPAREN = Suppress(")")
SEP = Suppress(Literal(".."))

#recognize numbers w. < & > uncertainty specs, then strip the <> chars to make it fixed
gbIndex = Word(nums+"<>").setParseAction(lambda s,l,t: int(t[0].replace("<","").replace(">","")) )
SimpleSlice=Group(gbIndex + SEP + gbIndex) | Group(gbIndex).setParseAction(lambda s,l,t: [[t[0][0],t[0][0]]])

#recursive def for nested function syntax:  f( g(), g() )
complexSlice = Forward()
complexSlice << (Literal("complement") | Literal("join")) + LPAREN + ( delimitedList(complexSlice) | delimitedList(SimpleSlice) ) + RPAREN 
featLocation = Group( SimpleSlice | complexSlice)

def parseGBLoc(s,l,t):
    """retwingles parsed genbank location strings, assumes no joins of RC and FWD sequences """
    strand = 1
    locationlist=[]
    
    #see if there are any complement operators
    for entry in t[0]:
        if entry=="complement": strand=-1

    for entry in t[0]:
        if type(entry)!=type("string"):
            locationlist.append([entry[0],entry[1]])
            
    #return locationlist and strand spec
    return [ ['location', locationlist ], ['strand',strand] ]

featLocation.setParseAction(parseGBLoc)

#==== Genbank Feature Key-Value Pairs
def strip_multiline(s,l,t):
    whitespace=re.compile("[\n]{1}[ ]+")
    return whitespace.sub(" ",t[0])
def toInt(s,l,t):
    return int(t[0])

# Quoted KeyVal:   /key="value"
QuoteFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') + QuotedString('"',multiline=True).setParseAction( strip_multiline ) )

# UnQuoted KeyVal: /key=value  (I'm assuming it doesn't do multilines this way? wrong! ApE does store long labels this way! sigh.)
#NoQuoteFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') + OneOrMore(CharsNotIn("\n")) )
keyvalspacedline = White(exact=21)+CharsNotIn("/")+OneOrMore(CharsNotIn("\n"))+LineEnd()
NoQuoteFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') + Combine(CharsNotIn("\n") + LineEnd() + ZeroOrMore( keyvalspacedline )))

# Special Case for Numerical Vals:  /bases=12  OR  /bases="12"
NumFeaturekeyval = Group( Suppress('/') + Word(alphas+nums+"_-") + Suppress('=') +\
                         (Suppress("\"") + Word(nums).setParseAction(toInt) + Suppress("\"") ) | \
                         (Word(nums).setParseAction(toInt) ) \
                         )

# Key Only KeyVal: /pseudo
# post-parse convert it into a pair to resemble the structure of the first three cases i.e. [pseudo, True]
FlagFeaturekeyval = Group(Suppress('/') + Word(alphas+nums+"_-")).setParseAction(lambda s,l,t: [[t[0][0],True]] )

Feature = Group( Word(alphas+nums+"_-").setParseAction(lambda s,l,t: [ ["type", t[0]] ] ) +\
                 featLocation.setResultsName("location") +\
                 OneOrMore( NumFeaturekeyval | QuoteFeaturekeyval | NoQuoteFeaturekeyval | FlagFeaturekeyval ) \
                 )

FeaturesEntry = Literal("FEATURES") + Literal("Location/Qualifiers") + Group(OneOrMore(Feature)).setResultsName("features")

#===============================================================================
# GenBank Sequence Parser

# sequence is just a column-spaced big table of dna nucleotides
# should it recognize full IUPAC alphabet?  NCBI uses n for unknown region
Sequence = OneOrMore(Suppress(Word(nums)) + OneOrMore(Word("ACGTacgtNn")))

# Group(  ) hides the setResultsName names def'd inside, such that one needs to first access this group and then access the dict of contents inside
SequenceEntry = Suppress(Literal("ORIGIN")) + Sequence.setParseAction(lambda s,l,t: "".join(t) ).setResultsName("sequence")


#===============================================================================
# Final GenBank Parser

#GB files with multiple records split by "//" sequence at beginning of line
GBEnd = Literal("//")

#Begin w. LOCUS, slurp all entries, then stop at the end!
GB = LocusEntry + OneOrMore(FeaturesEntry | SequenceEntry | GenericEntry) + GBEnd

#NCBI often returns sets of GB files
multipleGB = OneOrMore(Group(GB))

#===============================================================================
# End Genbank Parser
#===============================================================================


#===============================================================================
# Main JSON Conversion Routine

def strip_indent(str):
    whitespace=re.compile("[\n]{1}(COMMENT){0,1}[ ]+")
    return whitespace.sub("\n",str)

def concat_dict(dlist):
    """more or less dict(list of string pairs) but merges
    vals with the same keys so no duplicates occur
    """
    newdict={}
    for e in dlist:
        if e[0] in newdict.keys():
            newdict[e[0]]=(newdict[e[0]]+strip_indent(e[1]))
        else:
            newdict[e[0]]=strip_indent(e[1])
    return newdict

def toJSON(gbkstring):
    parsed = multipleGB.parseString(gbkstring)
    
    jseqlist=[]
    for seq in parsed:
        #Print to STDOUT some details (useful for long multi-record parses)
        print(seq['name'], ":  length:", len(seq['sequence']) , " #features:" , len(seq['features'].asList()))
        
        #build JSON object
        jseq = { "__format__" : "jseq v0.1",
                 "name" : seq["name"],
                 "type" : seq["seqtype"],
                 "date" : seq["date"],
                 "topology" : seq["topology"],
                 "sequence" : seq["sequence"],
                 "features" : map(dict,seq['features'].asList()),
                 "annotations" : concat_dict(seq['generics'])
                 }
        jseqlist.append(jseq)

    return jseqlist

# command-line json conversion
if __name__ == "__main__":
    #parse command line string
    usage = "usage: %prog [options] gbfile_in jseqfile_out"
    parser = OptionParser(usage)
    #parser.add_option("-o", "--output", dest="output",
    #                  help="output json file")
    (options, args) = parser.parse_args()

    #Load the GBK file and Build a JSON object out of the parse tree

    if len(args)>0:
        infile = open(args[0],'r').read()

        jseqlist=toJSON(infile)
        
        #Output to new JSON file
        if len(args)>1:
            outfile=open(args[1],'w')
        else:
            outfile=open(args[0].replace(".gb",".json"),'w')
        
        if len(jseqlist)>1:
            json.dump(jseqlist,outfile)
        else:
            json.dump(jseqlist[0],outfile)

