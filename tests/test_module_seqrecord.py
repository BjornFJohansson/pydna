import pytest

def test_add_feature():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as BSeqRecord
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.seqrecord  import SeqRecord
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(3,9)
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq="GGATCC")
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=Seq("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=Dseq("GGATCC"))    
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=SeqRecord("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=BSeqRecord("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    s.add_feature(seq=Dseqrecord("GGATCC"))
    assert s.extract_feature(0).seq == SeqRecord("GGATCC").seq
    s = SeqRecord("tttGGATCCaaa")
    with pytest.raises(TypeError):
        s.add_feature(seq=Dseqrecord("GGGGGG"))
    s = SeqRecord("tttATGaaaTAAggg")
    s.add_feature(3,12)
    assert s.features[0].qualifiers["label"] == ['orf9']

    
    from Bio.Seq import Seq 
    
    from pydna.seqrecord import SeqRecord 
    
    a=SeqRecord(Seq("atgtaa")) 
    
    a.add_feature(2,4) 
    
    assert a.list_features() == '+-----+---------------+-----+-----+-----+-----+------+------+\n| Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |\n+-----+---------------+-----+-----+-----+-----+------+------+\n|   0 | L:ft2         | --> | 2   | 4   |   2 | misc |  no  |\n+-----+---------------+-----+-----+-----+-----+------+------+'
    a.features[0].qualifiers
    del a.features[0].qualifiers["label"]
    assert a.list_features()=='+-----+---------------+-----+-----+-----+-----+------+------+\n| Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |\n+-----+---------------+-----+-----+-----+-----+------+------+\n|   0 | nd            | --> | 2   | 4   |   2 | misc |  no  |\n+-----+---------------+-----+-----+-----+-----+------+------+'
    a.features[0].qualifiers["note"]=["AwesomeFeature"]
    assert a.list_features()=='+-----+------------------+-----+-----+-----+-----+------+------+\n| Ft# | Label or Note    | Dir | Sta | End | Len | type | orf? |\n+-----+------------------+-----+-----+-----+-----+------+------+\n|   0 | N:AwesomeFeature | --> | 2   | 4   |   2 | misc |  no  |\n+-----+------------------+-----+-----+-----+-----+------+------+'
    
    
    
    

def test_stamp():
    from pydna.seqrecord  import SeqRecord
    from pydna.readers    import read
    from pydna.utils      import eq
    
    from Bio.Seq           import Seq
    from Bio.SeqRecord     import SeqRecord as Srec

    a=SeqRecord("attt")
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo" in a.definition
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    
    a=SeqRecord("attt")
    a.description = "something"
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo" in a.definition
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "something" in a.definition
    
    a=SeqRecord("attt")
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    a.seq = a.seq + "a"
    with pytest.raises(ValueError):
        a.stamp()

    a=SeqRecord(Seq("attt"))
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo" in a.definition
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    
    a=SeqRecord(Seq("attt"))
    a.description = "something"
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo" in a.definition
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert "something" in a.definition
    
    a=SeqRecord(Seq("attt"))
    assert a.stamp() == "SEGUID_ot6JPLeAeMmfztW1736Kc6DAqlo"
    a.seq = a.seq + "a"
    with pytest.raises(ValueError):
        a.stamp()
        
def test___hash__():
    from Bio.Seq import Seq
    from pydna.seqrecord  import SeqRecord

    s = SeqRecord(Seq("GGATCC"))
    t = SeqRecord(Seq("GGATCC"))
    u = SeqRecord(Seq("GGATCc"))
    assert hash(s) == hash(t) != hash(u)
    
    assert s==t
    assert s!=u                  
    assert s!="1"
    
    assert s<u
    assert not s>u
    
    with pytest.raises(TypeError):
        s>"3"
    with pytest.raises(TypeError):
        s<"3"


    
    
def test_olaps():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as BSeqRecord
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.seqrecord  import SeqRecord

    s = SeqRecord(Seq("GGATCC"))
    assert "GGATCC" == str(s.olaps("GGATCC", limit = 4)[0].seq)
    assert "GGATCC" == str(s.olaps(Seq("GGATCC"), limit = 4)[0].seq)
    assert "GGATCC" == str(s.olaps(BSeqRecord(Seq("GGATCC")), limit = 4)[0].seq)
    assert "GGATCC" == str(s.olaps(Dseq("GGATCC"), limit = 4)[0].seq)
    assert "GGATCC" == str(s.olaps(Dseqrecord(Dseq("GGATCC")), limit = 4)[0].seq)
    assert "GGATCC" == str(s.olaps(Dseqrecord("GGATCC"), limit = 4)[0].seq)


    
    
    
    
    
def test_format():
    from Bio.Seq import Seq
    from pydna.seqrecord  import SeqRecord

    s = SeqRecord(Seq("GGATCC"))
    s.format("gb")
    s.format("genbank")
    s.format("fasta")   

    
def test_seqrecord():
    import datetime
    import pydna
    from pydna import seqrecord, _PydnaWarning
    from Bio.Seq import Seq

    
    s = Seq("ATGAAATAA")
    
    obj = seqrecord.SeqRecord(s, name = "1234567890123456")
    
    with pytest.warns(None) as pdw:
        obj = seqrecord.SeqRecord(s, name = "12345678901234567")
    
    obj = seqrecord.SeqRecord(s, name = "<unknown name>")
    
    assert obj.name == "name"
    
    obj = seqrecord.SeqRecord(s, id   = "<unknown id>")
    
    assert obj.id == "id"
    
    obj = seqrecord.SeqRecord(s, description = "<unknown description>")
    
    assert obj.description == "description"
    
    obj = seqrecord.SeqRecord(s, annotations = {'date': '24-DEC-1970'})
    
    assert obj.annotations['date'] == '24-DEC-1970'
    
    
    obj = seqrecord.SeqRecord(s)
    
    assert obj.map_target == None
    
    assert obj.isorf()
    
    not_orf = Seq("TGAAATAAA")
    
    obj = seqrecord.SeqRecord(not_orf)
    
    assert not obj.isorf()
    
    s3 = Seq("aaaATGAAATAAttt")
    
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    
    loc = FeatureLocation(3,12)
    
    sf = SeqFeature(loc, type="misc", qualifiers={"label":"lbl"})
    
    obj = seqrecord.SeqRecord(s3, features = [sf])
    
    obj.add_colors_to_features_for_ape()
    
    assert obj.features[0].qualifiers['ApEinfo_fwdcolor'] == ['#66ffa3']
    assert obj.features[0].qualifiers['ApEinfo_revcolor'] == ['#66ffff']
    
    assert obj.gc() == 6.7
    
    repr(obj) == "SeqRecord(seq=Seq('aaaATGAAATAAttt'), id='id', name='name', description='description', dbxrefs=[])"
    
    obj.annotations = {'date': '24-DEC-1970'}
    
    st = str("ID: id\n"
              "Name: name\n"
              "Description: description\n"
              "Number of features: 1\n"
              "/date=24-DEC-1970\n"
              "Seq('aaaATGAAATAAttt', DNAAlphabet())")
    
    assert str(obj) == st
    
    
    obj.locus = "new123"
    assert obj.locus        == obj.name == "new123"
    obj.locus = "new456"
    assert obj.locus        == obj.name == "new456"
    
    obj.accession = "new123"
    assert obj.accession    == obj.id == "new123"
    obj.id = "new456"
    assert obj.accession    == obj.id == "new456"
    
    obj.definition = "new123"
    assert obj.definition    == obj.description == "new123"
    obj.description = "new456"
    assert obj.definition    == obj.description == "new456"
    
    
    with pytest.warns(None) as pdw:
        obj.locus = "12345678901234567"
    
    lf = str("+-----+---------------+-----+-----+-----+-----+------+------+\n"
             "| Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |\n"
             "+-----+---------------+-----+-----+-----+-----+------+------+\n"
             "|   0 | L:l b l       | --- | 3   | 12  |   9 | misc | yes  |\n"
             "+-----+---------------+-----+-----+-----+-----+------+------+")
    
    assert obj.list_features() == lf
    
    exft = obj.extract_feature(0)
    
    assert str(exft.seq) == 'ATGAAATAA'
    import textwrap
    gbf = textwrap.dedent('''LOCUS       1234567890123456          15 bp    DNA              UNK 24-DEC-1970
    DEFINITION  new456.
    ACCESSION   new456
    VERSION     new456
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
         misc            4..12
                         /label="lbl"
                         /ApEinfo_fwdcolor="#7f3f3f"
                         /ApEinfo_revcolor="#7f3f3f"
    ORIGIN
            1 aaaatgaaat aattt
    //'''.strip())
    
    #print()
    #print(gbf)
    
    #print(str(obj.format("gb")))
    
    #assert gbf+"\n" == str(obj.format("gb"))

    #print(obj.__hash__())

if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s","--cov=pydna","--cov-report=html"])




#    >>> import warnings
#    >>> from pydna import _PydnaWarning
#    >>> warnings.simplefilter('ignore', _PydnaWarning)

#obj.name = "1111111111111111"


#    def __init__(self, seq, id="<unknown id>", name="<unknown name>",
#                 description="<unknown description>", dbxrefs=None,
#                 features=None, annotations=None,
#                 letter_annotations=None):
#        """Create a SeqRecord.



