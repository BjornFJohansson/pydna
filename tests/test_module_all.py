#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

def test_repr():
    from pydna import all
    assert all.__all__==["Anneal", 
                         "pcr", 
                         "Assembly", 
                         "genbank",
                         "Genbank",
                         "download_text",  
                         "Dseqrecord",
                         "Dseq",
                         "read", 
                         "read_primer",  
                         "parse", 
                         "parse_primers", 
                         "ape", 
                         "primer_design", 
                         "assembly_fragments",
                         "circular_assembly_fragments",
                         "eq",
                         "gbtext_clean", 
                         "list_primers"]

if __name__ == '__main__':
    pytest.main([__file__, "-v", "-s", "--cov=pydna","--cov-report=html"])
