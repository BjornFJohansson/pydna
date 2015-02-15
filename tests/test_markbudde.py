#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test empty
'''

import unittest

#from pydna import ...

class test_empty(unittest.TestCase):

    def test_empty(self):
        ''' test mark budde'''
        import pydna           
        a = pydna.read('pGREG505.gb')
        self.assertTrue( a.name                 ,"pGREG505") 
        self.assertTrue( a.looped().name        ,"pGREG505")
        #self.assertTrue( a.annotations         ,"pGREG505")
        self.assertTrue( a.id                   ,"pGREG505")
        self.assertTrue( a.looped().id          ,"pGREG505")

        """
        >>> a = pydna.read('pGREG505.gb')
        >>> a.name
        'pGREG505'
        >>> a.looped().name
        'pGREG505'
        >>> a.annotations
        {'comment': '\nApEinfo:methylated:1', 'source': '', 'taxonomy': [], 'keywords': [''], 'accessions': ['pGREG505'], 'data_file_division': '   ', 'date': '15-DEC-2012', 'organism': '. .'}
        >>> a.id
        'pGREG505'
        >>> a.looped().id
        'pGREG505'
        >>> 
        """

if __name__ == '__main__':
    unittest.main()









