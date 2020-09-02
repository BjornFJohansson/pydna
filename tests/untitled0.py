#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sun Aug 30 19:39:16 2020 @author: bjorn
## In[] <-remove one # 

import pytest

class TestClass:
  
    def setup_class(cls):
        pass
        
    def test_buttons(self, data):
        # self.$attribute can be used, but not cls.$attribute?  
        pass
        print(456)
        
    def test_buttons2(self, data):
        # self.$attribute can be used, but not cls.$attribute?
        pass
        
    def teardown_class(cls):
        pass
    
    
if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s", "--cov=pydna","--cov-report=html"])
