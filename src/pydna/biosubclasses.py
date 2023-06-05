#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 08:08:03 2023

@author: bjorn
"""
from copy import copy as _copy
from copy import deepcopy as _deepcopy
from Bio.SeqFeature import SeqFeature as _SeqFeature
from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from Bio.SeqFeature import ExactPosition as _ExactPosition


def empty_copy(obj):
    class Empty(obj.__class__):
        def __init__(self):
            pass

    newcopy = Empty()
    newcopy.__class__ = obj.__class__
    return newcopy


# https://stackoverflow.com/questions/57181829/deepcopy-override-clarification
# https://stackoverflow.com/questions/1500718/how-to-override-the-copy-deepcopy-operations-for-a-python-object
# https://stackoverflow.com/questions/24756712/deepcopy-is-extremely-slow


class ExactPosition(_ExactPosition):
    def __copy__(self):
        c = empty_copy(self)
        c.__dict__ = self.__dict__
        return c

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, _deepcopy(v, memo))
        return result

    def __copy__(self):
        c = empty_copy(self)
        c.__dict__ = self.__dict__
        return c

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, _deepcopy(v, memo))
        return result


class SimpleLocation(_SimpleLocation):
    def __copy__(self):
        c = empty_copy(self)
        c.__dict__ = self.__dict__
        return c

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, _deepcopy(v, memo))
        return result


class SeqFeature(_SeqFeature):
    def __copy__(self):
        return 123

    def __deecopy__(self):
        return 123
