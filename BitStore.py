#!/usr/bin/env python
# Author: Rajendran Panda,   Jan 1, 2015
# Extended an on-line python receipe for bit field with more operators

class BitStore(object):
    def __init__(self,value=0):
        self._d = value

    def __getitem__(self, index):
        return (self._d >> index) & 1 

    def __setitem__(self,index,value):
        value    = (value&1L)<<index
        mask     = (1L)<<index
        self._d  = (self._d & ~mask) | value

    def __getslice__(self, start, end):
        mask = 2L**(end - start) -1
        return (self._d >> start) & mask

    def __setslice__(self, start, end, value):
        mask = 2L**(end - start) -1
        value = (value & mask) << start
        mask = mask << start
        self._d = (self._d & ~mask) | value
        return (self._d >> start) & mask

    def __rshift__(self, value):
	self._d = self._d >> value
	return self

    def __lshift__(self, value):
	self._d = self._d << value
	return self

    def __int__(self):
        return self._d

    def __and__(self, another):
	val = self._d & another._d
	return BitStore(val)

    def copy(self):
	return BitStore(self._d)
