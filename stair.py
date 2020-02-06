#!/usr/bin/env python

def numpaths(m, n) :
	if m == 0 : return 1
	if n == 0 : return 1
	return numpaths(m-1, n) + numpaths(m, n-1)

x = numpaths(16,12)

print x
