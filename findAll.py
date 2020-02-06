#! /usr/bin/env python

fn = raw_input("file? ")
fp = open(fn, 'r')
pattern = 'CTTGATCAT' 
text = fp.readline()

lText = len(text) - 1
lPattern = len(pattern)
print lPattern

p = pattern[0:lPattern]

for i in xrange(0, lText-lPattern):
	p1 = text[i:i+lPattern] 
	if(p1 == p):
		print i,
