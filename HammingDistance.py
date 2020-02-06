#!/usr/bin/env python


fn = raw_input("file? ")
fp = open(fn, 'r')
text1 = fp.readline()
text2 = fp.readline()

lText = len(text1)

dist = 0

for i in xrange(0,lText):
	if (text1[i] != text2[i]):
		dist += 1
print dist
