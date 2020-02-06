#!/usr/bin/env python


fn = raw_input("file? ")
fp = open(fn, 'r')
text = fp.readline()
#xxx = fp.readline().split(" ")
#yyy = [int(e) for e in xxx]
#k, L, t = yyy
# print k, L, t

k = 4
L = 30
t = 3

lText = len(text)

counts = dict()

for i in xrange(0,L-k):
	pattern = text[i:i+k]
	counts[pattern] = counts.get(pattern, 0) + 1

# print counts
numEntries = len(counts)
for i in xrange(0, numEntries):
	if(counts.values()[i] >= t):
		print counts.keys()[i]

for i in xrange(0, lText-L):
	if(text[i:i+k] != text[i+L-k:i+L]):
		pattern = text[i:i+k]
		counts[pattern] -= 1
		pattern = text[i+L-k:i+L]
		counts[pattern] = counts.get(pattern, 0) + 1
# 		print counts
		if (counts[pattern] == t):
			print pattern

	
