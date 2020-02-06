#! /usr/bin/env python

def HammingDist(pattern1, pattern2):
	l = len(pattern1)

	dist = 0	
	for i in xrange(0,l):
		if(pattern1[i] != pattern2[i]):
			dist += 1

	return dist

def CountPattern(text, pattern, distance=0):
	lText = len(text)
	lPattern = len(pattern)

	count = 0 
	for i in range (0, lText - lPattern + 1) :
        	p = text[i:i+lPattern] 
        	d = HammingDist(p, pattern) 
        	if (d <= distance) :
			print i,
                	count = count + 1
	return count

fn = raw_input('Name of file? :')

fp = open(fn, 'r')
text = fp.readline().rstrip()
pattern = fp.readline().rstrip()
dist = fp.readline().rstrip()
fp.close()

lText = len(text)
lPattern = len(pattern)
lDist = len(dist)

count = CountPattern(text, pattern, int(dist))
print 
print 'Count: ', count
