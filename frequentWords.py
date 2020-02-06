#! /usr/bin/env python

def HammingDist(pattern1, pattern2):
	l = len(pattern1)

	dist = 0	
	for i in xrange(0,l):
		if(pattern1[i] != pattern2[i]):
			dist += 1

	return dist

def CountPattern(text, pattern, distance):
	lText = len(text)
	lPattern = len(pattern)

	count = 0 
	for i in range (0, lText - lPattern + 1) :
        	p = text[i:i+lPattern] 
        	d = HammingDist(p, pattern) 
        	if (d <= distance) :
                	count = count + 1
	return count

def FrequentWords(text, k, d):
	n = {0:'A', 1:'G', 2:'T', 3:'C'}

	lText = len(text)
	max = 0
	count = {}

	for i in xrange(4**k):
		x = i
		pattern = ''	
		for j in xrange(k):
			r = x % 4
			x = x / 4
			pattern = n[r] + pattern 
		
		count[i] = CountPattern(text, pattern, d)
		if (count[i] > max):
			max = count[i]

	for i in xrange (4**k) :
		if (count[i] == max) :
			x = i
			pattern = ''	
			for j in xrange(k):
				r = x % 4
				x = x / 4
				pattern = n[r] + pattern 
			print pattern,

fn = raw_input('Name of file? :')

fp = open(fn, 'r')
text = fp.readline().rstrip()
x = fp.readline().rstrip().split(" ")
y = [int(e) for e in x]
k, d = y
fp.close()

#print text
#print k, d
words = FrequentWords(text, k, d)
