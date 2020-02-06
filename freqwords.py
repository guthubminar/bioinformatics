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
	lText = len(text)
	max = 0
	count = {}
	words = []

	for i in xrange (0, lText - k + 1) :
		pattern = text[i:i+k]
		count[i] = CountPattern(text, pattern, d)
		if (count[i] > max):
			max = count[i]

	for i in xrange (0, lText - k + 1) :
		pattern = text[i:i+k]
		if (count[i] == max) :
			print count[i], pattern
		if ((count[i] == max) and not (pattern in words)):
			words.append(pattern)
	return words

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
print words
