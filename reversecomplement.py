#! /usr/bin/env python

def ReverseComplement(text):
	lText = len(text)
	comp = {}
	comp['A'] = 'T'
	comp['T'] = 'A'
	comp['G'] = 'C'
	comp['C'] = 'G'

	revComp = []

	for i in xrange (0, lText - 1):
		revComp.append(comp[text[i]])

	revComp.reverse()
	print ''.join(revComp)


fn = raw_input("File?:")
fp = open(fn, 'r')
text = fp.readline()
print text[len(text)-1]
ReverseComplement(text)
