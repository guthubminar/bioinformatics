#!/usr/bin/env python

from BioInfoUtils import *

peptide = 'VKLFPWFNQY'
lPattern = 3 * len(peptide)

fn = raw_input("file? ")
fp = open(fn, 'r')

reminder = ''
count = 0
linenum = 0

while True:
	text = fp.readline().rstrip()

	if (text == '') :
		break

	linenum += 1
	text = reminder + text
#	print linenum, text
	lText = len(text)

	for i in xrange(lText - lPattern + 1):
		pattern = text[i:i+lPattern]
		revComp = ReverseComplement(pattern)
#		print pattern, revComp

		pat = transcribe(pattern)
#		print pat
		if (translate(pat) == peptide) :
			print pattern
			count += 1
			continue

		pat = transcribe(revComp)
#		print pat
		if (translate(pat) == peptide) :
			print pattern
			count += 1

	reminder = text[lText - lPattern + 1:lText]

print count
