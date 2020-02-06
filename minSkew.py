#!/usr/bin/env python


fn = raw_input("file? ")
fp = open(fn, 'r')
fp.readline()
text = fp.readline().rstrip()

lText = len(text)

skew = range(lText+1)
skew[0] = 0
min = lText

for i in xrange(0,lText):
	if (text[i] == 'G'):
		skew[i+1] = skew[i] + 1
	elif (text[i] == 'C'):
		skew[i+1] = skew[i] - 1
	else:
		skew[i+1] = skew[i]

	if (skew[i+1] < min):
		min = skew[i+1]

for i in xrange(1,lText+1):
	if skew[i] == min:
		print i,
