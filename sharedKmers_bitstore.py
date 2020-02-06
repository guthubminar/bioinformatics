#!/usr/bin/env python
# Author: Rajendran Panda, Jan 1, 2015

from BitStore import *

Compliment = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
agtc2int = {'A':0, 'G':1, 'T':2, 'C':3}
int2agtc = 'AGTC'

def int2kmer(x, k) :
	kmer = ''
	for j in xrange(k):
		r = x % 4
		x = x >> 2 
		kmer = int2agtc[r] + kmer
	return kmer

def kmer2int(kmer, k) :
	rslt = 0
	for i in range(k) :
		rslt = (rslt>>2) + agtc2int[kmer[i]]
	return rslt
		

# This function returns a bit field object of length= |text|
# Bits are set to 1 in positions where the pattern occurs in text
def FindOccurrences(text, pattern) :
	lText = len(text)
	k = len(pattern)
	bs = BitStore()
	for i in xrange(lText-k+1) :
		if text[i:i+k] == pattern : bs[i] = 1
	return bs
	
# This function returns occurrences for all 4**k k-mers in text
def Occurrence(text, k):
	BS = []		# occurrence bits storage

        lText = len(text)
        max = 0
        count = {}

        for i in xrange(4**k):
		kmer = int2kmer(i, k)
                bs = FindOccurrences(text, kmer)
		BS.append(bs)

	return BS

def SharedKmers(text1, text2, k) :
	# Break a k-mer into k1 number of 5-mers and a k2-mer, where k = 5*k1 + k2
	# This is to reduce the O(|text1| * |text2|) complexity to generating occurrences
	# of 5-mers and k2-mers in text2, then breaking every k-mer of text1 (and its
	# reverseCompliment) into 5-mers and a k2-mer and looking up the occurrence tables
	# to determine whether the k-mer is shared by the texts

	k1 = k/5
	k2 = k%5

	BS1 = Occurrence(text2, 5)
	if k2 > 0 :
		BS2 = Occurrence(text2, k2)

	# Find reverse compliment of text1
	t11 = list(text1)
	t11 = [Compliment[x] for x in t11]
	t11 = reversed(t11)
	t11 = ''.join(t11)

	m = len(text1); n = len(text2)

	for i in range(m - k + 1) :
		kmer = text1[i:i+k]
		rcKmer = t11[m-k-i:m-i]

		if k1 > 0 :
			subKmer = kmer[0:5]
			rcSubKmer = rcKmer[0:5]

			index1 = kmer2int(subKmer, 5)
			index2 = kmer2int(subKmer, 5)

			bs1 = BS1[index1].copy()
			bs2 = BS1[index2].copy()


			for j in range(1, k1) :
				subKmer = kmer[j*5 : j*5 + 5]
				rcSubKmer = rcKmer[j*5 : j*5 + 5]
				index1 = kmer2int(subKmer, 5)
				index2 = kmer2int(rcSubKmer, 5)
				tmp1 = BS1[index1].copy()
				tmp2 = BS1[index2].copy()

				bs1 = bs1 & (tmp1 >> (j*5))
				bs2 = bs2 & (tmp2 >> (j*5))


		if k2 > 0 :
			subKmer = kmer[k-k2:]
			rcSubKmer = rcKmer[k-k2:]
			index1 = kmer2int(subKmer, k2)
			index2 = kmer2int(rcSubKmer, k2)

			tmp1 = BS2[index1].copy()
			tmp2 = BS2[index2].copy()

			if k1 > 0 :
				bs1 = bs1 & (tmp1 >> (k1*5))
				bs2 = bs2 & (tmp2 >> (k1*5))
			else :
				bs1 = tmp1 
				bs2 = tmp2 

				
		for j in range(n) :
			if bs1[j] == 1 or bs2[j] == 1 : 
				coord = '(' + str(i) + ', ' + str(j) + ')'
				print coord

def SharedKmersBruteForce(text1, text2, k) :
	m = len(text1); n = len(text2)
	yloc = {} # hash table to store y locations for a kmer
	for i in range(m-k+1):
		pattern = text1[i:i+k]
		if pattern in yloc.keys():
			for j in yloc[pattern] :
				coord = '(' + str(i) + ', ' + str(j) + ')'
				print coord
		else:
			yloc[pattern] = []
			for j in range(n-k+1):
				if text2[j:j+k] == pattern :
					coord = '(' + str(i) + ', ' + str(j) + ')'
					print coord
					yloc[pattern].append(j) 
#------------------------------------------	 
fn = raw_input('File ? ')
fp = open(fn, 'r')

k = int(fp.readline().rstrip())
s1 = fp.readline().rstrip()
s2 = fp.readline().rstrip()
SharedKmersBruteForce(s2, s1, k)

