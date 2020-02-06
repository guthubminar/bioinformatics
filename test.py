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
		x = x / 4 
		kmer = int2agtc[r] + kmer
	return kmer

def kmer2int(kmer, k) :
	rslt = 0
	for i in range(k) :
		rslt = (rslt * 4) + agtc2int[kmer[i]]
	return rslt
		

def FindOccurrences(text, pattern) :
	lText = len(text)
	k = len(pattern)
	bs = []
	for i in xrange(lText-k+1) :
		if text[i:i+k] == pattern : bs.append(i)
	# print pattern, bs
	return bs
	
# This function returns occurrences for all 4**k k-mers in text
def Occurrence(text, k):
	BS = []		# occurrence storage

        lText = len(text)
        max = 0
        count = {}

        for i in xrange(4**k):
		kmer = int2kmer(i, k)
                bs = FindOccurrences(text, kmer)
		BS.append(bs)

	return BS

def SharedKmers(text1, text2, k, z) :
	# Break a k-mer into k1 number of z-mers and a k2-mer, where k = z*k1 + k2
	# This is to reduce the O(|text1| * |text2|) complexity to generating occurrences
	# of z-mers and k2-mers in text2, then breaking every k-mer of text1 (and its
	# reverseCompliment) into z-mers and a k2-mer and looking up the occurrence tables
	# to determine whether the k-mer is shared by the texts

	m = len(text1); n = len(text2)

	k1 = k/z; k2 = k%z

	if k1 > 0 : BS1 = Occurrence(text2, z)
	if k2 > 0 : BS2 = Occurrence(text2, k2)

	# print BS1; print BS2

	# Find reverse compliment of text1
	rcText1 = list(text1)
	rcText1 = [Compliment[x] for x in rcText1]
	rcText1 = reversed(rcText1)
	rcText1 = ''.join(rcText1)

	cache = {}

	for i in range(m - k + 1) :
		kmer = text1[i:i+k]
		if kmer in cache.keys() :
			for j in cache[kmer] :
				coord = '(' + str(i) + ', ' + str(j) + ')'
				print coord
			continue  # outer for loop

		rcKmer = rcText1[m-k-i:m-i]

		# print 'Start: ', kmer, rcKmer

		if k1 > 0 :
			subKmer = kmer[0:z]
			rcSubKmer = rcKmer[0:z]

			index1 = kmer2int(subKmer, z)
			index2 = kmer2int(rcSubKmer, z)

			bs1 = BS1[index1]
			bs2 = BS1[index2]

			# print bs1, bs2

			for j in range(1, k1) :
				if bs1 == [] : break
				subKmer = kmer[j*z : (j+1)*z]
				index1 = kmer2int(subKmer, z)
				tmp1 = BS1[index1]
				if tmp1 == [] : 
					bs1 = []
					break
				tmp1 = [x - j*z for x in tmp1]
				bs1 = list(set(bs1) & set(tmp1))

			for j in range(1, k1) :
				if bs2 == [] : break
				rcSubKmer = rcKmer[j*z : (j+1)*z]
				index2 = kmer2int(rcSubKmer, z)
				tmp2 = BS1[index2]
				if tmp2 == [] : 
					bs2 = []
					break
				tmp2 = [x - j*z for x in tmp2]
				bs2 = list(set(bs2) & set(tmp2))

		if k2 > 0 :
			subKmer = kmer[k-k2:]
			rcSubKmer = rcKmer[k-k2:]
			# print 'Check: ', subKmer, rcSubKmer
			index1 = kmer2int(subKmer, k2)
			index2 = kmer2int(rcSubKmer, k2)

			# print 'kmer & index: ', index1, subKmer
			# print 'kmer & index: ', index2, rcSubKmer

			if k1 > 0 :
				tmp1 = [x - k1*z for x in BS2[index1]]
				tmp2 = [x - k1*z for x in BS2[index2]]
				# print 'In : ', bs1, tmp1
				# print 'In : ', bs2, tmp2
				bs1 = list(set(bs1) & set(tmp1))
				bs2 = list(set(bs2) & set(tmp2))
				# print 'Out: ', bs1, bs2
			else :
				bs1 = BS2[index1]
				bs2 = BS2[index2]

		bs = bs1 + bs2
		bs = sorted(bs)
		cache[kmer] = bs
		# print i, kmer, bs
		# print

		for j in bs :
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
SharedKmers(s1, s2, k, 5)

