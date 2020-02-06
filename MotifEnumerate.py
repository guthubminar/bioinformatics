#!/usr/bin/env python

from BioInfoUtils import *
 
# Parse an integer x into a dna string of lenth k
def int2dna(x, k) :
	i2n = 'ACGT'
	dna = ''
	for i in xrange(k) :
		r = x % 4
		x = x / 4
		dna = i2n[r] + dna
	return dna

# Given a k-mer, generate motif patterns also of length k, with Hamming distance at most 'dist'

def GenerateMotifCandidates(kmer, dist) :
	candidates = []
	k = len(kmer)
	for i in xrange(0, 4**k) :
		dna = int2dna(i, k)
		if HammingDist(dna, kmer) <= dist :
			candidates.append(dna)
	return candidates


# Check if a given k-mer appears in the string dna with at most dist mismatches

def isMotif(dna, kmer, dist) :
	l = len(dna)
	k = len(kmer)

	for i in xrange (l-k+1) :
		if HammingDist(dna[i:i+k], kmer) <= dist : return True
	return False

# Enumerate all k-mer motifs in a set of DNA strings

def EnumerateMotifs(DnaSet, k, d) :
	Motifs = []
	for dna in DnaSet : # for each dna string in DnaSet
		l = len(dna)
		for i in xrange(l-k+1) : 
			kmer = dna[i:i+k]  # for each k-mer in dna
			Candidates = GenerateMotifCandidates(kmer, d)
			""" print 'Motif candiates for ', kmer, ' : ', Candidates """

			for candidate in Candidates :
				success = 1
				for DNA in DnaSet :
					if not isMotif(DNA, candidate, d) : 
						success = 0
						break
				if (success == 1) & (candidate not in Motifs) : Motifs.append(candidate)	
        return Motifs

#----------------------------------------------------------------------

def EnumerateMotifs_main(fp) :
	x = fp.readline().rstrip().split()
	k, d = int(x[0]), int(x[1])
	print k, d

	DnaSet = []
	while True :
		x = fp.readline()
		if not x : break

		DnaSet.append(x.rstrip())

	print 'DNA Set: ', DnaSet 

	motifs = EnumerateMotifs(DnaSet, k, d)
	for m in motifs :
		print m,
#------------------------------------------------------------------------
fn = raw_input('File : ')
fp = open(fn, 'r')

EnumerateMotifs_main(fp)
	
