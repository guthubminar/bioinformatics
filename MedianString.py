#!/usr/bin/env python

# Utilities for Finding Motifs 
# Author: Rajendran Panda
# Created: 11/28/2014

import random

# Find Hamming Distance between 2 equal length patterns 
def Distance(pattern1, pattern2):
	l = len(pattern1)

	dist = 0	
	for i in xrange(0,l):
		if(pattern1[i] != pattern2[i]):
			dist += 1

	""" print pattern1, ' --> ', pattern2, ' = ', dist """
	return dist

# Find distance between a pattern and a strand
def Pattern2StrandDistance(pattern, strand) :
	kp = len(pattern)
	ks = len(strand)

	minDist = kp
	index = 0

	for i in range (0, ks - kp + 1) :
		p = strand[i:i+kp]
		dist = Distance(pattern, p)
		if dist < minDist :
			minDist = dist
			index = i
	""" print pattern, ' --> ', strand, ' = ', minDist """
	return minDist, index

# Find distance between a pattern and a DNA
def Pattern2DnaDistance (pattern, dna) :
	Distance = 0
	Indices = []

	for x in dna :
		dist, index = Pattern2StrandDistance(pattern, x)
		Indices.append(index)
		Distance += dist
	""" print pattern, ' --> ', dna, ' = ', Distance """
	return Distance, Indices

# Find the k-mer Median String for a DNA
def MedianString(DNA, k) :
	int2nuc = 'ATGC'
	Distance = 2**20
	Indices = []
	winner = ''

        for i in xrange(4**k):
                x = i
                pattern = ''
                for j in xrange(k):
                        r = x % 4
                        x = x / 4
                        pattern = int2nuc[r] + pattern

		dist, indices = Pattern2DnaDistance(pattern, DNA)
		""" print pattern, dist """

		if dist < Distance:
			Distance = dist
			Indices = indices
			winner = pattern

		if dist == 0 :
			print pattern

	return Distance, winner, Indices

# Given a probablity profile for motifs, find the most probable k-mer in a strand
def ProbableKmer (Profile, Strand, k) :
	ks = len(Strand)
	bestIndex = 0
	maxProb = 0
	for i in xrange(ks - k + 1) :
		prob = 1

		for j in xrange (i, i+k) :
			prob *= (Profile[Strand[j]])[j-i]
			if prob == 0 : break

		if prob > maxProb :
			maxProb = prob
			bestIndex = i 

	return Strand[bestIndex:bestIndex+k], maxProb

# Given a probablity profile for all but one motifs, find the most probable k-mer in a strand
def GibbsKmer (Profile, Strand, k) :
	ks = len(Strand)

	prob = [1 for i in range(ks - k + 1)]
	sumProb = 0

	for i in xrange(ks - k + 1) :
		for j in xrange (i, i+k) :
			prob[i] *= (Profile[Strand[j]])[j-i] 
		sumProb += prob[i]

	rand = random.random()*sumProb
	for i in xrange(ks - k + 1) :
		rand -= prob[i]
		if rand <= 0 : break

	return Strand[i:i+k] 


# Build Profile for a set of motifs. Use pseudo counting
def GetProfile(motifs) :
	k = len(motifs[0])

	profile = {'A':[1 for i in range(k)], 'C':[1 for i in range(k)], 
				'G':[1 for i in range(k)], 'T':[1 for i in range(k)]}
	profileSize = 4

	for motif in motifs :
		profile, profileSize = UpdateFrequencyProfile(profile, profileSize, motif)
# Build Profile for a set of motifs. Use pseudo counting
def GetProfile(motifs) :
	k = len(motifs[0])

	profile = {'A':[1 for i in range(k)], 'C':[1 for i in range(k)], 
				'G':[1 for i in range(k)], 'T':[1 for i in range(k)]}
	profileSize = 4

	for motif in motifs :
		profile, profileSize = UpdateFrequencyProfile(profile, profileSize, motif)

	return profile, profileSize


# Update frequency profile with a new motif
def UpdateFrequencyProfile (profile, profileSize, motif) :
	k = len(motif)

	for i in xrange(k) :
		profile[motif[i]][i] += 1
	return profile, profileSize + 1

# Given the motifs and their frequency profile, score the motifs 
def Score(motifs, profile) :

	k = len(motifs[0])
	consensus = ''
	
	for i in xrange(k) :
		maxFreq = 0
		nuc = ''

		for x in ['A', 'C', 'G', 'T'] :
			if profile[x][i] > maxFreq :
				maxFreq = profile[x][i]
				nuc = x
		consensus = consensus + nuc

	score, indices = Pattern2DnaDistance(consensus, motifs)
	return score

				
# Greedy Search for Motifs
def GreedyMotifSearch(dna, k) :
	numStrands = len(dna)
	ks = len(dna[0])

	""" Initialize BestMotifs list with the first k-mers in each DNA strand """
	BestMotifs = []
	profile = {'A':[1 for i in range(k)], 'C':[1 for i in range(k)], 'G':[1 for i in range(k)], 'T':[1 for i in range(k)]}
	profileSize = 4    # Starting size is 4 due to pseudocounting

	for strand in dna :
		BestMotifs.append(strand[0:k])
		profile, profileSize = UpdateFrequencyProfile(profile, profileSize, strand[0:k])

	bestScore = Score(BestMotifs, profile)

	VisitedCache = []
	for i in xrange(ks - k + 1) :
		motifs = []
		profile = {'A':[1 for i in range(k)], 'C':[1 for i in range(k)], 'G':[1 for i in range(k)], 'T':[1 for i in range(k)]}
		profileSize = 4

		if(dna[0][i:i+k] in VisitedCache) : break

		motifs.append(dna[0][i:i+k])    # Add the seed motif
		VisitedCache.append(dna[0][i:i+k])    # Add this to VisitedCache

		profile, profileSize = UpdateFrequencyProfile(profile, profileSize, dna[0][i:i+k])

		for j in xrange(1, numStrands) :
			motif, prob = ProbableKmer(profile, dna[j], k)
			prob = prob/profileSize**k

			motifs.append(motif)
			profile, profileSize = UpdateFrequencyProfile(profile, profileSize, motif)

		score = Score(motifs, profile)
		if score < bestScore :
			bestScore = score
			BestMotifs = motifs

	return BestMotifs, bestScore 


# Randomized search for motif
def RandomMotifSearch(Dna, k) :
	numStrands = len(Dna)
	ks = len(Dna[0])

	globalBestMotifs = ['None']
	globalBestScore = 2**20       # A large positive integer

	random.seed(100)

	for i in range(1000) :  # Run random search 1000 times and take the best results

		bestMotifs = ['None']
		bestScore = 2**20       # A large positive integer

		motifs = []
		for strand in Dna :
			index = int(random.random()*(ks-k+1))
			motifs.append(strand[index:index+k])
		profile, profileSize = GetProfile(motifs)

		while True:	# Loop until no improvement

			score = Score(motifs, profile)
			if score >= bestScore: break	# Break the while loop as there is no improvement
	
			# else

			bestScore = score
			bestMotifs = motifs
			""" print 'Local best score : ', score
			print motifs """

			motifs = []
			for strand in Dna :
				motif, prob = ProbableKmer(profile, strand, k)
				motifs.append(motif)
			
			profile, profileSize = GetProfile(motifs)

		if score < globalBestScore :
			globalBestScore = score
			globalBestMotifs = motifs

	return globalBestMotifs, globalBestScore	


# Motif Search with Gibbs Importance Sampling
def GibbsMotifSearch (Dna, k, numIter) :
        numStrands = len(Dna)
        ks = len(Dna[0])

        random.seed(100)

	globalBestScore = 2**30
	globalBestMotifs = []

	for outerIter in range(20) :

		motifs = []
		for strand in Dna :
			index = int(random.random()*(ks-k+1))
			motifs.append(strand[index:index+k])
		profile, profileSize = GetProfile(motifs)
		score = Score(motifs, profile)

		bestScore = score
		bestMotifs = motifs

        	for iter in range(numIter) :  # Run Gibbs sampling search and take the best results

			# Randomly identify one strand and remove it from the motifs and profile
			i = int(random.random() * numStrands)

			motifs.remove(motifs[i])
			profile, profileSize = GetProfile(motifs)

			motif = GibbsKmer(profile, Dna[i], k)
			motifs.insert(i, motif)
			profile, profileSize = GetProfile(motifs)

			score = Score(motifs, profile)
			if score < bestScore :
				bestMotifs = motifs
				bestScore = score
				print 'Score: ', score, motifs

		if bestScore < globalBestScore :
			globalBestScore = bestScore
			globalBestMotifs = bestMotifs

	return globalBestMotifs, globalBestScore

#-----------------------------------------------------------------------
def MedianString_main(fp) :
	k = int(fp.readline().rstrip())
	Dna = []

	while True :
                x = fp.readline()
                if not x : break
                Dna.append(x.rstrip())
	Distance, Pattern, Indices = MedianString(Dna, k)
	print Pattern, Distance


def ProbableKmer_main(fp) :
	Strand = fp.readline().rstrip()
	k = int(fp.readline().rstrip())
        Profile = {}

        for x in ['A', 'C', 'G', 'T']:
                prof = fp.readline().rstrip().split()
                Profile[x] = [float(p) for p in prof]

	kmer, prob = ProbableKmer(Profile, Strand, k)
	print kmer, prob

def GreedyMotifSearch_main(fp) :
	y = fp.readline().rstrip().split()
	k, numStrands, numIter = [int(x) for x in y]
	Dna = []
	
	for i in xrange(numStrands) :
		Dna.append(fp.readline().rstrip())

	BestMotifs, score = GreedyMotifSearch(Dna, k)

	for motif in BestMotifs :
		print motif

	print
	print 'Score : ', score


def RandomMotifSearch_main(fp) :
	y = fp.readline().rstrip().split()
	k, numStrands, numIter = [int(x) for x in y]
	Dna = []
	
	for i in xrange(numStrands) :
		Dna.append(fp.readline().rstrip())

	BestMotifs, score = RandomMotifSearch(Dna, k)

	for motif in BestMotifs :
		print motif


	print
	print 'Score : ', score

def GibbsMotifSearch_main(fp) :
        y = fp.readline().rstrip().split()
        k, numStrands, numIter = [int(x) for x in y]
        Dna = []

        for i in xrange(numStrands) :
                Dna.append(fp.readline().rstrip())

        BestMotifs, score = GibbsMotifSearch(Dna, k, numIter)

        for motif in BestMotifs :
                print motif


        print
        print 'Score : ', score

#-----------------------------------------------------------------------

fn = raw_input('Input file: ')
fp = open(fn, 'r')
MedianString_main(fp)
# ProbableKmer_main(fp)
# GreedyMotifSearch_main(fp)
# RandomMotifSearch_main(fp)
# GibbsMotifSearch_main(fp)
