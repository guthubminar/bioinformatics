#!/usr/bin/env python

from itertools import ifilterfalse
from collections import Counter
from collections import namedtuple

AminoAcidMass = {'G': 57, 'Z': 65, 'A': 71, 'S': 87, 'P': 97, 'X': 98, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q':128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

AAMass = [57, 65, 71, 87, 97, 98, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

Peptide = namedtuple("p", "peptide mass score")
	

# Test if list1 is a subset of list2 
# When the lists have repeated elements, the counts of elements need to be considered

def counterSubset(list1, list2):
	c1, c2 = Counter(list1), Counter(list2)
	for k, n in c1.items():
		if n > c2[k]:
			return False
        return True



# Return the mass of a peptide

def Mass(peptide) :
	l = len(peptide)
	m = 0
	for i in xrange(l) :
		m += i
	return m



# Inline function to return parent (max) mass of a mass spectrum
ParentMass = lambda spectrum : spectrum[len(spectrum)-1]



# Return linear/cyclic spectrum of Peptide if isCyclic is 0/1

def Spectrum(Peptide, isCyclic) :
	lPeptide = len(Peptide)
	PrefixMass = range(lPeptide+1)
	PrefixMass[0] = 0
	
	for i in xrange(1, lPeptide+1) :
		PrefixMass[i] = PrefixMass[i-1] + Peptide[i-1]

	peptideMass = PrefixMass[lPeptide]

	Spectrum = []
	Spectrum.append(0)

	for i in xrange(lPeptide) :
		for j in xrange(i+1, lPeptide+1) :
			Spectrum.append(PrefixMass[j] - PrefixMass[i])

            		if ((isCyclic == 1) and (i > 0) and (j < lPeptide)) :
				Spectrum.append(peptideMass - (PrefixMass[j] - PrefixMass[i]))

	return sorted(Spectrum)




# Check if a peptide's linear/cyclic spectrum is consistent with the given spectrum

def Consistent(peptide, spectrum, isCyclic)  :
	s = Spectrum(peptide, isCyclic)
	""" print s, spectrum """
	return counterSubset(s, spectrum)


# Return the scoring for a peptide's linear/cyclic mass spectrum w.r.to the given experimental spectrum

def PeptideScoring(peptide, spectrum, isCyclic)  :
	s = Spectrum(peptide, isCyclic)
	c1, c2 = Counter(s), Counter(spectrum)

	score = 0
	for k, v in c1.items():
		if not k in c2.keys():
			continue
		elif v <= c2[k] :
			score += v
		else :
			score += c2[k]
	return score


# Return all peptides whose cyclic spectrum is consistent with given Spectrum

def CycloPeptideSequencing(Spectrum) :
        parentMass = ParentMass(Spectrum)
        """ print parentMass """

        rightMassList = []
        goodList = []
        Peptides = []
        Peptides.append('')
        numPep = 1

        while((numPep) != 0) :
                """" Expand every peptide in Peptides in both directions by one amino acid """
                newlist = []
                for peptide in Peptides :
                        """ print 'Expanding: ', peptide """
                        for aa in AminoAcidMass.keys() :
                                pep = peptide + aa
                                if (not pep in newlist) :
                                        newlist.append(pep)

                                pep = aa + peptide
                                if (not pep in newlist) :
                                        newlist.append(pep)

                Peptides = newlist
                """ print 'After expansion: ', Peptides """

		""" Remove peptides whose linear spectrum is inconsistent with the given Spectrum """
		Peptides[:] = [peptide for peptide in Peptides if not Consistent(peptide, Spectrum, 0)]
		""" print 'Consitent ones: ', Peptides """

                """ Find the peptides which have mass = parentMass and given spectrum """
                rightMassList[:] = [peptide for peptide in Peptides if Mass(peptide) == parentMass]
                """ print 'Right mass: ', rightMassList """

                for peptide in rightMassList:
                        if (Consistent(peptide, Spectrum, 1)) :
                                goodList.append(peptide)

                """ Remove peptides matching parentMass from the list before next iteration """
                Peptides[:] = [peptide for peptide in Peptides if not Mass(peptide) == parentMass]
                """ print Peptides """

                numPep = len(Peptides)
        return goodList



# Return peptides with atmost N top scores

def CycloPeptideApproxSequencing(Spectrum, N) :
	parentMass = ParentMass(Spectrum)
	""" print parentMass """

	goodList = []
	rightMassList = []
	candidates = []

	emptyPeptide = Peptide('', 0, 0)
	candidates.append(emptyPeptide)
	best = emptyPeptide
	maxScore = 0
	numPep = 1

	while((numPep) != 0) :
		"""" Expand every peptide in Peptides in both directions by one amino acid """
		newlist = []
		for item in candidates :
			""" print 'Expanding: ', peptide """
			for aa in AminoAcidMass.keys() :
				pep = item.peptide + aa
				mass = item.mass + AminoAcidMass[aa]
				score = PeptideScoring(pep, Spectrum, 0)
	
				p = Peptide(pep, mass, score)
				if (not p in newlist) :
					newlist.append(p)

				pep = aa + item.peptide 
				mass = item.mass + AminoAcidMass[aa]
				score = PeptideScoring(pep, Spectrum, 0)

				p = Peptide(pep, mass, score)
				if (not p in newlist) :
					newlist.append(p)

		candidates = newlist
		""" print 'After expansion: ', candidates """

		""" Calculate Cyclic Peptide Score for peptides equalling parentMass """
		for x in candidates:
			if (x.mass == parentMass) :
				s = PeptideScoring(x.peptide, Spectrum, 1)
				x = Peptide(x.peptide, x.mass, s)
				if x.score > maxScore:
					maxScore = x.score
					best = x
					goodList.append(x)
				
		
		""" Keep only peptides of mass less than parentMass """
		candidates[:] = [x for x in candidates if x.mass < parentMass]
		""" print 'After removing heavier ones: ', candidates """


		""" Trim the candidates list to keep only the top N scorers """
		candidates.sort(key = lambda tup: tup[2], reverse=True)
		""" print 'Sorted : ', candidates """

		numPep = len(candidates)

		if numPep > N :
			numPep = N

			"""
			minScore = candidates[N-1].score

			for x in candidates[N:] :
				if x.score == minScore:
					numPep += 1
				else :
					break

			"""

			candidates = candidates[:numPep]	
			""" print 'Trimmed : ', candidates """

	goodList.sort(key = lambda tup: tup[2], reverse=True)
	return goodList

def Convolution(spectrum, min, max) :
	convolution = []
	n = len(spectrum)
	for i in xrange (n) :
		for j in xrange (i) :
			x = spectrum[i] - spectrum[j]
			if ( (x >= min) & (x <= max) ):
				convolution.append(x)
	return sorted(convolution)

def CPSequence(spectrum, m, n) :
	spectrum[:] = sorted(spectrum)
	convolution = Convolution(spectrum, 57, 200)

	convolution = Counter(convolution)
	
	counts = sorted(convolution.values(), reverse=True)
	minFrequency = counts[m-1]

	masses = []
	for k in convolution.keys() :
		if convolution[k] >= minFrequency :
			masses.append(k)

	parentMass = ParentMass(spectrum)

        goodList = []
        candidates = []

	for mass in masses :
		peptide = [mass]
		score = PeptideScoring(peptide, spectrum, 0)
		pep = Peptide(peptide, mass, score)
		if (score > maxScore) :
			maxScore = score
			best = pep
		candidates.append(pep)

        numPep = len(candidates)

        while((numPep) != 0) :
                """" Expand every peptide in candidates in both directions by one amino acid """
                newlist = []
                for pep in candidates :
                        for mass in masses :
				pepmass = pep.mass + mass
				if (pepmass > parentMass):
					break

				p1 = []
				p2 = []
				p1[:] = pep.peptide
				p2[:] = pep.peptide
                                p1.insert(0, mass)
                                p2.append(mass)

				if (pepmass == parentMass) :
					score1 = PeptideScoring(p1, spectrum, 1)
					score2 = PeptideScoring(p2, spectrum, 1)
				else :
					score1 = PeptideScoring(p1, spectrum, 0)
					score2 = PeptideScoring(p2, spectrum, 0)

				pep1 = Peptide(p1, pepmass, score1)
				pep2 = Peptide(p2, pepmass, score2)

				if (pepmass == parentMass) :
					if (not pep1 in goodList) :
						goodList.append(pep1)

					if (not pep2 in goodList) :
						goodList.append(pep2)
				else :
					if (not pep1 in newlist) :
						newlist.append(pep1)

					if (not pep2 in newlist) :
						newlist.append(pep2)

                candidates = newlist

                """ Trim the candidates list to keep only the top 'n' scorers """
                candidates.sort(key = lambda tup: tup[2], reverse=True)

                numPep = len(candidates)

                if numPep > n :
                        numPep = n
                        minScore = candidates[n-1].score

                        for x in candidates[n:] :
                                if x.score == minScore:
                                        numPep += 1
                                else :
                                        break

                        candidates = candidates[:numPep]

        goodList.sort(key = lambda tup: tup[2], reverse=True)
        return goodList


# -----------------------------------------------
# Below are main functions for class assignments
# -----------------------------------------------

def CPSequencing_main(fp) :
	inlist = fp.readline().rstrip().split(" ")
	spectrum = [int(x) for x in inlist]

	goodList = CycloPeptideSequencing(spectrum)

	for peptide in goodList:
		s = ''
		l = len(peptide)
		s = str(AminoAcidMass[peptide[0]])
		for i in xrange(1, l):
			s = s + '-' + str(AminoAcidMass[peptide[i]])
	print s,	

def CPScoring_main(fp) :
	peptide = fp.readline().rstrip()
	inlist = fp.readline().rstrip().split(" ")
	spectrum = [int(x) for x in inlist]
	print PeptideScoring(peptide, spectrum, 1)

def CPApproxSequencing_main(fp) :
	N = int(fp.readline().rstrip())
	print N

	inlist = fp.readline().rstrip().split(" ")
	print inlist
	spectrum = [int(x) for x in inlist]
	print spectrum

	goodList = CycloPeptideApproxSequencing(spectrum, N)

	maxScore = goodList[1].score

	for x in goodList :
		if x.score < maxScore:
			break
		s = ''
		l = len(x.peptide)
		s = str(AminoAcidMass[x.peptide[0]])
		for i in xrange(1, l):
			s = s + '-' + str(AminoAcidMass[x.peptide[i]])
		print s,	

def Convolution_main(fp) :
	inlist = fp.readline().rstrip().split(" ")
	print inlist
	spectrum = [int(x) for x in inlist]
	print spectrum
	
	ans = Convolution(sorted(spectrum))
	for i in ans:
		print i,

def CPSequence_main(fp) :
	m = int(fp.readline().rstrip())
	n = int(fp.readline().rstrip())
	inlist = fp.readline().rstrip().split(" ")
	spectrum = [int(x) for x in inlist]
	print spectrum

	goodList = CPSequence(spectrum, m, n)
	maxScore = goodList[0].score

	for pep in goodList :
		if pep.score < maxScore :
			break

		l = len(pep.peptide)
		s = str(pep.peptide[0])
		for i in xrange(1, l) :
			s = s + '-' + str(pep.peptide[i])
		print s,

#---------------------------------------------------
fn = raw_input('file? ')
fp = open(fn, 'r')

# CPSequencing_main(fp)
# CPScoring_main(fp)
# CPApproxSequencing_main(fp)
# Convolution_main(fp)
CPSequence_main(fp)
