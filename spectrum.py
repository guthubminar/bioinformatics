#!/usr/bin/env python

AminoAcidMass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

def Spectrum(Peptide, isCyclic) :
	lPeptide = len(Peptide)
	PrefixMass = range(lPeptide+1)
	PrefixMass[0] = 0
	
	for i in xrange(1, lPeptide+1) :
		PrefixMass[i] = PrefixMass[i-1] + AminoAcidMass[Peptide[i-1]]

	peptideMass = PrefixMass[lPeptide]

	Spectrum = []
	Spectrum.append(0)

	for i in xrange(lPeptide) :
		for j in xrange(i+1, lPeptide+1) :
			Spectrum.append(PrefixMass[j] - PrefixMass[i])

            		if ((isCyclic == 1) and (i > 0) and (j < lPeptide)) :
				Spectrum.append(peptideMass - (PrefixMass[j] - PrefixMass[i]))

	return sorted(Spectrum)

peptide = 'KQNWMDMHAVLTG'
l = Spectrum(peptide, 0)
for i in xrange(len(l)) :
	print l[i],

print 
print 

c = Spectrum(peptide, 1)
for i in xrange(len(c)) :
	print c[i],
