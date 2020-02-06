#!/usr/bin/env python

AAMass = [57 71 87 97 99 101 103 113 114 115 128 129 131 137 147 156 163 186]

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
