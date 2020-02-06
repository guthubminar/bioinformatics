# Utilities for Bioinformatic Processing
# Author: Rajendran
# Created: 11/13/2014

# Find Hamming Distance between 2 equal length strings
def HammingDist(pattern1, pattern2):
	l = len(pattern1)

	dist = 0	
	for i in xrange(0,l):
		if(pattern1[i] != pattern2[i]):
			dist += 1

	return dist

# Count the number of times 'pattern' appears in 'text' with Hamming Distance <= distance
def CountPattern(text, pattern, distance=0):
	lText = len(text)
	lPattern = len(pattern)

	count = 0 
	for i in range (0, lText - lPattern + 1) :
        	p = text[i:i+lPattern] 
        	d = HammingDist(p, pattern) 
        	if (d <= distance) :
			print i,
                	count = count + 1
	return count

# Print all k-mers occuring within 'text' with Hamming Distance not more than d
def FrequentWords(text, k, d):
	n = {0:'A', 1:'G', 2:'T', 3:'C'}

	lText = len(text)
	max = 0
	count = {}

	for i in xrange(4**k):
		x = i
		pattern = ''	
		for j in xrange(k):
			r = x % 4
			x = x / 4
			pattern = n[r] + pattern 
		
		count[i] = CountPattern(text, pattern, d)
		if (count[i] > max):
			max = count[i]

	for i in xrange (4**k) :
		if (count[i] == max) :
			x = i
			pattern = ''	
			for j in xrange(k):
				r = x % 4
				x = x / 4
				pattern = n[r] + pattern 
			print pattern,

# Print the positions of minimum skew in the DNA string 'text' 
def minSkew(text):
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


# Print Reverse Complement of DNA string 'text'
def ReverseComplement(text):
	lText = len(text)
	comp = {}
	comp['A'] = 'T'
	comp['T'] = 'A'
	comp['G'] = 'C'
	comp['C'] = 'G'

	revComp = []

	for i in xrange (0, lText):
		revComp.append(comp[text[i]])

	revComp.reverse()
	return ''.join(revComp)


TCAG = 'TCAG'
UCAG = 'UCAG'

tcag2int = {'T':0, 'C':1, 'A':2, 'G':3}
ucag2int = {'U':0, 'C':1, 'A':2, 'G':3}

# Amino Acid
# Below, o, a and p are stop codons (ocher, amber, and opal respectively); other letters are amino acids
AA = 'FFLLSSSSYYoaCCpWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

# int2nuc = ['T', 'C', 'A', G']
# nuc2int = {'T':0, 'C':1, 'A':2, 'G':3}

# aa = ['F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', 'Ocher', 'Amber', 'C', 'C', 'Opal', 'W',
#      'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R',
#      'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R',
#      'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G']

# Convert n (0 to 63) to codon
def int2codon (n):
	codon = ''
	for i in xrange(3):
		r = n % 4
		n = n / 4
		codon = UCAG[r] + codon 
	return codon

# convert a codon to integer (0 to 63)
codon2int = lambda codon: 16*ucag2int[codon[0]] + 4*ucag2int[codon[1]] + ucag2int[codon[2]]


# Look up AminoAcid (AA) for a codon
codon2aa = lambda codon: AA[codon2int(codon)]

# Transcribe a 'dna' string to 'rna'
def transcribe(dna):
	trans = {}
	trans['A'] = 'A'
	trans['T'] = 'U'
	trans['G'] = 'G'
	trans['C'] = 'C'

	l = len(dna)
	rna = []
	for i in xrange (0, l):
                rna.append(trans[dna[i]])

	return ''.join(rna)


# Translate an 'rna' string to a 'peptide' chain 
def translate(rna) :
	lRna = len(rna)
	peptide = ''
	for i in xrange(0, lRna, 3):
		codon = rna[i:i+3]
		aa = codon2aa(codon)
		if ( (aa == 'o') |  (aa == 'a') |  (aa == 'p')) :
			break
		else:
			peptide += aa
	return peptide
