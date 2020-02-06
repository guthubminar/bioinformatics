#!/usr/bin/env python

from graph import *
from CycloPeptide import *
import sys

AAMass = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

AA2Mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q':128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
# AA2Mass = {'X':4, 'Z':5}

Mass2AA = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'I', 114: 'N', 115: 'D', 128: 'K', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}


"""
Construct a directed graph from a peptide's spectrum
	o Spectrum is a vector of masses of prefixes and suffixes of a peptide
	o The graph will have one node for each entry in the spectrum, plus a node for 0 mass (empty prefix). The mass is stored as node's data.
	o It will have an edge between every pair of nodes whose mass difference equals the mass of an aminoacid.
"""

def Spectrum2Graph(spectrum) :
	numNodes = len(spectrum) + 1

	G = Graph(DIRECTED)

	node = G.addNode()
	node.data = 0

	for i in range(1, numNodes) :
		node = G.addNode()
		node.data = spectrum[i-1]

	for i in range(numNodes) :
		lNode = G.Nodes[i]
		for j in range(i+1, numNodes) :
			rNode = G.Nodes[j]
			mass = rNode.data - lNode.data

			if mass in Mass2AA :
				edge = G.addEdge(lNode, rNode)
				edge.data = Mass2AA[mass]

	return G


"""
Construct a directed graph from a spectral vector
	o Suppose [s1, s2, ...., si, ......sm] is the Spectral vector. 
	o In the vector: position i is the mass of a prefix peptide; value si is the likelihood that the prefix of a (unknown) peptide will annotate the intensity of spectrum at mass i.
	o The graph will have one node for each entry in the spectrum, plus a node for 0 mass (empty prefix). The node id (i) denotes the mass. Weight si is stored as data of node i.
	o It will have an edge between every pair of nodes whose mass (id) difference equals the mass of an aminoacid.
"""

def SpectralVector2Graph(sVector) :
	numNodes = len(sVector) + 1

	G = Graph(DIRECTED)

	node = G.addNode()
	node.data = 0

	for i in range(1, numNodes) :
		node = G.addNode()
		node.data = sVector[i-1]

	for i in range(numNodes) :
		lNode = G.Nodes[i]
		for j in range(i+1, numNodes) :
			rNode = G.Nodes[j]
			mass = j - i

			if mass in Mass2AA :
				edge = G.addEdge(lNode, rNode)
				edge.data = Mass2AA[mass]

	return G



# Generate a spectrum of masses of prefixes or suffixes or both of a peptide
FULL_SPECTRUM = 0
PREFIXES_ONLY = 1
SUFFIXES_ONLY = 2
def Peptide2Spectrum(peptide, scope = FULL_SPECTRUM) :
	lPeptide = len(peptide)

	if scope == FULL_SPECTRUM or scope == PREFIXES_ONLY :
		prefixSpectrum = range(lPeptide-1)
		prefixSpectrum[0] = AA2Mass[peptide[0]]
		for i in range(1, lPeptide-1) :
			prefixSpectrum[i] = prefixSpectrum[i-1] + AA2Mass[peptide[i]]

	if scope == FULL_SPECTRUM or scope == SUFFIXES_ONLY :
		suffixSpectrum = range(lPeptide)
		suffixSpectrum[0] = AA2Mass[peptide[lPeptide-1]]
		for i in range(1, lPeptide) :
			suffixSpectrum[i] = suffixSpectrum[i-1] + AA2Mass[peptide[lPeptide-1-i]]

	if scope == PREFIXES_ONLY :
		prefixSpectrum.append(prefixSpectrum[-1] + AA2Mass[peptide[lPeptide-1]])
		return sorted(prefixSpectrum)

	if scope == SUFFIXES_ONLY :
		suffixSpectrum.append(0)
		return sorted(suffixSpectrum)

	# else return FULL_SPECTRUM
	return sorted(prefixSpectrum + suffixSpectrum)

# Generate a PeptideVector of a given Peptide
def Peptide2Vector(peptide) :
	spectrum = Peptide2Spectrum(peptide, PREFIXES_ONLY)
	# print sorted(spectrum)
	vector = [0 for i in range(spectrum[-1])]
	for i in spectrum :
		vector[i-1] = 1 
	return vector

# Find the peptide from a peptide vector
def Vector2Peptide(vector) :
	peptide = ''
	mass = 0
	for i in vector :
		mass += 1
		if i == 1: 
			peptide += Mass2AA[mass]
			mass = 0
	return peptide

"""
Find the peptide represented by the given spectrum
Algorithm:
	o Construct a directed graph, G, from the spectrum
	o Iterate over all paths from the source (node with 0 mass) to the sink (node with full mass) of G :
		o Get the peptide given by the edges of the path
		o Find the spectrum of this peptide
		o If this spectrum is identical to the given spectrum, we have the peptide we are looking for.
"""
def DecodeSpectrum(spectrum) :
	G = Spectrum2Graph(spectrum)
	numNodes = len(G.Nodes)
	source = G.Nodes[0]
	sink = G.Nodes[-1]
	for path in DFS_paths(G, source, sink) :
		peptide = ''
		for i in range(len(path)-1) :
			lNode = G.Nodes[path[i].id]
			rNode = G.Nodes[path[i+1].id]
			edge = G.getEdge(lNode, rNode)
			peptide += edge.data
		# print peptide

		idealSpectrum = Peptide2Spectrum(peptide, FULL_SPECTRUM)
		# print idealSpectrum
		if sorted(spectrum) == sorted(idealSpectrum) :
			return peptide
	
"""
Find the peptide that 'best' represents a given spectral vector.
The best peptide is given by the max (node) weighted path in the graph of the spectral vector.

Algorithm:
	o Input: Spectral Vector, sVector = [s1, s2, ....., sm]
	o Create a directed acyclic graph, G, from the spectral vector, sVector. 
		G has nodes with id 0 through m, and node[i].data = si
  	o Initialize an array, 'paths' = [ (d0, p0), (d1, p1), ......(dm, pm) ], 
		where is the prefix peptide at x, and dx is the total weight of nodes from source to node x.
		Initial values are dx = -infinity (except d0 = 0), and px = ''.
	o Starting from source, update (dx, px) in topological order.
		Note, the order of nodes, viz. 0, 1, 2, ...., m is already a topo-sort.
	o Output: pm 
"""
def DecodeSpectralVector(sVector) :
	m = len(sVector)
	G = SpectralVector2Graph(sVector)
	source = G.Nodes[0]
	sink = G.Nodes[m]

	paths = [ (-sys.maxint, '') for i in range(m+1) ]
	paths[0] = (0, '')

	for i in range(m) :
		di, pi = paths[i]
		print
		print i, ':', di, pi
		print '-------------'
		for e in G.Nodes[i].outEdges :
			j = e.rNode.id
			dj, pj = paths[j]
			d = di + e.rNode.data
			if d > dj : 
				paths[j] = (d, pi + e.data)
				print j, ':', paths[j]

	score, peptide = paths[m]
	print score, peptide
	return peptide
#---------------------------------------------------------------------------------------
# Given a spectral vector, sVector, identify the peptide (subsequence of a given proteome)
#	that scores maximum with the sVector and has mass equal to the length of the spectrum

def IdentifyPeptide(sVector, proteome) :
	sVector = [0] + sVector

	m = len(sVector)-1; n = len(proteome)
	print proteome; print sVector
	print 'Proteome Length:', n, 'Full mass:', m

	maxScore = -sys.maxint
	peptide = ''
	for i in range(n) :
		mass = AA2Mass[proteome[i]]
		if mass > m : continue
		score = sVector[mass]

		print
		for j in range(i+1, n) :
			delta = AA2Mass[proteome[j]]
			if mass + delta > m : break

			mass += delta
			score += sVector[mass]
			print i, ':', proteome[i:j+1], 'Mass:', mass, 'Score:', score,
			if score >= maxScore and mass == m :
				maxScore = score
				peptide = proteome[i:j+1]
				print '<- Best'
			else: print

	print 'Best :', peptide, maxScore
	return peptide, maxScore
#---------------------------------------------------------------------------------------
# Find the Peptide-Spectrum Match (PSM) set
def PSMSearch(sVectors, proteome, threshold) :
	PSMSet = []
	for sVector in sVectors :
		peptide, score = IdentifyPeptide(sVector, proteome)
		print score, peptide
		if score >= threshold : PSMSet.append(peptide)

	return PSMSet
#---------------------------------------------------------------------------------------
# Calculate the Size (or Probability) of a Spectral Dictionary
SIZE = 1
PROBABILITY = 0
def SpectralDictionarySize(sVector, minThresh, maxThresh, sizeOrProb = SIZE) :
	sVector = [0] + sVector
	fullMass = len(sVector) - 1
	cache = {}

	def Size(mass, score, depth) :	
		# for x in range(depth) : print ' ',
		# print 'size(', mass, score, ')'

		if mass == 0 and score == 0 : return 1
		if mass <= 0 or score < 0 : return 0
		if (mass, score) in cache : return cache[(mass, score)]

		size = 0
		for aa in AA2Mass :
			# print aa,
			size += Size(mass - AA2Mass[aa], score - sVector[mass], depth+1) 

		if sizeOrProb == PROBABILITY : size *= 0.05
		cache[(mass, score)] = size
		# print 'Size(', mass, score, ') = ', size
		return size

	count = 0
	for t in range(minThresh, maxThresh+1) :
		count += Size(fullMass, t, 0)

	# print len(cache)
	return count
#---------------------------------------------------------------------------------------
"""
Spectral Alignment:
Input: A peptide, a spectral vector, and maxChanges
Output: A mutated peptide that maximizes the score w.r.to the spectral vector.
Constraints: 
	o Number of aminoacids in the peptide whose mass changes <= maxChanges
	o Total mass of the mutated peptide = Length of the spectral vector
How this is done: 
	o Using dynamica programming.
	o Propagate solutions in topo order from source to sink nodes
	o Drop on the way all infeasible and inferior solutions
Details:
	o Solutions is an array of size equal to the size of spectral vector.
	o At each node, there are multiple solutions, stored in dictionaries.
	o A solution is a dictionary entry: (size, changes) : (score, path)
		size = size of peptide built so far
		changes = number of mutations so far
		score, path = score and path of a solution
"""
def SpectralAlignment(peptide, sVector, maxChanges) :
	#----------------------------------------------------------------------
	def UpdateSolution(index, size, changes, score, path) :
		if changes > maxChanges or size > S : return
		if index == N and size != S : return

		sol = Solutions[index]
		for (s, c) in sorted(sol) :
			if (s == size) :
				if (c <= changes) and (sol[(s,c)][0] >= score) : return
				if (c > changes) and (sol[(s,c)][0] <= score) : del sol[(s,c)]
		sol[(size,changes)] = (score, path)
	#----------------------------------------------------------------------
		
	N = len(sVector) 
	sVector = [0] + sVector
	S = len(peptide)

	Solutions = [{} for i in range(N+1)]
	Solutions[0] = {(0,0) : (0, '')}

	for i in range(N) :
		for (size, changes) in Solutions[i] :
			for j in range(i+1, N+1) :
				if size == S-1 : 
					if j != N : continue
					if (changes == maxChanges) and (N-i != AA2Mass[peptide[size]]) : continue
				elif j == N : continue

				s = size + 1
				score, path = Solutions[i][(size, changes)]
				score += sVector[j]
				path += peptide[size]
				delta = j - i - AA2Mass[peptide[size]]

				if delta == 0 :
					c = changes
				else : 
					c = changes + 1
					if delta > 0 : path += '(+' + str(delta) + ')'
					else : path += '(' + str(delta) + ')'
				UpdateSolution(j, s, c, score, path)
	
	bestScore = -sys.maxint
	sol = Solutions[N]
	for key in sol :
		score, path = sol[key]
		if score > bestScore :
			bestScore = score
			alignment = path
			
	return alignment
#---------------------------------------------------------------------------------------
# Graph construction from a peptide's spectrum
def p11_1(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	spectrum = fp.readline().strip().split()
	spectrum = [int(x) for x in spectrum]

	G = Spectrum2Graph(spectrum)

	sys.stdout = fpOut
	for id in G.Edges :
		print str(G.Edges[id].lNode.data) + '->' + str(G.Edges[id].rNode.data) + ':' +  G.Edges[id].data

#---------------------------------------------------------------------------------------
# Decoding a given spectrum to its peptide
def p11_2(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	spectrum = fp.readline().strip().split()
	spectrum = [int(x) for x in spectrum]

	peptide = DecodeSpectrum(spectrum)
	sys.stdout = fpOut
	print peptide
#---------------------------------------------------------------------------------------
# Constructing a peptide vector
def p11_3(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	peptide = fp.readline().strip()
	vector = Peptide2Vector(peptide)
	sys.stdout = fpOut
	for i in vector :
		print i, 
#---------------------------------------------------------------------------------------
# Constructing a peptide from a peptide vector
def p11_4(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	vector = fp.readline().strip().split()
	vector = [int(i) for i in vector]
	peptide = Vector2Peptide(vector)
	sys.stdout = fpOut
	print peptide 
#---------------------------------------------------------------------------------------
# Decode a spectral vector to a peptide that scores best against the spectral vector
def p11_5(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	sVector = fp.readline().strip().split()
	sVector = [int(i) for i in sVector]
	peptide = DecodeSpectralVector(sVector)
	sys.stdout = fpOut
	print peptide 
#---------------------------------------------------------------------------------------
# Identify a peptide (subsequence of the given proteome) that scores best with a given spectral vector 
def p12_1(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	sVector = fp.readline().strip().split()
	sVector = [int(i) for i in sVector]
	proteome = fp.readline().strip()
	peptide, score = IdentifyPeptide(sVector, proteome)
	sys.stdout = fpOut
	print peptide 
##---------------------------------------------------------------------------------------
# Peptide-Spectrum Matching
def p12_2(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	sVectors = []
	while True:
		data = fp.readline().strip().split()
		if len(data) > 1 :
			sVector = [int(i) for i in data]
			sVectors.append(sVector)
		else :
			proteome = data[0]
			break
	threshold = int(fp.readline().strip())

	PSMSet = PSMSearch(sVectors, proteome, threshold)
	sys.stdout = fpOut
	for peptide in PSMSet : print peptide 
#---------------------------------------------------------------------------------------
def p12_3(infile, outfile, sizeOrProbability = SIZE) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	sVector = fp.readline().strip().split()
	sVector = [int(x) for x in sVector]
	minThresh = int(fp.readline().strip())
	maxThresh = int(fp.readline().strip())

	size = SpectralDictionarySize(sVector, minThresh, maxThresh, sizeOrProbability)
	print size
	sys.stdout = fpOut
	print size
#---------------------------------------------------------------------------------------
def p12_5(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	peptide = fp.readline().strip()
	sVector = fp.readline().strip().split()
	sVector = [int(x) for x in sVector]
	maxChanges = int(fp.readline().strip())

	print len(sVector)
	mass = 0 
	for aa in peptide : mass += AA2Mass[aa]
	print 'Peptide Mass:', mass

	alignment = SpectralAlignment(peptide, sVector, maxChanges)
	print alignment
	sys.stdout = fpOut
	print alignment
#---------------------------------------------------------------------------------------
# p11_1('./Regression/11.1.1', './Solutions/11.1.1')
# p11_1('./Regression/11.1.2', './Solutions/11.1.2')
# p11_1('./Regression/11.1.3', './Solutions/11.1.3')
# p11_2('./Regression/11.2.1', './Solutions/11.2.1')
# p11_2('./Regression/11.2.2', './Solutions/11.2.2')
# p11_2('./Regression/11.2.3', './Solutions/11.2.3')
# p11_3('./Regression/11.3.1', './Solutions/11.3.1')
# p11_3('./Regression/11.3.2', './Solutions/11.3.2')
# p11_3('./Regression/11.3.3', './Solutions/11.3.3')
# p11_4('./Regression/11.4.1', './Solutions/11.4.1')
# p11_4('./Regression/11.4.2', './Solutions/11.4.2')
# p11_5('./Regression/11.5.1', './Solutions/11.5.1')
# p11_5('./Regression/11.5.2', './Solutions/11.5.2')
# p11_5('./Regression/11.5.3', './Solutions/11.5.3')
# p12_1('./Regression/12.1.1', './Solutions/12.1.1')
# p12_1('./Regression/12.1.2', './Solutions/12.1.2')
# p12_1('./Regression/12.1.3', './Solutions/12.1.3')
# p12_2('./Regression/12.2.1', './Solutions/12.2.1')
# p12_2('./Regression/12.2.2', './Solutions/12.2.2')
# p12_2('./Regression/12.2.3', './Solutions/12.2.3')
# p12_3('./Regression/12.3.1', './Solutions/12.3.1')
# p12_3('./Regression/12.3.2', './Solutions/12.3.2')
# p12_3('./Regression/12.3.3', './Solutions/12.3.3')
# p12_3('./Regression/12.4.1', './Solutions/12.4.1', PROBABILITY)
# p12_3('./Regression/12.4.2', './Solutions/12.4.2', PROBABILITY)
# p12_3('./Regression/12.4.3', './Solutions/12.4.3', PROBABILITY)
# p12_5('./Regression/12.5.1', './Solutions/12.5.1')
# p12_5('./Regression/12.5.2', './Solutions/12.5.2')
# p12_5('./Regression/12.5.3', './Solutions/12.5.3')
p12_5('./Regression/12.5.4', './Solutions/12.5.4')
