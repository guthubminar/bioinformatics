#!/usr/bin/env python

from graph import *
from BioInfoUtils import *
from copy import *
import sys

def wait() :
	raw_input("Press enter to continue")

def GenerateDistanceMatrix(G, numLeaves) :
	leaves = []
	leaves = [G.Nodes[id] for id in xrange(numLeaves)]

	matrix = {}
	for i in xrange(numLeaves) : matrix[(i,i)] = 0

	for start in leaves :
		G.clearData(NODES)
		start.data = 0
		print 'From Leaf:', start.id
		dfs_visit = DFS(G, start.id)
		dfs_visit.next()  # Ignore the first node visited, which is the start node.
		for n, e in dfs_visit :
			print 'Visiting: ', n.id
			if e.lNode == n : n.data = e.rNode.data + e.data
			else : n.data = e.lNode.data + e.data
			if n in leaves :
				matrix[(start.id, n.id)] = n.data
	return matrix

# Find Limb Length for leaf L from the Distance Matrix, D of size N x N
def LimbLength(D, N, L) :
	limbLength = sys.maxint
	for i in xrange(N) :
		if i == L : continue
		for j in xrange(N) :
			if j == L or j == i : continue
			length = (D[L][i] + D[L][j] - D[i][j])/2
			if length < limbLength : limbLength = length
	return limbLength

# Construct an Additive Phylogeny Tree from a Distance Matrix, D, of size n x n
def AdditivePhylogeny(D, G, n) :
	print 'At n : ', n

	if n == 2 :
		edge = G.addEdge(G.Nodes[0], G.Nodes[1])
		edge.data = D[0][1]
		edge = G.addEdge(G.Nodes[1], G.Nodes[0])
		edge.data = D[1][0]
		return G

	limbLength = LimbLength(D, n, n-1) 
	print 'Limb: ', limbLength

	for j in xrange(n-1) :
		D[j][n-1] -= limbLength
		D[n-1][j] = D[j][n-1]

	found = 0
	for i in xrange(n-1) :
		for k in xrange(n-1) :
			if k == i : continue
			if D[i][k] == D[i][n-1] + D[n-1][k] :
				found = 1
				break
		if found : break
	if not found : 
		print 'Error: AdditivePhylogeny failed.'
		return

	G = AdditivePhylogeny(D, G, n-1)

	print 'Back at n: ', n
	print 'Graph:'
	G.print_adjacency(1)

	paths = list(DFS_paths(G, G.Nodes[i], G.Nodes[k]))
	path = paths[0]	# Assuming there is a unique path

	nodesInPath = path
	path = []
	for x in xrange(len(nodesInPath)-1) :
		edge = G.getEdge(nodesInPath[x], nodesInPath[x+1])
		path.append(edge)

	print  'Path:'
	for e in path : print e.id, e.data
	print 

	distance = D[i][n-1]
	for edge in path :
		if edge.data < distance :
			distance -= edge.data
			continue
		elif edge.data == distance :
			newEdge = G.addEdge(edge.rNode, G.Nodes[n-1])
			newEdge.data = limbLength
			newEdge = G.addEdge(G.Nodes[n-1], edge.rNode)
			newEdge.data = limbLength
			break
		else :
			newNode = G.addNode()
			newEdge = G.addEdge(edge.lNode, newNode)
			newEdge.data = distance
			newEdge = G.addEdge(newNode, edge.lNode)
			newEdge.data = distance
			newEdge = G.addEdge(edge.rNode, newNode)
			newEdge.data = edge.data - distance
			newEdge = G.addEdge(newNode, edge.rNode)
			newEdge.data = edge.data - distance
			newEdge = G.addEdge(newNode, G.Nodes[n-1])
			newEdge.data = limbLength
			newEdge = G.addEdge(G.Nodes[n-1], newNode)
			newEdge.data = limbLength
			G.removeEdge(edge)
			G.removeEdgeBetween(edge.rNode, edge.lNode)
			break
				
	return G

# Construct an evolutionary tree from a distance/age matrix using Unweighted Pair Group Method with Arithmetic Means 
def UPGMA(matrix) :
	n = len(matrix)
	G = Graph(DIRECTED)
	for i in xrange(n) : node = G.addNode()
	clusterSize = []
	clusterSize = [1 for i in xrange(n)]
	numClusters = n

	while numClusters > 1 :
		print 'Matrix:', matrix
		minDist = sys.maxint
		for i in xrange(n) :
			if clusterSize[i] == 0 : continue
			for j in xrange(i+1, n) :
				if clusterSize[j] == 0 : continue
				dist = matrix[i][j]
				if dist < minDist :
					minDist = dist
					I, J = i, j

		print 'I,J: ', I, J
		node = G.addNode()
		node.data = minDist/2

		edge = G.addEdge(node, G.Nodes[I])
		edge.data = node.data - G.Nodes[I].data
		edge = G.addEdge(G.Nodes[I], node)
		edge.data = node.data - G.Nodes[I].data

		edge = G.addEdge(node, G.Nodes[J])
		edge.data = node.data - G.Nodes[J].data
		edge = G.addEdge(G.Nodes[J], node)
		edge.data = node.data - G.Nodes[J].data

		for k in xrange(n) :
			dist = matrix[I][k]*clusterSize[I] 
			dist += matrix[J][k]*clusterSize[J]
			dist = dist/(clusterSize[I] + clusterSize[J])
			matrix[k].append(dist)
		clusterSize.append(clusterSize[I] + clusterSize[J])
		clusterSize[I] = 0
		clusterSize[J] = 0
		newRow = []
		newRow = [matrix[k][n] for k in xrange(n)]
		newRow.append(0)
		matrix.append(newRow)
		n += 1
		numClusters -= 1

	return G

# Compute Neighbor-Joining Matrix for a distance matrix, D
def NeighborJoiningMatrix(D) :
	n = len(D)
	totalDistance = []
	totalDistance = [sum(list(D[i])) for i in xrange(n)]
	NJMat = []
	for i in xrange(n) :
		row = []
		for j in xrange(n) :
			if i == j : row.append(0)
			else : row.append((n-2)*D[i][j] - totalDistance[i] - totalDistance[j])
		NJMat.append(row)
	return NJMat

# Find the neighbors to be joined from a distance Matrix, D
def NeighborsToJoin(D, totalDistance) :
	n = len(D)

	minMetric = sys.maxint
	for i in xrange(n) :
		for j in xrange(n) :
			if i != j : 
				x = (n-2)*D[i][j] - totalDistance[i] - totalDistance[j]
				if x < minMetric :
					minMetric = x
					I, J = i, j
	return I, J

# Construct an Evolutionary Tree using Neighbor Joining Algorithm
def NeighborJoining(distanceMatrix, evolutionaryTree, nodeArray) :
	print distanceMatrix
	for x in nodeArray : print x.id,
	print

	n = len(distanceMatrix)
	if n == 2 :
		edge = evolutionaryTree.addEdge(nodeArray[0], nodeArray[1])
		edge.data = distanceMatrix[0][1]
		edge = evolutionaryTree.addEdge(nodeArray[1], nodeArray[0])
		edge.data = distanceMatrix[0][1]
		return evolutionaryTree

	totalDistance = []
	totalDistance = [sum(list(distanceMatrix[i])) for i in xrange(n)]
	I, J = NeighborsToJoin(distanceMatrix, totalDistance)
	node1 = nodeArray[I]
	node2 = nodeArray[J]

	delta = (totalDistance[I] - totalDistance[J])/float(n-2)
	limbLength1 = (distanceMatrix[I][J] + delta)/2.0
	limbLength2 = (distanceMatrix[I][J] - delta)/2.0

	p, q = min(I,J), max(I,J)
	nodeArray = nodeArray[:p] + nodeArray[p+1:q] + nodeArray[q+1:]
	node = evolutionaryTree.addNode()
	nodeArray.append(node)

	for k in xrange(n) :
		distanceMatrix[k].append((distanceMatrix[k][I] + distanceMatrix[k][J] - distanceMatrix[I][J])/2.0)
	newRow = []
	newRow = [distanceMatrix[k][n] for k in xrange(n)]
	newRow.append(0)
	distanceMatrix.append(newRow)

	for k in xrange(n+1) :
		distanceMatrix[k] = distanceMatrix[k][:p] + distanceMatrix[k][p+1:q]+distanceMatrix[k][q+1:]
	distanceMatrix = distanceMatrix[:p] + distanceMatrix[p+1:q] + distanceMatrix[q+1:]

	evolutionaryTree = NeighborJoining(distanceMatrix, evolutionaryTree, nodeArray)

	edge = evolutionaryTree.addEdge(node, node1)
	edge.data = limbLength1
	edge = evolutionaryTree.addEdge(node1, node)
	edge.data = limbLength1
		
	edge = evolutionaryTree.addEdge(node, node2)
	edge.data = limbLength2
	edge = evolutionaryTree.addEdge(node2, node)
	edge.data = limbLength2

	return evolutionaryTree

def Score(NodeData, id, sonId, daughterId) :
	for ch1 in 'ACGT' :
		minSonScore = sys.maxint
		minDaughterScore = sys.maxint
	
		for ch2 in 'ACGT' :
			if ch1 == ch2 : dist = 0
			else : dist = 1

			s1 = NodeData[sonId].minScore[ch2] + dist 
			if s1 < minSonScore :
				minSonScore = s1
				NodeData[id].minSon[ch1] = ch2

			s2 = NodeData[daughterId].minScore[ch2] + dist 
			if s2 < minDaughterScore :
				minDaughterScore = s2
				NodeData[id].minDaughter[ch1] = ch2

		NodeData[id].minScore[ch1] = minSonScore + minDaughterScore

def PropagateSolution(NodeData, node, ch) :
	# print "PropapagateSolution - node, ch :", node.id, ch

	son = node.outEdges[0].rNode
	if son.outEdges :
		ch1 = NodeData[node.id].minSon[ch]
		son.data += ch1
		PropagateSolution(NodeData, son, ch1)

	daughter = node.outEdges[1].rNode
	if daughter.outEdges :
		ch1 = NodeData[node.id].minDaughter[ch]
		daughter.data += ch1
		PropagateSolution(NodeData, daughter, ch1)

	return

# Given an evolutionary tree T with dna strings at the leaves, find dna strings for other nodes such that T is parsimonious
# That is, the total number of required mutations in T is minimum.

def SmallParsimony(T, dnasize) :
	N = len(T.Nodes)
	R = T.root()

	class data :
		def __init__(self) :
			self.minScore = {'A':0, 'C':0, 'G':0, 'T':0}
			self.minSon = {'A':'', 'C':'', 'G':'', 'T':''}
			self.minDaughter = {'A':'', 'C':'', 'G':'', 'T':''}

	for id in T.Nodes :
		node = T.Nodes[id]
		if not T.isLeaf(node) :	node.data = ''

	for id in T.Edges : T.Edges[id].data = 0

	total, totalScore = 0, 0

	for index in xrange(dnasize) :
		NodeData = [data() for i in range(N)]		

		numVisited = 0
		for id in T.Nodes :
			node = T.Nodes[id]

			if T.isLeaf(node) :	
				node.flag = 1
				numVisited += 1

				for ch in 'ACGT' :
					if node.data[index] == ch :
						NodeData[id].minScore[ch] = 0
					else :
						NodeData[id].minScore[ch] = sys.maxint
			else : node.flag = 0

		while numVisited < N :
			for id in T.Nodes :
				node = T.Nodes[id]
				if node.flag : continue

				visitedNeighbors = [n[0] for n in node.neighbors() if n[0].flag]
				if len(visitedNeighbors) < 2 : continue
				son = visitedNeighbors[0]
				daughter = visitedNeighbors[1]

				"""
				son = node.outEdges[0].rNode
				daughter = node.outEdges[1].rNode
				if not (son.flag and daughter.flag) : continue
				"""

				Score(NodeData, id, son.id, daughter.id)

				"""
				if id > 191 :
					print id, "->", son.id, son.data, NodeData[son.id].minScore
					print id, "->", daughter.id, daughter.data, NodeData[daughter.id].minScore
					print "      ", NodeData[id].minScore, NodeData[id].minSon, NodeData[id].minDaughter
					wait()
				"""

				node.flag = 1
				numVisited += 1

		"""
		for i in xrange(N) :
			print NodeData[i].minScore.values(), NodeData[i].minSon.values(), NodeData[i].minDaughter.values()
		print
		"""

		minScore = sys.maxint
		for key in NodeData[R.id].minScore :
			if NodeData[R.id].minScore[key] < minScore :
				minScore = NodeData[R.id].minScore[key] 
				ch = key

		total += minScore

		R.data += ch
		# print 'Root:', R.id
		PropagateSolution(NodeData, R, ch)

	# print "Total", total

	for id in T.Edges :
		e = T.Edges[id]
		e.data = HammingDist(e.lNode.data, e.rNode.data)
		totalScore += e.data

	# print "Total Score ", totalScore
	# print R.data
	return totalScore, T
#-------------------------------------------------------------------------------
# Parsimony string assignment to internal nodes of an UNROOTED tree
def SmallParsimonyUnrooted(G, dnasize) :

	numNodes = len(G.Nodes)
	root = G.Nodes[numNodes-1]

	# Pick arbitrarily an edge (We pick the edge incident to the first node) 
	e = G.Nodes[0].inEdges[0]
	son = e.lNode
	daughter = e.rNode

	e1 = G.addEdge(root, son)
	e2 = G.addEdge(root, daughter)
	G.removeEdge(e)

	totalScore, G = SmallParsimony(G, dnasize)

	# Remove the 2 edges we added and reinstate the original edge
	e = G.addEdge(son, daughter)
	e.data = e1.data + e2.data
	G.removeEdge(e1)
	G.removeEdge(e2)
	
	return totalScore, G
#-------------------------------------------------------------------------------
## Find the Nearest Neighbors for a tree
def NearestNeighbor(tree, nodeId1, nodeId2, configuration) :

	a = tree.Nodes[nodeId1]
	b = tree.Nodes[nodeId2]
	e = tree.getEdge(a, b)

	neighbors = [n[0] for n in a.neighbors() if n[0] is not b]
	w = neighbors[0]
	x = neighbors[1]

	neighbors = [n[0] for n in b.neighbors() if n[0] is not a]
	y = neighbors[0]
	z = neighbors[1]

	if configuration == 1 :
		tree.addEdge(a, y) 
		tree.addEdge(b, x) 
		tree.removeEdge(tree.getEdge(a, x))
		tree.removeEdge(tree.getEdge(b, y))
		return tree, a, b, x, y

	else :
		tree.addEdge(a, z) 
		tree.addEdge(b, x) 
		tree.removeEdge(tree.getEdge(a, x))
		tree.removeEdge(tree.getEdge(b, z))
		return tree, a, b, x, z

#-------------------------------------------------------------------------------
# Rewire to generate or revert nearest neighbor configurations
def Rewire(tree, a, b, x, y) :
	tree.addEdge(a, y) 
	tree.addEdge(b, x) 
	tree.removeEdge(tree.getEdge(a, x))
	tree.removeEdge(tree.getEdge(b, y))
	return tree 
#-------------------------------------------------------------------------------
# Explore different tree topologies and internal node assignments to find minimum parsimony.
# We explore the nearest neighbor configurations for every edge to greedily improve the parsimony score.
 
def LargeParsimony(G, numLeaves, dnasize, fpOut) :

	def printSolution() :
		print score
		for e in G.Edges.values() :
			print('%s->%s' % (e.lNode.id, e.rNode.id))
			print('%s->%s' % (e.rNode.id, e.lNode.id))
		print

		sys.stdout = fpOut

		print score
		for e in G.Edges.values() :
			print('%s->%s' % (e.lNode.id, e.rNode.id))
			print('%s->%s' % (e.rNode.id, e.lNode.id))
			print

	# score, G = SmallParsimonyUnrooted(G, dnasize)
 	# printSolution() 

	minScore = sys.float_info.max
 
	while True:
		solution = []
		for edge in G.Edges.values() :
			a = edge.lNode
			b = edge.rNode

			if a.id < numLeaves or b.id < numLeaves : continue
	
			# print 'a & b :', a.id, '->', b.id

			neighbors = [n[0] for n in a.neighbors() if n[0] is not b]
			w = neighbors[0]
			x = neighbors[1]

			neighbors = [n[0] for n in b.neighbors() if n[0] is not a]
			y = neighbors[0]
			z = neighbors[1]

			G = Rewire(G, a, b, x, y)	# Reconfigure topology to a nearest neighbor
			# print "Rewired Graph:"
			# G.print_adjacency(0)
			score, G = SmallParsimonyUnrooted(G, dnasize)
			G = Rewire(G, a, b, y, x)	# Revert to orig config


			if score < minScore :
				minScore = score
				solution = [a, b, x, y]

			G = Rewire(G, a, b, x, z)	# Reconfigure to next nearest neighbor
			# print "Rewired Graph:"
			# G.print_adjacency(0)
			score, G = SmallParsimonyUnrooted(G, dnasize)
			G = Rewire(G, a, b, z, x)	# Revert to orig config

			# print "Rewired Graph:"
			# G.print_adjacency(0)
			if score < minScore :
				minScore = score
				solution = [a, b, x, z]
	
		if solution == [] : break
		else :
			G = Rewire(G, solution[0], solution[1], solution[2], solution[3])
			printSolution()
				
#-------------------------------------------------------------------------------
#Prob 4.1 Generate the distance matrix for a graph
def p4_1(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numLeaves = int(fp.readline().rstrip())

	G = Graph(DIRECTED)
	while True :
		line = fp.readline()
		if not line : break

		edge, distance = line.rstrip().split(':')
		lNodeId, rNodeId = edge.split('->')

		lNode = G.addNode(int(lNodeId))
		rNode = G.addNode(int(rNodeId))
		edge = G.addEdge(lNode, rNode)
		edge.data = int(distance)
		# print lNode.id, '->', rNode.id, ':', edge.data

	matrix = GenerateDistanceMatrix(G, numLeaves)
		
	for i in xrange(numLeaves) :
		for j in xrange(numLeaves-1) :
			fpOut.write('%d ' % matrix[(i,j)])
		fpOut.write('%d\n' % matrix[(i,numLeaves-1)])
#-------------------------------------------------------------------------------
# Prob 4.2 Find the Limb Length for a Leaf from the Distance Matrix
def p4_2(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	numLeaves = int(fp.readline().rstrip())
	leafId = int(fp.readline().rstrip())

	matrix = []
	while True :
		line = fp.readline()
		if not line : break

		row = line.rstrip().split()
		row = [int(x) for x in row]
		matrix.append(row)
	# print matrix

	limbLength = LimbLength(matrix, numLeaves, leafId)
	
	fpOut.write('%d' % limbLength)
	print 'Limb Length: ', limbLength
#-------------------------------------------------------------------------------
# Prob 4.3 Construct a Phylogeny Tree from an additive distance matrix
def p4_3(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

        numLeaves = int(fp.readline().rstrip())

        matrix = []
        while True :
                line = fp.readline()
                if not line : break

                row = line.rstrip().split()
                row = [int(x) for x in row]
                matrix.append(row)

	G = Graph()
	for i in xrange(numLeaves) :
		G.addNode()

	G = AdditivePhylogeny(matrix, G, numLeaves)
	G.print_adjacency(1)
	for id in G.Edges :
		e = G.Edges[id]
		fpOut.write('%s->%s:%s\n' % (e.lNode.id, e.rNode.id, e.data))
#-------------------------------------------------------------------------------
# Prob 5.1 Construct an evolutionary tree from age matrix using UPGMA
def p5_1(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

        numLeaves = int(fp.readline().rstrip())

        matrix = []
        while True :
                line = fp.readline()
                if not line : break

                row = line.rstrip().split()
                row = [float(x) for x in row]
                matrix.append(row)

	G = UPGMA(matrix)

	G.print_adjacency(1)
	sys.stdout = fpOut
	G.print_adjacency(1)
#-------------------------------------------------------------------------------
# Prob 5.2 Construct an evolutionary tree from a distance matrix
def p5_2(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

        numLeaves = int(fp.readline().rstrip())

        matrix = []
        while True :
                line = fp.readline()
                if not line : break

                row = line.rstrip().split()
                row = [int(x) for x in row]
                matrix.append(row)

	G = Graph()
	for i in xrange(numLeaves) :
		G.addNode()

	nodeArray = []
	nodeArray = [G.Nodes[id] for id in xrange(numLeaves)]

	G = NeighborJoining(matrix, G, nodeArray)
	G.print_adjacency(1)
	for id in G.Edges :
		e = G.Edges[id]
		fpOut.write('%s->%s:%s\n' % (e.lNode.id, e.rNode.id, e.data))
##-------------------------------------------------------------------------------
# Prob 5.3 Construct parsimony tree for a given tree topology
def p5_3(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

        numLeaves = int(fp.readline().rstrip())
	G = Graph(DIRECTED)
	for i in xrange(2*numLeaves - 1) :
		G.addNode()

	count = 0
        while True :
                line = fp.readline()
                if not line : break

                left, right = line.rstrip().split('->')
		if count < numLeaves :
			G.Nodes[count].data = right
			dnasize = len(right)
			right = count
			count += 1
		else :
			right = int(right)

		lNode = G.Nodes[int(left)]
		rNode = G.Nodes[right]
		G.addEdge(lNode, rNode)
		print int(left), right, lNode.id, "->", rNode.id

	"""
	for id in G.Edges :
		e = G.Edges[id]
		print('%s->%s:%s\n' % (e.lNode.data, e.rNode.data, e.data))
	"""

	totalScore, G = SmallParsimony(G, dnasize)

	print totalScore
	for id in G.Edges :
		e = G.Edges[id]
		print('%s->%s:%s' % (e.lNode.data, e.rNode.data, e.data))
		print('%s->%s:%s' % (e.rNode.data, e.lNode.data, e.data))

	sys.stdout = fpOut

	print totalScore
	for id in G.Edges :
		e = G.Edges[id]
		print('%s->%s:%s' % (e.lNode.data, e.rNode.data, e.data))
		print('%s->%s:%s' % (e.rNode.data, e.lNode.data, e.data))
##-------------------------------------------------------------------------------
# Prob 5.4 Construct parsimony tree for a given UNROOTED tree topology
def p5_4(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

        numLeaves = int(fp.readline().rstrip())
	N = 2*numLeaves - 1 

	G = Graph(DIRECTED)
	for i in xrange(N) : G.addNode()

	count = 0
        while True :
                line = fp.readline()	# Ignore the edge in one of the directions
                line = fp.readline()
                if not line : break

                left, right = line.rstrip().split('->')
		if count < numLeaves :
			G.Nodes[count].data = right
			dnasize = len(right)
			right = count
			count += 1
		else :
			right = int(right)

		lNode = G.Nodes[int(left)]
		rNode = G.Nodes[right]
		G.addEdge(lNode, rNode)
		# print int(left), right, lNode.id, "->", rNode.id

	# Pick the last-added edge and split into 2 edges to create a root
	e = G.Edges[N-3]
	e1 = G.addEdge(G.Nodes[N-1], e.lNode)
	e2 = G.addEdge(G.Nodes[N-1], e.rNode)
	G.removeEdge(e)

	"""
	for id in G.Edges :
		e = G.Edges[id]
		print('%s->%s:%s\n' % (e.lNode.data, e.rNode.data, e.data))
	"""

	totalScore, G = SmallParsimony(G, dnasize)

	# Remove the 2 edges we added and reinstate the original edge
	e = G.addEdge(e.lNode, e.rNode)
	e.data = e1.data + e2.data
	G.removeEdge(e1)
	G.removeEdge(e2)

	print totalScore
	for id in G.Edges :
		e = G.Edges[id]
		print('%s->%s:%s' % (e.lNode.data, e.rNode.data, e.data))
		print('%s->%s:%s' % (e.rNode.data, e.lNode.data, e.data))

	sys.stdout = fpOut

	print totalScore
	for id in G.Edges :
		e = G.Edges[id]
		print('%s->%s:%s' % (e.lNode.data, e.rNode.data, e.data))
		print('%s->%s:%s' % (e.rNode.data, e.lNode.data, e.data))
##-------------------------------------------------------------------------------
# Prob 5.5 Generate Nearest Neighbors for a given tree topology
def p5_5(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

	G1 = Graph(UNDIRECTED)
	G2 = Graph(UNDIRECTED)

	a, b = fp.readline().rstrip().split()

	count = 0
        while True :
                line = fp.readline()	# Ignore the edge in one of the directions
                line = fp.readline()
                if not line : break
                left, right = line.rstrip().split('->')

		for G in [G1, G2] :
			lNode = G.addNode(int(left))
			rNode = G.addNode(int(right))
			G.addEdge(lNode, rNode)
			# print int(left), right, lNode.id, "->", rNode.id

		"""
		for id in G.Edges :
			e = G.Edges[id]
			print('%s->%s:%s\n' % (e.lNode.id, e.rNode.id, e.id))
		"""

	G1 = NearestNeighbor(G1, int(a), int(b), 1)
	G2 = NearestNeighbor(G2, int(a), int(b), 2)

	for id in G1.Edges :
		e = G1.Edges[id]
		print('%s->%s' % (e.lNode.id, e.rNode.id))
		print('%s->%s' % (e.rNode.id, e.lNode.id))
	print

	for id in G2.Edges :
		e = G2.Edges[id]
		print('%s->%s' % (e.lNode.id, e.rNode.id))
		print('%s->%s' % (e.rNode.id, e.lNode.id))

	sys.stdout = fpOut

	for id in G1.Edges :
		e = G1.Edges[id]
		print('%s->%s' % (e.lNode.id, e.rNode.id))
		print('%s->%s' % (e.rNode.id, e.lNode.id))
	print

	for id in G2.Edges :
		e = G2.Edges[id]
		print('%s->%s' % (e.lNode.id, e.rNode.id))
		print('%s->%s' % (e.rNode.id, e.lNode.id))

##-------------------------------------------------------------------------------
# Large Parsimony problem
def p5_6(infile, outfile) :
        fp = open(infile, 'r')
        fpOut = open(outfile, 'w')

	G = Graph(DIRECTED)

	numLeaves = int(fp.readline().rstrip())
	numNodes = 2*numLeaves - 1

	for i in xrange(numNodes) : G.addNode()

	count = 0
        while True :
                line = fp.readline()
                if not line : break

                left, right = line.rstrip().split('->')

		if not left.isdigit() : continue
		elif not right.isdigit() :
			G.Nodes[count].data = right
			dnasize = len(right)
			G.addEdge(G.Nodes[int(left)], G.Nodes[count])
			count += 1
		else :
			left, right = int(left), int(right)
			if left > right :
				G.addEdge(G.Nodes[left], G.Nodes[right])
	
	print "Input Graph:"
	G.print_adjacency(0)

	LargeParsimony(G, numLeaves, dnasize, fpOut)

#------------------------------------------------------------
"""
p4_1('./Regression/4.1.1', './Solutions/4.1.1')
p4_1('./Regression/4.1.2', './Solutions/4.1.2')
p4_1('./Regression/4.1.3', './Solutions/4.1.3')
p4_2('./Regression/4.2.1', './Solutions/4.2.1')
p4_2('./Regression/4.2.2', './Solutions/4.2.2')
p4_2('./Regression/4.2.3', './Solutions/4.2.3')
p4_3('./Regression/4.3.1', './Solutions/4.3.1')
p4_3('./Regression/4.3.2', './Solutions/4.3.2')
p4_3('./Regression/4.3.3', './Solutions/4.3.3')
p5_1('./Regression/5.1.1', './Solutions/5.1.1')
p5_1('./Regression/5.1.2', './Solutions/5.1.2')
p5_2('./Regression/5.2.1', './Solutions/5.2.1')
p5_2('./Regression/5.2.2', './Solutions/5.2.2')
p5_2('./Regression/5.2.3', './Solutions/5.2.3')
p5_3('./Regression/5.3.1', './Solutions/5.3.1')
p5_3('./Regression/5.3.2', './Solutions/5.3.2')
p5_3('./Regression/5.3.3', './Solutions/5.3.3')
p5_3('./Regression/5.3.4', './Solutions/5.3.4')
p5_4('./Regression/5.4.1', './Solutions/5.4.1')
p5_4('./Regression/5.4.2', './Solutions/5.4.2')
p5_4('./Regression/5.4.3', './Solutions/5.4.3')
p5_5('./Regression/5.5.1', './Solutions/5.5.1')
p5_5('./Regression/5.5.2', './Solutions/5.5.2')
p5_6('./Regression/5.6.1', './Solutions/5.6.1')
"""
p5_6('./Regression/5.6.2', './Solutions/5.6.2')
