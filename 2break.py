#!/usr/bin/env python
# Author: Rajendran Panda, Dec. 31, 2014

from graph import * # from my graph.py implementation

RED = 0; BLUE = 1

def Chromosome2Graph(chromosome, G) :

# A chromosome is specified as a sequence of synteny blocks (signed integers) where the sign
# represents the block orientation. 
# eg. of a chromosome: (+4 -3 -1 +2)

# Given a chromosome and a graph G (possibly an empty graph), this routine adds edges to G 
# corresponding to the connections between the synteny blocks. The chromosome is treated
# as cyclic, so an edge between the last block to the first block is also added. 
# In the graph, each synteny block N is represented by 2 nodes with labels 2N (head) and 
# 2N-1 (tail).

	numBlocks = len(chromosome)
	for blk in chromosome :
		G.addNode(abs(blk)*2)
		G.addNode(abs(blk)*2 - 1)

	for i in range(numBlocks) :
		j = i + 1
		if j == numBlocks : j = 0

		thisBlk = chromosome[i]
		nextBlk = chromosome[j]

		if thisBlk > 0 : n1 = G.findNode(abs(thisBlk)*2)
		else : n1 = G.findNode(abs(thisBlk)*2 - 1)

		if nextBlk > 0 : n2 = G.findNode(abs(nextBlk)*2 - 1)
		else : n2 = G.findNode(abs(nextBlk)*2) 
	
		G.addEdge(n1, n2, str(n1.name)+'-'+str(n2.name))

	# G.print_adjacency()

	return G

def Graph2Chromosomes(G) :
	G.clearFlags(NODES_N_EDGES)

	nodesToVisit = len(G.Nodes)
	chromo = ''

	while nodesToVisit > 0 :
		node = G.findNode(1)

		if node.flag :
			for node in G.Nodes :
				if node.flag == 0 : break

		startNode = node

		chromo = chromo + '('

		while True :
			blk = (node.name + 1)/2
			if node.name%2 :
				chromo = chromo + '+' + str(blk)
				exitNode = G.findNode(node.name + 1)
			else :
				chromo = chromo + '-' + str(blk)
				exitNode = G.findNode(node.name - 1)

			node.flag = 1; exitNode.flag = 1
			# print 'Visited :', node.name, exitNode.name
			nodesToVisit -= 2

			if exitNode.Edges[0].lNode == exitNode :
				node = exitNode.Edges[0].rNode
			else :
				node = exitNode.Edges[0].lNode

			if node == startNode :
				chromo = chromo + ')'
				break    # from while
			else :
				chromo = chromo + ' '

	print chromo

#------------------------------------------------------------------------------
fn = raw_input('Input file : ')
fp = open(fn, 'r')

P = Graph(UNDIRECTED)
line = fp.readline().rstrip()
length = len(line)
line = line[1:length-1].split(')(')
for chromo in line :
	chromo = chromo.split(' ')
	chromo = [int(x) for x in chromo]
	P = Chromosome2Graph(chromo, P)
Graph2Chromosomes(P)
P.colorEdges(RED)

G = Graph(UNDIRECTED)
line = fp.readline().rstrip()
length = len(line)
line = line[1:length-1].split(')(')
for chromo in line :
	chromo = chromo.split(' ')
	chromo = [int(x) for x in chromo]
	G = Chromosome2Graph(chromo, G)
Graph2Chromosomes(G)
G.colorEdges(BLUE)

H = mergeGraphs(P, G)

"""
for node in H.Nodes :
	print node.name,
print
for edge in H.Edges :
	print edge.name,
print

"""

Cycles = H.extractCycles()
print 'Num blocks: ', len(G.Nodes)/2
print 'Num cycles: ', len(Cycles)
print '2-break distance: ', len(G.Nodes)/2 - len(Cycles)

H.clearFlags(NODES_N_EDGES)
while True :
	# find a blue edge that is not visited
	for edge in H.Edges :
		if (edge.flag == 0) & (edge.data == BLUE) : break

	if (edge.flag != 0) & (edge.data != BLUE) : break # from while loop

	# print '*****', edge.name

	for edge1 in edge.lNode.Edges :
			if edge1.data == RED : break
	for edge2 in edge.rNode.Edges :
			if edge2.data == RED : break
	if edge1 == edge2 :
		edge.flag = 1
		continue	# the while loop

	i = edge.lNode; j = edge.rNode

	if edge1.lNode == i :
		k = edge1.rNode
	else :
		k = edge1.lNode

	if edge2.lNode == j :
		l = edge2.rNode
	else :
		l = edge2.lNode

	# print 'ijkl :', i.name, j.name, k.name, l.name
	H.removeEdgeBetween(i, k)
	H.removeEdgeBetween(j, l)
	e = H.addEdge(i, j)
	e.data = RED
	e = H.addEdge(k, l)
	e.data = RED

	i = P.findNode(i.name)
	j = P.findNode(j.name)
	k = P.findNode(k.name)
	l = P.findNode(l.name)

	# print 'ijkl :', i.name, j.name, k.name, l.name

	P.removeEdgeBetween(i, k)
	P.removeEdgeBetween(j, l)
	e = P.addEdge(i, j)
	e.data = RED
	e = P.addEdge(k, l)
	e.data = RED

	edge.flag = 1
	Graph2Chromosomes(P)
	
Cycles = H.extractCycles()
print 'Num blocks: ', len(G.Nodes)/2
print 'Num cycles: ', len(Cycles)
print '2-break distance: ', len(G.Nodes)/2 - len(Cycles)

""" for cycle in Cycles :
	cycle = [node.name for node in cycle]
	print cycle """
