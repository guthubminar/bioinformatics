#!/usr/bin/env python
# Author: Rajendran Panda, Dec. 2014

import copy 

UNDIRECTED = 0; DIRECTED = 1
NODES = 0; EDGES = 1; NODES_N_EDGES = 2; EDGES_N_NODES = 2

def sortByName(obj1, obj2) :
	return obj1.name > obj2.name

def adjacent(n1, n2) :
	l1 = len(n1)
	if len(n2) < l1 - 1:
		return False
	if (n1[1:l1] == n2[0:l1-1]) :
		return True
	return False

class Node :
	def __init__(self, nodeName = None) :
		self.name = nodeName
		self.Edges = []
		self.inEdges = []
		self.outEdges = []
		self.flag = 0
		self.data = 0

	def degree(self) :
		return len(self.Edges)
		
	def inDegree(self) :
		return len(self.inEdges)
		
	def outDegree(self) :
		return len(self.outEdges)
		
class Edge :
	def __init__(self, edgeName = None) :
		self.name = edgeName
		self.lNode = None
		self.rNode = None
		self.flag = 0
		self.weight = 0

class Graph :
        def __init__(self, type = DIRECTED) :
		self.type = type
		self.Nodes = []
		self.Edges = []
                self.numNodes = 0
		self.numEdges = 0

	def clearFlags(self, ofWhat = NODES_N_EDGES) :
		if ofWhat == EDGES or ofWhat == NODES_N_EDGES :
			for edge in self.Edges :
				edge.flag = 0

		if ofWhat == NODES or ofWhat == NODES_N_EDGES :
			for node in self.Nodes :
				node.flag = 0

	def findNode(self, nodeName) :
		for node in self.Nodes :
			if node.name == nodeName :
				return node
		return 0

	def findEdge(self, edgeName) :
		for edge in self.Edges :
			if edge.name == edgeName :
				return edge
		return 0

        def addEdge(self, lNode, rNode, edgeName = None) :
		edge = Edge(edgeName)
		edge.lNode = lNode
		edge.rNode = rNode

		if self.type == DIRECTED :	
			lNode.outEdges.append(edge)
			rNode.inEdges.append(edge)
		else :
			lNode.Edges.append(edge)
			rNode.Edges.append(edge)
		
		self.Edges.append(edge)
		self.numEdges += 1
		return edge

	def removeEdge(self, edge) :
		if self.type == DIRECTED :
			edge.lNode.outEdges.remove(edge)
			edge.rNode.inEdges.remove(edge)
			self.Edges.remove(edge)
		else :
			edge.lNode.Edges.remove(edge)
			edge.rNode.Edges.remove(edge)
			self.Edges.remove(edge)

	def removeEdgeBetween(self, lNode, rNode) :
		if self.type == DIRECTED :
			for edge in lNode.outEdges :
				if edge.rNode == rNode : break
			lNode.outEdges.remove(edge)
			rNode.inEdges.remove(edge)
			self.Edges.remove(edge)
		else :
			for edge in lNode.Edges :
				if edge.lNode == rNode or edge.rNode == rNode : 
					lNode.Edges.remove(edge)
					rNode.Edges.remove(edge)
			self.Edges.remove(edge)
		
        def print_adjacency(self) :
		if self.type == DIRECTED :
			ARROW = '->'
		else :
			ARROW = '-'


		for edge in self.Edges :
			print edge.lNode.name+ARROW+edge.rNode.name+':'+edge.name

	def addNode(self, nodeName) :
		node = self.findNode(nodeName)
		if not node :
			node = Node(nodeName)
			self.Nodes.append(node)
			self.numNodes += 1
		return node

	def colorEdges(self, color) :
		for edge in self.Edges :
			edge.data = color


	def extractCycle(self, node) :
		if node.flag == 1 : return []
		startNode = node

		cycle = []
		while True :
			for edge in node.Edges :
				if edge.flag == 0 :
					cycle.append(node)
					if edge.rNode == node :
						nextNode = edge.lNode
					else :
						nextNode = edge.rNode
					edge.flag = 1
					node.data -= 1	
					if node.data == 0 :
						node.flag = 1
					nextNode.data -= 1
					if nextNode.data == 0 :
						nextNode.flag = 1

					if nextNode == startNode :
						return cycle
					break
			node = nextNode

	def extractCycles(self) :
		self.clearFlags()
		for node in self.Nodes :
			node.data = len(node.Edges)

		Cycles = []
		while True :
			for node in self.Nodes :
				if node.flag == 0 : break

			if node.flag == 1 : return Cycles

			cycle = self.extractCycle(node)
			Cycles.append(cycle)
	
def ExtractCycle(node, numToVisit) :
	startNode = node
	Cycle = [node]
	
	# print 'Start : ', node.name, node.flag
	while True :   # Loop terminates after a cycle is extracted
		# Get the next node to the current (or starting) node
		for edge in node.outEdges :
			# print edge.rNode.name, edge.rNode.flag
			if edge.flag == 0 :
				node = edge.rNode
				edge.flag = 1
				# print 'Next : ', node.name
				break # from the for loop

		Cycle.append(node)
		node.flag -= 1 # because another inEdge is used

		if node.flag == 0:
			numToVisit -= 1  # one less node to visit

		if node == startNode :  # we completed a cycle
			break

	return Cycle, numToVisit


# Return an Eulerian Cycle of the graph
def EulerianCycle(graph) :
	Nodes = graph.Nodes
	numToVisit = len(Nodes)  # number of nodes to be visited

	# Clear flags of all nodes and edges
	for edge in graph.Edges :
		edge.flag = 0

	for node in Nodes :
		node.flag = len(node.inEdges) 
		if node.flag != len(node.outEdges) :
			print node.name, node.flag, len(node.outEdges)
			print 'Graph is not Eulerian'
			exit()

	# Extract the first subcycle 
	eCycle, numToVisit = ExtractCycle(Nodes[0], numToVisit)

	# Now extract all other subcycles iteratively
	while numToVisit > 0 :  # Loop terminates when all sub cycles are extracted
		cycleLength = len(eCycle)
		for i in range(cycleLength) :
			if eCycle[i].flag > 0 : break
		cycle, numToVisit = ExtractCycle(eCycle[i], numToVisit)

		eCycle = eCycle[:i] + cycle[:] + eCycle[i+1:]
	return eCycle


# Return an Eulerian Path of the graph
def EulerianPath(graph) :
        Nodes = graph.Nodes
        numToVisit = len(Nodes)  # number of nodes to be visited

	found = 0	
        for node in Nodes :
                node.flag = len(node.inEdges)
		numOutEdges = len(node.outEdges)
                if node.flag == numOutEdges + 1 :
			sink = node
			print 'Sink : ', sink.name
			found += 1
			if found == 2 : break
		elif node.flag == numOutEdges - 1 :
			src = node
			print 'Src : ', src.name
			found += 1
			if found == 2 : break
		elif node.flag != numOutEdges :
			print 'The graph has no Euler Path'
			exit() 

	# Add a virtual edge from sink to src
	edge = graph.addEdge(sink, src)
	src.flag += 1


	eCycle = EulerianCycle(graph)

	for j in range(len(eCycle)) :
		if ((eCycle[j] == sink) and (eCycle[j+1] == src)) : break
	ePath = eCycle[j+1:] + eCycle[1:j+1] 

        return ePath

def Assemble(str1, str2) :
	str2 = list(str2)
	m = len(str2)
	str1 = list(str1)
	str1.append(str2[m-1])
	n = len(str1)
	for i in range(m) :
		if str1[n-m+i] == '-' :
			str1[n-m+i] = str2[i]
		elif str1[n-m+i] != str2[i] and str2[i] != '-' :
			print 'Can not assemble ', str1, ' and ', str2
			return 0, "".join(str1)
	return 1, "".join(str1)
		
def ReadPairsPath2Sequence(path, k, d) :
	gap = "".join(['-' for i in range(d)])
	prefix = path[0].name[:k]
	suffix = path[0].name[k+1:]
	print prefix, gap, suffix
	seq = prefix + gap + suffix
	for node in path[1:] :
		prefix = node.name[:k]
		suffix = node.name[k+1:]
		text = prefix + gap + suffix
		# print 'Assembling: ', seq, text
		success, seq = Assemble(seq, text)
		if not success :
			return 0, None
	return 1, seq

def Contigs(G) :
	srcSink = []
	for node in G.Nodes :
		if len(node.inEdges) == 1 and len(node.outEdges) == 1 : continue
		# else
		node.flag = len(node.outEdges)
		srcSink.append(node)

	edgesToVisit = len(G.Edges)
	contigs = []
	while edgesToVisit > 0 :
		for node in srcSink :
			if node.flag == 0 : continue
			path = []
			path.append(node)
			node.flag -= 1
			while True :
				for edge in node.outEdges :
					if edge.flag == 1 : continue
					edge.flag = 1
					edgesToVisit -= 1
					node = edge.rNode
					path.append(node)
					break
				if node in srcSink : break
			contigs.append(path)
	return contigs

def mergeGraphs(P, Q) :
	R = copy.deepcopy(P)

	for edge in Q.Edges :
		lNode = R.findNode(edge.lNode.name)
		rNode = R.findNode(edge.rNode.name)
		newEdge = R.addEdge(lNode, rNode, edge.name)
		newEdge.data = edge.data
	return R
#---------------------------------------------------------------------------
def deBruijnGraph(fp) :
	g = Graph()

	k = int(fp.readline().rstrip())
	text = fp.readline().rstrip()
	ks = len(text)

	for i in range(ks-k+1) :
		lNode = g.addNode(text[i:i+k-1])
		rNode = g.addNode(text[i+1:i+k])
		edge = g.addEdge(lNode, rNode, text[i:i+k])
		""" print lNode.name, '->', edge.name, '->', rNode.name """

	g.print_adjacency()

def dBGfromPatterns(fp) :
	g = Graph()

	while True:
		text = fp.readline().rstrip()
		if not text : break

		k = len(text)
                lNode = g.addNode(text[:k-1])
                rNode = g.addNode(text[1:])
                edge = g.addEdge(lNode, rNode, text)

	# g.print_adjacency()
	return g

def dBGfromReadPairs(fp) :
        g = Graph()

	k, d = fp.readline().rstrip().split()
	k = int(k); d = int(d)

        while True:
                text = fp.readline().rstrip()
                if not text : break

                lNode = g.addNode(text[:k-1] + '/' + text[k+1:2*k])
                rNode = g.addNode(text[1:k]+ '/' + text[k+2:])
                edge = g.addEdge(lNode, rNode, text)

        # g.print_adjacency()
        return g, k, d

def findEulerianPathOrCycle(fp, isPath) :
	g = Graph()

	while True :
		line = fp.readline()
		if not line : break

		[lName, dummy, rNames] = line.rstrip().split()
		rNames = rNames.split(',')

		lNode = g.addNode(lName)
		for name in rNames :
			rNode = g.addNode(name)
			g.addEdge(lNode, rNode)

	if isPath :  
		eCycle = EulerianPath(g)
	else : 
		eCycle = EulerianCycle(g)

	print eCycle[0].name,
	for node in eCycle[1:] :
		print '->%s' % node.name, 
	print

def findUniversalString(k) :
	g = Graph()

	for i in xrange (2**k) :
		x = i
		num = ''
		for j in xrange(k):
			r = x % 2
			x = x / 2
			num = str(r) + num
		lNode = g.addNode(num[:k-1])
               	rNode = g.addNode(num[1:])
               	edge = g.addEdge(lNode, rNode, num)
		print 'Graph: ', lNode.name, rNode.name, edge.name

        g.print_adjacency()
        cycle = EulerianCycle(g)
	uniNum = cycle[0].name
	for node in cycle[1:2**k-k+2] :
		uniNum = uniNum + node.name[k-2]

	print uniNum


#----------------------------------------------------------------------------
