#!/usr/bin/env python 
# Author: Rajendran Panda, Mar. 7, 2015

from graph import *
# from multiprocessing import Queue, Process

VERBOSE = 1

# Data structure to store multiple strings
# This is a tree structure derived from DIRECTED Graph

class Trie(Graph) :
	def __init__(self) :
		Graph.__init__(self)
		self.root = self.addNode('0')
		self.nodeNum = 1

	# 'child' here is a successor of 'parent' through the edge with name 'name'
	def findChild(self, parent, name) :
		for edge in parent.outEdges :
			if edge.name == name :
				return edge.rNode
		return None

	# Function to store a string on the Trie
	def addString(self, string) :
		parent = self.root
		l = len(string)

		for i in range(l) :
			child = self.findChild(parent, string[i])
			if not child : break
			parent = child

		for j in range(i,l) :
			child = self.addNode(str(self.nodeNum))
			self.nodeNum += 1
			self.addEdge(parent, child, string[j])
			parent = child

	# Function that finds the locations on the 'text' where 
	# the strings stored on the Trie match
	def match(self, text) :
		l = len(text)
		for i in range(0, l) :
			parent = self.root
			for j in range(i, l+1) :
				if not parent.outEdges : 
					print i,
					break
				else :
					symbol = text[j]
					child = self.findChild(parent, symbol)
					if child : parent = child
					else : break

	# When the Trie stores all the suffixes of a text, it is called a Suffix tree.
	# Using a Suffix tree, this function finds the longest pattern that occurs more than once in the text.
	def longestRepeat(self) :
		self.root.data = 0
		unvisited = [self.root]
		longest = None
		maxLength = 0

		# The code below does a BFS of the graph, visiting only nodes with outDegree > 0
		# During the BFS, the length of strings from root to a node is stored as node data
		# We keep track of the node providing the longest string from the root
		while unvisited :
			node = unvisited.pop()
			for edge in node.outEdges :
				child = edge.rNode
				if (child.outDegree() > 0) :
					child.data = node.data + len(edge.name)
					unvisited.append(child)
					if child.data > maxLength :
						maxLength = child.data
						longest = child

		# We now back-trace from the 'longest' node to the 'root'
		# to obtain the longest pattern that repeats more than once
		print maxLength, longest.name	
		currNode = longest
		pattern = ''
		while currNode != self.root :
			pattern = currNode.inEdges[0].name + pattern
			currNode = currNode.inEdges[0].lNode
			print pattern

		print pattern

	
	# Function to find the longest shared substring between 'text' and the string of the suffix tree
	def sharedSubString(self, text, verbose=0) :
		t = text
		parent = self.root
		while t :
			deadEnd = 1
			if verbose == 1 : print t
			for edge in parent.outEdges :
				l = commonPrefixLength(t, edge.name)
				if l == 0 : continue
				if l == len(edge.name) :
					deadEnd = 0
					t = t[l:]
					parent = edge.rNode
					break	# for loop
				if l < len(edge.name) :
					deadEnd = 1
					t = t[l:]
					break	# for loop
					
			if deadEnd : break	# while loop
		return text[:(len(text) - len(t))]


# Function to provide the length of the common prefix between 2 strings	
def commonPrefixLength(string1, string2) :
	l1 = len(string1)
	l2 = len(string2)

	if l1 > l2 : l = l2
	else : l = l1

	count = 0
	for i in xrange(l) :
		if string1[i] == string2[i] : count += 1
		else : break

	return count

# Function to build a Suffix tree
# A suffix tree is a Trie where in all the suffixes of a 'text' are stored
def SuffixTree(text) :
	T = Trie()
	length = len(text)

	for i in xrange(length) :
		currNode = T.root
		subText = text[i:]
		
		while subText :
			deadEnd = 1
			for edge in currNode.outEdges :
				l = commonPrefixLength(subText, edge.name)

				if l == 0 : continue

				if l == len(edge.name) :
					deadEnd = 0
					subText = subText[l:]
					currNode = edge.rNode
					break  # for loop

				# else
				deadEnd = 1
				newNode = T.addNode(str(T.nodeNum))
				T.nodeNum += 1
				T.addEdge(currNode, newNode, edge.name[:l])
				# print 'Adding1 : ', currNode.name, '-->', newNode.name, ':', edge.name[:l]
				T.addEdge(newNode, edge.rNode, edge.name[l:])
				# print 'Adding2 : ', newNode.name, '-->', edge.rNode.name, ':', edge.name[l:]
				currNode = newNode
				subText = subText[l:]
				T.removeEdge(edge)
				# print 'Removing :', edge.lNode.name, '-->', edge.rNode.name, ':', edge.name
				break  # for loop

			if deadEnd and subText :
				newNode = T.addNode(str(T.nodeNum))
				T.nodeNum += 1
				T.addEdge(currNode, newNode, subText)
				# print 'Adding3 : ', currNode.name, '-->', newNode.name, ':', subText 
				break  # while loop
	return T

# Find the longest shared substring between 2 strings
# Each string is given by its suffix tree
def LongestSharedSubString(st1, st2) :
	# Keep st1 as the tree with fewer nodes.
	if len(st1.Nodes) > len(st2.Nodes) :
		st = st2
		st2 = st1
		st1 = st


	print 'st1: ', len(st1.Nodes)
	print 'st2: ', len(st2.Nodes)

	# For a node N in st1, let 'text' be the string from root to N.
	# Search st2 for the longest common prefix, call it the 'sharedText', between 'text' and any suffix string in st2
	# We do this search, starting from the root of st1 in a BFS manner.
	# If 'sharedText' is a proper subset of 'text', then the subtree rooted at N need not be processed.

	LSSS = ''        # longest shared sub string
	maxLength = 0
	
	st1.root.data = ''
	unvisited = [st1.root]

	while unvisited :
		parent = unvisited.pop()
		print 'Process : ', parent.name

		for edge in parent.outEdges :
			child = edge.rNode
			text = parent.data + edge.name
			sharedText = st2.sharedSubString(text)

			if len(sharedText) == len(text) :
				child.data = text
				unvisited.append(child)
			
			if len(sharedText) > maxLength :
				maxLength = len(sharedText)
				LSSS = sharedText

			parent.data = ''     # hope the memory for parent.data is freed by garbage collection

	print 'Longest Shared Substring is :'
	print LSSS
	print 'Length :', maxLength

# Find the shortest non-shared substring between 2 strings
# Each string is given by its suffix tree
def ShortestNonSharedSubString(st1, st2) :
	print 'st1: ', len(st1.Nodes)
	print 'st2: ', len(st2.Nodes)

	# For a node N in st1, let 'text' be the string from root to N.
	# Search st2 for the longest common prefix, call it the 'sharedText', between 'text' and any suffix string in st2
	# We do this search, starting from the root of st1 in a BFS manner.
	# If 'sharedText' is a proper subset of 'text', then the subtree rooted at N need not be processed.

	SNSS = ''        # shortest non-shared substring
	minLength = 1000000
	
	st1.root.data = ''
	unvisited = [st1.root]

	while unvisited :
		parent = unvisited.pop(0)
		if len(parent.data) >= minLength -1 :
			# print 'Not processing : ', parent.name, parent.data
			parent.data = ''
			continue

		# print 'Processing children of : ', parent.name, parent.data
		for edge in parent.outEdges :
			child = edge.rNode
			text = parent.data + edge.name
			sharedText = st2.sharedSubString(text)
			# print 'Processed: ', child.name, text
			# print 'Shared : ', sharedText
			# print

			if len(sharedText) == len(text) :
				if len(sharedText) < minLength - 1 :
					child.data = text
					unvisited.append(child)
			else :	
				l = len(sharedText) + 1
				if l < minLength :
					minLength = l
					SNSS = text[:l]
					# print sharedText
					# print SNSS
					# print

		parent.data = ''     # hope the memory for parent.data is freed by garbage collection

	print 'Shortest Non-Shared Substring is :'
	print SNSS
	print 'Length :', minLength

#------------Below are functions each executing one of the home work problems-------
# Prob 1: Given a bunch of strings, build a Trie
def ch7_1(fp) :
	T = Trie()

	while True :
		line = fp.readline()
		if not line : break

		string = line.rstrip()
		T.addString(string)

	for edge in T.Edges :
		print edge.lNode.name + '-->' + edge.rNode.name + ':' + edge.name
#----------------------------------------------------------------
# Prob 2: Build a Trie with the given set of strings and find the locations in the text where they match the text
def ch7_2(fp) :
	T = Trie()

	text = fp.readline().rstrip()

	while True :
		line = fp.readline()
		if not line : break

		string = line.rstrip()
		T.addString(string)

	T.match(text)
#----------------------------------------------------------------
# Prob 3: Given a text, build a suffix tree with all the suffixes of text
def ch7_3(fp) :
	
	text = fp.readline().rstrip()

	T = SuffixTree(text)

	for edge in T.Edges :
		print edge.name
#----------------------------------------------------------------
# Prob 4: Build a suffix tree of the given text and use it to find out
# the longest subsequence that repeats more than once
def ch7_4(fp) :
	text = fp.readline().rstrip()

	T = SuffixTree(text)
	for edge in T.Edges :
		print edge.lNode.name + '-->' + edge.rNode.name + ':' + edge.name

	T.longestRepeat()
#----------------------------------------------------------------
# Prob 5: Find the longest shared substring between 2 strings
def ch7_5(fp) :
	s1 = fp.readline().rstrip()
	s2 = fp.readline().rstrip()

	st1 = SuffixTree(s1)
	st2 = SuffixTree(s2)

	LongestSharedSubString(st1, st2)
#----------------------------------------------------------------
# Prob 6: Find the shortest non-shared substring between 2 strings
def ch7_6(fp) :
#	results = Queue()

	s1 = fp.readline().rstrip()
	s2 = fp.readline().rstrip()

#	p1 = Process(target = SuffixTree, args = (s1,results))	
#	p2 = Process(target = SuffixTree, args = (s2,results))

#	p1.start()
#	p2.start()

#	p1.join()
#	p2.join()

#	st1 = results.get()
#	st2 = results.get()

	st1 = SuffixTree(s1)
	st2 = SuffixTree(s2)
	ShortestNonSharedSubString(st1, st2)
#----------------------------------------------------------------

fn = raw_input("Input? :")
fp = open(fn, 'r')

#ch7_1(fp)
#ch7_2(fp)
#ch7_3(fp)
#ch7_4(fp)
#ch7_5(fp)
ch7_6(fp)
