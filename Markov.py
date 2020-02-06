#!/usr/bin/env python
# Rajndran Panda	May 2015

import sys

#-------------------------------------------------------------
def PrintTransitionMatrix(states, Matrix) :
	print '\t',
	for s in states :
		print s, '\t',
	print

	for i in range(len(states)) :
		print states[i], '\t',
		for j in range(len(states)) :
			x = Matrix[i][j]
			if x == 0.0 : x = '0'
			else : x = "%0.3f" % x
			print x, '\t',
		print
#-------------------------------------------------------------
def PrintEmissionMatrix(states, symbols, Matrix) :
	print '\t',
	for s in symbols :
		print s, '\t',
	print

	for i in range(len(states)) :
		print states[i], '\t',
		for j in range(len(symbols)) :
			x = Matrix[i][j]
			if x == 0.0 : x = '0'
			else : x = "%0.3f" % x
			print x, '\t',
		print
#-------------------------------------------------------------
# Find probability of a path,  given the transition probability matrix
def PathProbability(path, tranProbs) :
	size = len(path)

	pathProb = 0.5
	for i in range(size-1) :
		pathProb *= tranProbs[(path[i], path[i+1])]

	return pathProb
#-------------------------------------------------------------
# Find emission probability of a path,  given the emission probability matrix
def PathEmissionProbability(path, hiddenPath, emissionProbs) :
	size = len(path)

	pathProb = 1.0
	for i in range(size) :
		pathProb *= emissionProbs[(hiddenPath[i], path[i])]

	return pathProb
#-------------------------------------------------------------
# Find hidden path maximizing the emission of a given string
def HiddenPath(str, symbols, states, transProbs, emissionProbs) :

	M = len(states)
	N = len(str)

        symbol2index = { symbols[i]:i for i in range(len(symbols)) }

	# Create an M x N matrix to hold scores of nodes 
	scores = [ [0.0 for j in range(N)] for i in range(M) ]

	# Create an M x N matrix to hold pointers from each node to the optimal state (in previous column)
	optStates = [ [0 for j in range(N)] for i in range(M) ]

	# Initialize the scores for the first column
	for i in range(M) :
		scores[i][0] = emissionProbs[i][symbol2index[str[0]]] 
		# print scores[i][0],
	# print
	
	for j in range(1, N) :
		for i in range(M) :
			maxScore = 0.0
			for k in range(M) :
				score = scores[k][j-1] * transProbs[k][i]
				if score > maxScore :
					maxScore = score
					optStates[i][j] = k
			scores[i][j] = maxScore * emissionProbs[i][symbol2index[str[j]]]
		
	# print scores
	# print optStates

	maxScore = 0.0
	for i in range(M) :
		if scores[i][N-1] > maxScore :
			maxScore = scores[i][N-1] 
			I = i

	hiddenPath = states[I]

	for j in range(N-1,0,-1) :
		I = optStates[I][j]
		hiddenPath = states[I] + hiddenPath

	return hiddenPath

#-------------------------------------------------------------
# Find the likelihood of a given outcome
def OutcomeLikelihood(str, symbols, states, transProbs, emissionProbs) :

	M = len(states)
	N = len(str)
	symbol2index = { symbols[i]:i for i in range(len(symbols)) }

	# Create an M x N matrix to hold outcome likelihood of nodes, called forward 
	forward = [ [0.0 for j in range(N)] for i in range(M) ]

	# Initialize the scores for the first column
	for i in range(M) :
		forward[i][0] = emissionProbs[i][symbol2index[str[0]]] / M 
		# print forward[i][0],
	# print
	
	for j in range(1, N) :
		for i in range(M) :
			score = 0.0
			for k in range(M) :
				score += forward[k][j-1] * transProbs[k][i]
			forward[i][j] = score * emissionProbs[i][symbol2index[str[j]]]
		
	# print forward

	likelihood = 0.0
	for i in range(M) :
		likelihood += forward[i][N-1] 

	return likelihood

#-------------------------------------------------------------
# Path probability problem
def p6_6(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	str = fp.readline().strip()
	fp.readline()
	states = fp.readline().strip().split()
	numStates = len(states)
	fp.readline()
	toStates = fp.readline().strip().split()
	tranProbs = {}

	for i in range(numStates) :
		line = fp.readline().strip().split()
		for j in range(numStates) :
			tranProbs[(line[0], states[j])] = float(line[j+1])

	prob = PathProbability(str, tranProbs)
	print prob
	sys.stdout = fpOut
	print prob
#-------------------------------------------------------------
def p6_7(infile, outfile) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	str = fp.readline().strip()
	fp.readline()
	symbols = fp.readline().strip().split()
	numSymbols = len(symbols)
	fp.readline()

	hiddenPath = fp.readline().strip()
	fp.readline()
	states = fp.readline().strip().split()
	numStates = len(states)
	fp.readline()

	emittedSymbols = fp.readline().strip().split()
	emissionProbs = {}

	for i in range(numStates) :
		line = fp.readline().strip().split()
		for j in range(numSymbols) :
			emissionProbs[(line[0], emittedSymbols[j])] = float(line[j+1])

	print emissionProbs

	prob = PathEmissionProbability(str, hiddenPath, emissionProbs)
	print prob
	sys.stdout = fpOut
	print prob
#-------------------------------------------------------------
FIND_HIDDEN_PATH = 0
FIND_OUTCOME_LIKELIHOOD = 1
FIND_CONDITIONAL_PROBS = 2

def p6_8(infile, outfile, PROBLEM_TYPE) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	str = fp.readline().strip()
	fp.readline()

	symbols = fp.readline().strip().split()
	fp.readline()

	states = fp.readline().strip().split()
	fp.readline()

	fp.readline().strip().split()
	
	numSymbols = len(symbols)
	numStates = len(states)

	transProbs = []
	for i in range(numStates) :
		line = fp.readline().strip().split()
		row = [float(x) for x in line[1:]]
		transProbs.append(row)

	# print transProbs

	fp.readline()
	fp.readline().strip().split()
	emissionProbs = []

	for i in range(numStates) :
		line = fp.readline().strip().split()
		row = [float(x) for x in line[1:]]
		emissionProbs.append(row)

	# print emissionProbs

	if PROBLEM_TYPE == FIND_HIDDEN_PATH :
		hiddenPath = HiddenPath(str, symbols, states, transProbs, emissionProbs)

		# print hiddenPath
		sys.stdout = fpOut
		print hiddenPath

	if PROBLEM_TYPE == FIND_OUTCOME_LIKELIHOOD :
		likelihood = OutcomeLikelihood(str, symbols, states, transProbs, emissionProbs)

		# print likelihood
		sys.stdout = fpOut
		print likelihood

	if PROBLEM_TYPE == FIND_CONDITIONAL_PROBS :
		conditionalProbs = SoftDecoding(str, symbols, states, transProbs, emissionProbs)

		# print conditional probs
		sys.stdout = fpOut		
		for state in states :
			print state, '\t',
		print

		for i in range(len(str)) :
			for j in range(numStates) :
				print '%0.4f' % conditionalProbs[i][j], '\t',
			print
#-------------------------------------------------------------
# Construct a Hidden Markov Model for the given (multiple) Alignment
# 	Alignment is multiple texts aligned using indels. insertionThresh is used to classify a column as insertion.
#	HMM uses states: S, I0, M1, D1, I1, M2, D2, I2, ........., E
#  	where S is Start, Mi are normal, Di are Deletion, Ii are Insertion and E is End states.
#	Transition arcs exist as follows: 
#		S->I0, S->M1, S->D1
#		Ii->Mj, Ii->Ii, Ii->Dj, where j = i+1
#		Di->Mj, Di->Ij, Di->Dj, where j = i+1

def HMM(Alignment, symbols, insertionThresh, pseudocount, fpOut) :
	numTexts = len(Alignment); textLength = len(Alignment[0])
	mIndex2col = {} 	# index to column number dictionary for M columns

	index = 0
	for i in range(textLength) :
		count = 0 
		for text in Alignment :
			if (text[i] == '-') : count += 1
		if (float(count)/float(numTexts) < insertionThresh) : 
			index += 1; mIndex2col[index] = i
			colType.append('M'); colIndex.append(index)
		else :
			colType.append('I'); colIndex.append(index)

	maxIndex = colIndex[-1]
	Insertions = [[0 for i in range(maxIndex+1)] for j in range(numTexts)]

	for j in range(textLength) :
		if colType[j] == 'I' : 
			for i in range(numTexts) :
				if Alignment[i][j] != '-': Insertions[i][colIndex[j]] += 1

	# Form an array of states
	states = [('S', 0), ('I', 0)] 
	for i in range(1, maxIndex+1) : states += [('M', i), ('D', i), ('I', i)]
	states += [('E', maxIndex+1)] 

	# Create and initialize a transitions matrix 
	Transitions = {}
	for s1 in states:
		for s2 in states:
			Transitions[(s1, s2)] = 0

	# S2* transitions
	if colType[0] == 'I':
		for i in range(numTexts) :
			if Insertions[i][0] > 0: Transitions[('S0', 'I0')] += 1
			else :
				k = mIndex2col[1]
				if maxIndex >= k :
					if Alignment[i][k] != '-': Transitions[('S', 'M1')] += 1
					else : Transitions[('S0', 'D1')] += 1
	else : # colType[0] = 'M'
		for i in range(numTexts) :
			if Alignment[i][1] != '-': Transitions[('S0', 'M1')] += 1
			else : Transitions[('S0', 'D1')] += 1

	
	# I2I, I2M, I2D transitions
	for j in range(maxIndex + 1) :
		sI = ('I', j)
		for i in range(numTexts) :
			if Insertions[i][j] > 0 : Transitions[(sI,sI)] += Insertions[i][j] - 1

		if j+1 in mIndex2col : 
			k = mIndex2col[j+1]
			sM = ('M', j+1); sD = ('D', j+1)
			for i in range(numTexts) :
				if Insertions[i][j] > 0 :
					if Alignment[i][k] == '-' : Transitions[(sI, sD)] += 1
					else : Transitions[(sI, sM)] += 1

	# *2E transitions
	sI = ('I', maxIndex)
	sM = ('M', maxIndex)
	sD = ('D', maxIndex)
	sE = ('E', maxIndex+1)

	if colType[-1] == 'I' :
		for i in range(numTexts) :
			if Insertions[i][-1] > 0 : Transitions[(sI, sE)] += 1

	k = mIndex2col[maxIndex]	
		for i in range(numTexts) :
			if Alignment[i][-1] == '-' : Transitions[(sD, sE)] += 1
			else : Transitions[(sM, sE)] += 1


	# M2* transitions
	for i in range(numTexts) :
		if colType[j] == 'I' : continue

		if Alignment[i][j] == '-' :
			Transitions[(states[-4], states[-1])] += 1
		else :
			Transitions[(states[-3], states[-1])] += 1

	for j in range(textLength) :

		if colType[j] == 'I': continue

		if j == textLength - 1 :	# M2E and D2E transitions
			for i in range(numTexts) :
				if Alignment[i][j] != '-' :
					Transitions[(states[-4], states[-1])] += 1
				else :
					Transitions[(states[-3], states[-1])] += 1
		else : # M2* and D2* transitions
			jj = colIndex[j]
			nextM = index2mColNum[jj+1]
			for i in range(numTexts) :
				if Alignment[i][j] == '-' :
					if Alignment[i][nextM] == '-' : 
					Transitions[()] += 1
				else : Transitions[(states[0], states[3]) += 1
	else :
		for i in range(numTexts) :
			if Alignment[i][1] != '-': Transitions[(states[0], states[2]) += 1
			else : Transitions[(states[0], states[3]) += 1
	

		type = colType[j]
		index = colIndex[j]

		if colType[j] == 'M' :
			for i in range(numTexts) :
				if Alignment[i][j] != '-' :
					if Insertions[k] > 0 :
						transProbs[(fromM, toI) += 1
					if j+2 < textLength :
					if Insertions[colIndex[j+1]] > 0 :
						

	for i in range(-1, textLength-1) :
		print 'Column:', i; print '----------'

		m2m, m2i, m2d, m2e = 0, 0, 0, 0
		d2m, d2i, d2d, d2e = 0, 0, 0, 0
		
		if i == -1 :	# this is start node, which is taken to be in normal (M) state 
			for text in Alignment :
				for j in range(i+1, textLength) :
					if text[j] != '-' :
						if colType[j] == 'M' : m2m += 1
						else : 
							m2i += 1; break
					if colType[j] == 'M' : 
						m2d += 1; break

			total = m2i + m2m + m2d
			transProbs[('S', 'I0')] = float(m2i) / total 
			transProbs[('S', 'M1')] = float(m2m) / total 
			transProbs[('S', 'D1')] = float(m2d) / total

		else : # i = 0, 1, ....textLength-2
			if (i == 0) or (colType[i-1] == 'M') : i2m, i2i, i2d, i2e = 0, 0, 0, 0
			for text in Alignment :
				if text[i] == '-' :
					if colType[i] == 'I' : continue

					# Column type is 'M'
					end = 1
					for j in range(i+1, textLength) :
						if text[j] != '-' : 
							if colType[j] == 'M' : 
								d2m += 1
								print text[i], text[j], 'd2m:', d2m
							else : 
								d2i += 1
								print text[i], text[j], 'd2i:', d2i
							end = 0
							break
						if colType[j] == 'M': 
							d2d += 1
							print text[i], text[j], 'd2d:', d2d
							end = 0
							break
					if end :
						d2e += 1
						print text[i], 'd2e:', d2e
				else : 
					end = 1
					for j in range(i+1, textLength) :
						if text[j] != '-' : 
							if colType[i] == 'M' :
								if colType[j] == 'M' : 
									m2m += 1
									print text[i], text[j], 'm2m:', m2m
								else : 
									m2i += 1
									print text[i], text[j], 'm2i:', m2i
							else :
								if colType[j] == 'M' : 
									i2m += 1
									print text[i], text[j], 'i2m:', i2m
								else : 
									i2i += 1
									print text[i], text[j], 'i2i:', i2i
							end = 0
							break

						if colType[j] == 'M': 
							if colType[i] == 'M' :
								m2d += 1
								print text[i], text[j], 'm2d:', m2d
							else :
								i2d += 1
								print text[i], text[j], 'i2d:', i2d
							end = 0
							break

					if end :
						if colType[i] == 'M' :
							m2e += 1		
							print text[i], 'm2e:', m2e
						else :
							i2e += 1		
							print text[i], 'i2e:', i2e
			
			# print 'm2d:', m2d, 'm2i:', m2i, 'm2m:', m2m, 'm2e:', m2e
			# print 'd2d:', d2d, 'd2i:', d2i, 'd2m:', d2m, 'd2e:', d2e
			# print 'i2d:', i2d, 'i2i:', i2i, 'i2m:', i2m, 'i2e:', i2e

			if i == 0 and colType[0] == 'M'  :	# Set I0* when the first column is not I0
				transProbs[('I'+str(0), 'D'+str(1))] = 0
				transProbs[('I'+str(0), 'I'+str(0))] = 0
				transProbs[('I'+str(0), 'M'+str(1))] = 0
				
		
			k = suffix[i]

			if i == textLength-2 :
				if colType[i] == 'M' :
					if colType[i+1] == 'M' :
						m2e = m2m + d2m
						d2e = m2d + d2d

						if m2e > 0: transProbs[('M'+str(k+1), 'E')] = 1.0
						if d2e > 0: transProbs[('D'+str(k+1), 'E')] = 1.0
				

						total = d2d + d2m 
						if total > 0 :
							d2d = float(d2d) / total
							d2m = float(d2m) / total
						transProbs[('D'+str(k), 'D'+str(k+1))] = d2d
						transProbs[('D'+str(k), 'M'+str(k+1))] = d2m
				

						total = m2d + m2m 
						if total > 0.0 :
							m2d = float(m2d) / total
							m2m = float(m2m) / total
						transProbs[('M'+str(k), 'D'+str(k+1))] = m2d
						transProbs[('M'+str(k), 'M'+str(k+1))] = m2m

					else : # colType[i+1] = 'I'
						transProbs[('I'+str(k), 'E')] = 1.0
				
						total = d2i + d2e
						if total > 0 :
							d2i = float(d2i) / total
							d2e = float(d2e) / total
						transProbs[('D'+str(k), 'I'+str(k))] = d2i
						transProbs[('D'+str(k), 'E')] = d2e
						
						total = m2i + m2e
						if total > 0.0 :
							m2i = float(m2i) / total
							m2e = float(m2e) / total
						transProbs[('M'+str(k), 'I'+str(k))] = m2i
						transProbs[('M'+str(k), 'E')] = m2e

				if colType[i] == 'I' :
					if colType[i+1] == 'I' :
						for text in Alignment :
							if text[i+1] != '-' : i2e += 1

						total = i2i + i2e 
						if total > 0 :
							i2i = float(i2i) / total 
							i2e = float(i2e) / total 
						transProbs[('I'+str(k), 'I'+str(k))] = i2i
						transProbs[('I'+str(k), 'E')] = i2e
					else : # colType[i+1] = 'M' :
						for text in Alignment :
							if text[i+1] == '-' : d2e += 1
							else : m2e += 1
						if m2e > 0 : transProbs[('M'+str(k+1), 'E')] = 1.0
						if d2e > 0 : transProbs[('D'+str(k+1), 'E')] = 1.0

						total = i2d + i2m 
						if total > 0 :
							i2d = float(i2d) / total 
							i2m = float(i2m) / total 
						transProbs[('I'+str(k), 'D'+str(k+1))] = i2d
						transProbs[('I'+str(k), 'M'+str(k+1))] = i2m

			else :
				if colType[i] == 'M' :
					total = m2m + m2d + m2i + m2e
					if total > 0 :
						m2m = float(m2m) / total
						m2d = float(m2d) / total
						m2i = float(m2i) / total
						m2e = float(m2e) / total

						transProbs[('M'+str(k), 'M'+str(k+1))] = m2m
						transProbs[('M'+str(k), 'D'+str(k+1))] = m2d
						transProbs[('M'+str(k), 'I'+str(k))] = m2i
						transProbs[('M'+str(k), 'E')] = m2e

					total = d2m + d2d + d2i + d2e
					if total > 0 :
						d2m = float(d2m) / total
						d2d = float(d2d) / total
						d2i = float(d2i) / total
						d2e = float(d2e) / total

						transProbs[('D'+str(k), 'M'+str(k+1))] = d2m
						transProbs[('D'+str(k), 'D'+str(k+1))] = d2d
						transProbs[('D'+str(k), 'I'+str(k))] = d2i
						transProbs[('D'+str(k), 'E')] = d2e

				else : # colType = 'I'
					if colType[i+1] == 'M' :
						total = i2m + i2d + i2i
						i2m = float(i2m) / total 
						i2d = float(i2d) / total 
						i2i = float(i2i) / total 
						transProbs[('I'+str(k), 'M'+str(k+1))] = i2m
						transProbs[('I'+str(k), 'D'+str(k+1))] = i2d
						transProbs[('I'+str(k), 'I'+str(k))] = i2i

	# Adjust with pseudocount
	if pseudocount > 0.0 :
		for state in states :
			if state == 'E' : continue

			type, index = state[0], state[1:]
			if index == '' : index = 0
			else : index = int(index)

			sm = 'M' + str(index+1)
			sd = 'D' + str(index+1)
			si = 'I' + str(index)

			if index == maxSuffix : 
				if (state, si) in transProbs :
					transProbs[(state, si)] += pseudocount
				else :
					transProbs[(state, si)] = pseudocount
				if (state, 'E') in transProbs :
					transProbs[(state, 'E')] += pseudocount
				else :
					transProbs[(state, 'E')] = pseudocount

				total = transProbs[(state, si)] + transProbs[(state, 'E')] 

				transProbs[(state, si)] /= total
				transProbs[(state, 'E')] /= total
			else :
				if (state, sm) in transProbs :
					transProbs[(state, sm)] += pseudocount
				else : 
					transProbs[(state, sm)] = pseudocount

				if (state, sd) in transProbs :
					transProbs[(state, sd)] += pseudocount
				else : 
					transProbs[(state, sd)] += pseudocount

				if (state, si) in transProbs :
					transProbs[(state, si)] += pseudocount
				else :
					transProbs[(state, si)] = pseudocount

				total = transProbs[(state, sm)] + transProbs[(state, sd)] + transProbs[(state, si)] 

				transProbs[(state, sm)] /= total
				transProbs[(state, sd)] /= total
				transProbs[(state, si)] /= total

	sys.stdout = fpOut

	# Print the Transition Matrix

	print '\t',
	for s1 in states :
		print s1, '\t',
	print

	for s1 in states :
		print s1,
		for s2 in states :
			if (s1, s2) in transProbs:
				x = transProbs[(s1, s2)]
				if x == 0.0 : x = '0'
				else : x = "%0.3f" % transProbs[(s1, s2)]
				print '\t', x, 
			else :
				print '\t0',
		print

	print '--------'

	# Calculate Emission Matrix
	
	symbol2index = {}
	for i in range(len(symbols)) :
		symbol2index[symbols[i]] = i

	for i in range(textLength) :
		key = colType[i]+str(suffix[i])

		if key not in emissionProbs :
			emissionProbs[key] = [0 for j in range(len(symbols))]

		for text in Alignment :
			if (text[i] != '-') : emissionProbs[key][symbol2index[text[i]]] += 1

	for key in emissionProbs :
		total = sum(emissionProbs[key])
		emissionProbs[key] = [float(x)/total for x in emissionProbs[key]] 
		emissionProbs[key] = [x + pseudocount for x in emissionProbs[key]] 
		total = sum(emissionProbs[key])
		emissionProbs[key] = [float(x)/total for x in emissionProbs[key]] 

	if pseudocount > 0.0 :
		for i in range(maxSuffix+1) :
			key = 'I' + str(i)
			if key not in emissionProbs :
				emissionProbs[key] = [1.0/len(symbols) for j in range(len(symbols))]

	# Print Emission Matrix
	
	for ch in symbols :
		print '\t', ch,
	print

	for key in states :
		print key,
		if key in emissionProbs :
			for x in emissionProbs[key] : 
				if x == 0.0 : x = '0'
				else : x = "%0.3f" % x
				print '\t', x,
		else :
			for i in range(len(symbols)) : print '\t0',
		print
#-------------------------------------------------------------
def EstimateHMM(text, path, symbols, states) :
	numStates = len(states)
	numSymbols = len(symbols)
	pathLength = len(path)

	symbol2index = { symbols[i]:i for i in range(numSymbols) }
	state2index = { states[i]:i for i in range(numStates) }

	Transition = [[0 for i in range(numStates)] for j in range(numStates)]
	Emission = [[0 for i in range(numSymbols)] for j in range(numStates)]

	# Compute Transition Matrix
	for i in range(pathLength-1) :
		row = state2index[path[i]]
		col = state2index[path[i+1]]
		Transition[row][col] += 1

	for row in Transition :
		total = sum(row)
		for i in range(numStates) :
			if total > 0 : row[i] = float(row[i])/total
			else : row[i] = 1.0/numStates
	
	# Compute Emission Matrix
	for i in range(pathLength) :
		row = state2index[path[i]]
		col = symbol2index[text[i]]
		Emission[row][col] += 1

	for row in Emission :
		total = sum(row)
		for i in range(numSymbols) :
			if total > 0 : row[i] = float(row[i])/total
			else : row[i] = 1.0/numSymbols
		
	return Transition, Emission
#-------------------------------------------------------------
def ViterbiLearning(text, symbols, states, transProbs, emissionProbs, maxIterations) :
	while maxIterations :
		path = HiddenPath(text, symbols, states, transProbs, emissionProbs)
		transProbs, emissionProbs = EstimateHMM(text, path, symbols, states)
		maxIterations -= 1

	return transProbs, emissionProbs

#-------------------------------------------------------------
def SoftDecoding(str, symbols, states, transProbs, emissionProbs) :

	M = len(states)
	N = len(str)
	symbol2index = { symbols[i]:i for i in range(len(symbols)) }

	# Create M x N matrices to hold forward and backward scores of nodes 
	forward = [ [0.0 for j in range(N)] for i in range(M) ]
	backward = [ [0.0 for j in range(N)] for i in range(M) ]

	conditionalProbs = [ [0.0 for j in range(M)] for i in range(N) ]

	# Initialize the scores for the first and last column
	for i in range(M) :
		forward[i][0] = emissionProbs[i][symbol2index[str[0]]] 
		# backward[i][N-1] = emissionProbs[i][symbol2index[str[N-1]]] 
		backward[i][N-1] = 1.0
	
	for j in range(1, N) :
		for i in range(M) :
			forward[i][j] = 0.0
			for k in range(M) :
				forward[i][j] += forward[k][j-1] * transProbs[k][i] 
			forward[i][j] *= emissionProbs[i][symbol2index[str[j]]]
	
	forwardOfSink = sum([forward[i][N-1] for i in range(M)])
 	print 'forward(sink) =', forwardOfSink
	
	for j in range(N-2, -1, -1) :
		for i in range(M) :
			backward[i][j] = 0.0
			for k in range(M) :
				backward[i][j] += backward[k][j+1] * transProbs[i][k] * emissionProbs[k][symbol2index[str[j+1]]]

	for i in range(N) :
		for j in range(M) :
			conditionalProbs[i][j] = forward[j][i] * backward[j][i] / forwardOfSink
		"""
		total = sum(conditionalProbs[i])
		for j in range(M) :
			conditionalProbs[i][j] /= total
		"""

	return conditionalProbs
#-------------------------------------------------------------
PSEUDOCOUNT = 1

def p10_1(infile, outfile, pseudocount = 0) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	if pseudocount :
		threshold, epsilon = fp.readline().strip().split()
		threshold = float(threshold)
		epsilon = float(epsilon)
	else :
		threshold = float(fp.readline().strip())

	fp.readline()
	symbols = fp.readline().strip().split()
	fp.readline()

	Alignment = []
	while True:
		text = fp.readline()
		if not text : break
		Alignment.append(text.strip())

	if pseudocount :
		HMM(Alignment, symbols, threshold, epsilon, fpOut)
	else :
		HMM(Alignment, symbols, threshold, 0.0, fpOut)
#-------------------------------------------------------------
ESTIMATE_HMM = 0
VITERBI_LEARNING = 1

def p10_4(infile, outfile, problem) :
	fp = open(infile, 'r')
	fpOut = open(outfile, 'w')

	if problem == VITERBI_LEARNING :
		maxIter = int(fp.readline().strip())
		fp.readline()

	text = fp.readline().strip()
	fp.readline()
	symbols = fp.readline().strip().split()
	fp.readline()

	if problem == ESTIMATE_HMM :
		path = fp.readline().strip()
		fp.readline()

	states = fp.readline().strip().split()
	fp.readline()

	if problem == ESTIMATE_HMM :
		transProbs, emissionProbs = EstimateHMM(text, path, symbols, states)

	else : # problem = VITERBI_LEARNING
		numSymbols = len(symbols)
		numStates = len(states)

		transProbs = []
		fp.readline()
		for i in range(numStates) :
			line = fp.readline().strip().split()
			row = [float(x) for x in line[1:]]
			transProbs.append(row)

		fp.readline()
		fp.readline()
		emissionProbs = []

		for i in range(numStates) :
			line = fp.readline().strip().split()
			row = [float(x) for x in line[1:]]
			emissionProbs.append(row)


		transProbs, emissionProbs = ViterbiLearning(text, symbols, states, transProbs, emissionProbs, maxIter)

	sys.stdout = fpOut
	PrintTransitionMatrix(states, transProbs)
	print '--------'
	PrintEmissionMatrix(states, symbols, emissionProbs) 
#-------------------------------------------------------------
"""
p6_6('./Regression/6.6.3', './Solutions/6.6.3')
p6_7('./Regression/6.7.1', './Solutions/6.7.1')
p6_7('./Regression/6.7.2', './Solutions/6.7.2')

p6_8('./Regression/6.8.1', './Solutions/6.8.1', FIND_HIDDEN_PATH)
p6_8('./Regression/6.8.2', './Solutions/6.8.2', FIND_HIDDEN_PATH)
p6_8('./Regression/6.8.3', './Solutions/6.8.3', FIND_HIDDEN_PATH)

p6_8('./Regression/6.9.1', './Solutions/6.9.1', FIND_OUTCOME_LIKELIHOOD)
p6_8('./Regression/6.9.2', './Solutions/6.9.2', FIND_OUTCOME_LIKELIHOOD)
"""
p10_1('./Regression/10.1.1', './Solutions/10.1.1') 
"""
p10_1('./Regression/10.1.2', './Solutions/10.1.2') 
p10_1('./Regression/10.1.3', './Solutions/10.1.3') 
p10_1('./Regression/10.1.4', './Solutions/10.1.4') 
p10_1('./Regression/10.1.5', './Solutions/10.1.5') 

p10_1('./Regression/10.2.1', './Solutions/10.2.1', PSEUDOCOUNT) 
p10_1('./Regression/10.2.2', './Solutions/10.2.2', PSEUDOCOUNT) 
p10_1('./Regression/10.2.3', './Solutions/10.2.3', PSEUDOCOUNT) 

p10_4('./Regression/10.4.1', './Solutions/10.4.1', ESTIMATE_HMM) 
p10_4('./Regression/10.4.2', './Solutions/10.4.2', ESTIMATE_HMM) 
p10_4('./Regression/10.4.3', './Solutions/10.4.3', ESTIMATE_HMM) 
p10_4('./Regression/10.5.1', './Solutions/10.5.1', VITERBI_LEARNING) 
p10_4('./Regression/10.5.2', './Solutions/10.5.2', VITERBI_LEARNING) 
p10_4('./Regression/10.5.3', './Solutions/10.5.3', VITERBI_LEARNING) 

p6_8('./Regression/10.6.1', './Solutions/10.6.1', FIND_CONDITIONAL_PROBS)
p6_8('./Regression/10.6.2', './Solutions/10.6.2', FIND_CONDITIONAL_PROBS)
p6_8('./Regression/10.6.3', './Solutions/10.6.3', FIND_CONDITIONAL_PROBS)
"""
