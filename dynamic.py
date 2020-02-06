#!/usr/bin/env python
import sys
sys.setrecursionlimit(100000)

def DPChange(money, Coins) :
	minCoins = {0:0}
	for m in range(1, money+1) :
		minCoins[m] = 10000
		for i in range(len(Coins)) :
			if m < Coins[i] : continue
			numCoins = minCoins[m - Coins[i]] + 1
			if numCoins < minCoins[m] :
				minCoins[m] = numCoins
	return minCoins[money]

def ManhattanTourist(n, m, Down, Right) :
	S = []
	S.append([0 for x in range(m+1)])
	
	for i in range(1, n+1) :
		S.append([0 for x in range(m+1)])
		S[i][0] = S[i-1][0] + Down[i-1][0]
	for j in range(1, m+1) :
		S[0][j] = S[0][j-1] + Right[0][j-1]
	# print S
	for i in range(1, n+1) :
		for j in range (1, m+1) :
			S[i][j] = max(S[i-1][j] + Down[i-1][j], S[i][j-1] + Right[i][j-1])
	# print S
	return S[n][m]

def LCSBacktrack(v, w) :
	South = 1
	East = 2
	SE = 3

	n = len(v)
	m = len(w)

	# Create n+1 by m+1 matrices, S and Backtrack
	S = []
	Backtrack = []
	for i in range(n+1) :
		S.append([0 for x in range(m+1)])
		Backtrack.append([0 for x in range(m+1)])

	for i in range (n+1) :
		S[i][0] = 0

	for j in range (m+1) :
		S[0][j] = 0

	for i in range(1,n+1) :
		for j in range(1,m+1) :
			S[i][j] = max(S[i-1][j], S[i][j-1])
			if v[i-1] == w[j-1] :
				S[i][j] = max(S[i][j], S[i-1][j-1] + 1)
			if S[i][j] == S[i-1][j] :
				Backtrack[i][j] = South
			elif S[i][j] == S[i][j-1] :
				Backtrack[i][j] = East
			else :	Backtrack[i][j] = SE
        return Backtrack

def OutputLCS(Backtrack, v, i, j) :
	# print 'i,j = ', i, j
	South = 1; East = 2; SE = 3

	if i==0 or j==0 : return ''
	if Backtrack[i][j] == South :
		return OutputLCS(Backtrack, v, i-1, j)
	elif Backtrack[i][j] == East :
		return OutputLCS(Backtrack, v, i, j-1)
        else :  
		s = OutputLCS(Backtrack, v, i-1, j-1)+v[i-1]
		return s

def MinimumEdit(v, w) :
        South = 1; East = 2; SE = 3

        n = len(v); m = len(w)

        # Create n+1 by m+1 matrices, S and Backtrack
        S = []; Backtrack = []
        for i in range(n+1) :
                S.append([0 for x in range(m+1)])
                Backtrack.append([0 for x in range(m+1)])

        for i in range (n+1) :
                S[i][0] = i 

        for j in range (m+1) :
                S[0][j] = j 

        for i in range(1,n+1) :
                for j in range(1,m+1) :
                        S[i][j] = min(S[i-1][j], S[i][j-1]) + 1
                        if v[i-1] == w[j-1] :
                                S[i][j] = min(S[i][j], S[i-1][j-1])
			else: 	S[i][j] = min(S[i][j], S[i-1][j-1] + 1)
                        if S[i][j] == S[i-1][j] + 1 :
                                Backtrack[i][j] = South
                        elif S[i][j] == S[i][j-1] + 1 :
                                Backtrack[i][j] = East
                        else :  Backtrack[i][j] = SE
        return S[n][m], Backtrack

def OutputLCS2(Backtrack, text1, text2, i, j) :
        Origin = 0; South = 1; East = 2; SE = 3

        # print i, j, Backtrack[i][j]
        if (i==0 and j==0) or Backtrack[i][j] == Origin :
                s1 = ''
                s2 = ''
        elif Backtrack[i][j] == South :
                s1, s2 = OutputLCS2(Backtrack, text1, text2, i-1, j)
                s1 = s1 + text1[i-1]
                s2 = s2 + '-'
        elif Backtrack[i][j] == East :
                s1, s2 = OutputLCS2(Backtrack, text1, text2, i, j-1)
                s1 = s1 + '-'
                s2 = s2 + text2[j-1]
        else :
                s1, s2 = OutputLCS2(Backtrack, text1, text2, i-1, j-1)
                s1 = s1 + text1[i-1]
                s2 = s2 + text2[j-1]
        return s1, s2

#--------------------------------------------------------------------------
def DPChange_main(fp) :
	money = int(fp.readline().rstrip())
	tmp = fp.readline().rstrip().split(',')
	Coins = [int(x) for x in tmp]
	print DPChange(money, Coins)

def ManhattanTourist_main(fp) : 
	Down = []
	Right = []
	n, m = fp.readline().rstrip().split()
	n, m = int(n), int(m)

	for row in range(1, n+1) :
		Down.append([int(item) for item in fp.readline().rstrip().split()])

	fp.readline()

	for row in range(0, n+1) :
		Right.append([int(item) for item in fp.readline().rstrip().split()])

	print ManhattanTourist(n,m,Down,Right)

def AlignSequences(fp) :
	v = fp.readline().rstrip()
	w = fp.readline().rstrip()
	Backtrack = LCSBacktrack(v,w)
	common = OutputLCS(Backtrack, v, len(v), len(w))
	print common

def MinimumEdit_main(fp) :
        v = fp.readline().rstrip()
        w = fp.readline().rstrip()
        cost, Backtrack = MinimumEdit(v,w)
        print cost
	s1, s2 = OutputLCS2(Backtrack, v, w, len(v), len(w))
	print s1
	print s2


def LongestPath_main(fp) :
	src = int(fp.readline().rstrip())
	sink = int(fp.readline().rstrip())

	
#--------------------------------------------------------------------------
fn = raw_input('fn? ')
fp = open(fn, 'r')

# AlignSequences(fp)
MinimumEdit_main(fp)
