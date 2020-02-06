#!/usr/bin/env python

def MultiSeqAlignment(text1, text2, text3) :
	l = len(text1); m = len(text2); n = len(text3)

	# Initialize a l+1 by m+1 by n+1 matrix
	# Each matrix entry has a Cost and Path edge
	S = []
	for i in range(l+1) :
		if i>0 : 
			iDir = '0'
		else:	iDir = '1'
	
		T = []
		for j in range(m+1) :
			if j>0 : 
				ijDir = iDir + '0'
			else :	ijDir = iDir + '1'

			U = [[0, '000'] for k in range(n+1)] 
			for k in range(n+1) :
				if k>0 : 
					U[k][1] = ijDir + '0'
				else :	U[k][1] = ijDir + '1'
			T.append(U)
		S.append(T)

	for i in range(1, l+1) :
		for j in range(1, m+1) :
			for k in range(1, n+1) :
				if(text1[i-1] == text2[j-1]) and (text2[j-1] == text3[k-1]) :
					aligned = 1
				else :	aligned = 0

				if S[i-1][j-1][k-1][0]+aligned > S[i-1][j-1][k][0] :
					max1 = S[i-1][j-1][k-1][0]+aligned 
					edge1 = '000'
				else :	
					max1 = S[i-1][j-1][k][0] 
					edge1 = '001'

				if S[i-1][j][k-1][0] > S[i-1][j][k][0] :
					max2 = S[i-1][j][k-1][0]
					edge2 = '010'
				else :
					max2 = S[i-1][j][k][0]
					edge2 = '011'

				if S[i][j-1][k-1][0] > S[i][j-1][k][0] :
					max3 = S[i][j-1][k-1][0]
					edge3 = '100'
				else :
					max3 = S[i][j-1][k][0]
					edge3 = '101'

				if max1 > max2 :
					max4 = max1
					edge4 = edge1
				else :
					max4 = max2
					edge4 = edge2

				if max3 > S[i][j][k-1][0] :
					max5 = max3
					edge5 = edge3
				else :
					max5 = S[i][j][k-1][0]
					edge5 = '110'

				if max4 > max5 :
					score = max4
					edge = edge4
				else :
					score = max5
					edge = edge5

				S[i][j][k] = [score, edge]

	return S

def Backtrack(text1, text2, text3, S, x, y, z) :
	print 'x, y, z = ', x, y, z
	if x == 0 and y == 0 and z == 0 :
		return '', '', ''

	edge = S[x][y][z][1]

	if edge == '000' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x-1, y-1, z-1)
		return s1+text1[x-1], s2+text2[y-1], s3+text3[z-1]

	if edge == '001' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x-1, y-1, z)
		return s1+text1[x-1], s2+text2[y-1], s3+'-'

	if edge == '010' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x-1, y, z-1)
		return s1+text1[x-1], s2+'-', s3+text3[z-1]

	if edge == '011' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x-1, y, z)
		return s1+text1[x-1], s2+'-', s3+'-'

	if edge == '100' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x, y-1, z-1)
		return s1+'-', s2+text2[y-1], s3+text3[z-1]
	
	if edge == '101' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x, y-1, z)
		return s1+'-', s2+text2[y-1], s3+'-'

	if edge == '110' :
		s1, s2, s3 = Backtrack(text1, text2, text3, S, x, y, z-1)
		return s1+'-', s2+'-', s3+text3[z-1]
#------------------------------------------------------
fn = raw_input('File ? ')
fp = open(fn, 'r')

text1 = fp.readline().rstrip()
text2 = fp.readline().rstrip()
text3 = fp.readline().rstrip()

l = len(text1); m = len(text2); n = len(text3)

S = MultiSeqAlignment(text1, text2, text3)
s1, s2, s3 = Backtrack(text1, text2, text3, S, l, m, n)
print S[l][m][n][0]
print s1
print s2
print s3

